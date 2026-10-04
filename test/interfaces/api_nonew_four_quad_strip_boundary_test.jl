if !isdefined(@__MODULE__,:QuadTriNoNewFourQuadStripCertificates)
    include("../geometry/quadtri_nonew_four_quad_strip_certificates.jl")
end
if !isdefined(@__MODULE__,:QuadTriNoNewTriangleCertificates)
    include("../geometry/quadtri_nonew_triangle_certificates.jl")
end

module NoNewFourQuadStripBoundaryTests

using Test,Tessella,SHA,TOML
using Tessella.Elements: msh_dimension
using ..QuadTriNoNewFourQuadStripCertificates
using ..QuadTriNoNewTriangleCertificates
const API=Tessella.API
const Cert=QuadTriNoNewFourQuadStripCertificates
const PROFILES=(:one,:three,:graded)

function with_model(f,action)
    API.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"four_quad_strip_boundary.geo");write(path,f.source)
            opened=API.open_geo!(path;mesh_dim=0);out=Int.(opened.lists["sweep"])
            return action((;source=f.surface,top=out[1],volume=out[2],laterals=out[3:end]))
        end
    finally
        API.finalize()
    end
end

function snapshot()
    model=API.CURRENT[]
    return (;model,value=repr(model),params=deepcopy(model.curve_params),
        cache=API.LAST_MESH[],class=API.LAST_MESH_CLASS[],
        overlay=API.LAST_MESH_HIGH_ORDER[],edges=API.LAST_MESH_EDGES[],faces=API.LAST_MESH_FACES[],
        payload=(API.mesh.get_nodes(),API.mesh.get_elements()),
        history=(API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[]),
        options=copy(API.OPTIONS),callback=model.meshing.size_callback)
end

function unchanged(before)
    @test API.CURRENT[]===before.model && repr(before.model)==before.value
    @test before.model.curve_params==before.params
    @test API.LAST_MESH[]===before.cache && API.LAST_MESH_CLASS[]===before.class
    @test API.LAST_MESH_HIGH_ORDER[]===before.overlay
    @test API.LAST_MESH_EDGES[]===before.edges && API.LAST_MESH_FACES[]===before.faces
    @test (API.mesh.get_nodes(),API.mesh.get_elements())==before.payload
    @test (API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[])==before.history
    @test API.OPTIONS==before.options && before.model.meshing.size_callback===before.callback
end

function actual_carriers(carriers)
    model=API.CURRENT[]
    source=Tessella.Model.mesh_model_surface(model,carriers.source)
    volume=Tessella.Model.mesh_model_volume(model,carriers.volume)
    projected=Tessella.Model.model_to_mixed(model,volume,3,carriers.volume)
    tags,xyz,_=API.mesh.get_nodes(3,carriers.volume,true)
    points=reshape(xyz,3,:)
    positions=Dict(Tuple(p)=>tag for (tag,p) in zip(tags,eachcol(points)))
    @test length(positions)==length(tags)
    remap=[positions[Tuple(projected.coords[:,i])] for i in axes(projected.coords,2)]
    edges=Dict{Tuple,Tuple{Int,Int}}();quads=Dict{Tuple,Tuple{Int,Int}}()
    for (index,block) in enumerate(projected.blocks)
        dim=msh_dimension(block.msh);dim in (1,2) || continue
        owners=projected.elementary_entities[index]
        for (column,cell) in enumerate(eachcol(block.nodes))
            carrier=(dim,Int(owners[column]))
            if block.msh==1
                edges[Cert.key(remap[v] for v in cell)]=carrier
            else
                for i in eachindex(cell)
                    support=Cert.key((remap[cell[i]],remap[cell[mod1(i+1,length(cell))]]))
                    haskey(edges,support) && edges[support][1]<dim && continue
                    @test !haskey(edges,support) || edges[support]==carrier
                    edges[support]=carrier
                end
                block.msh==3 && (quads[Cert.key(remap[v] for v in cell)]=carrier)
            end
        end
    end
    return (;edges,quads,volume=carriers.volume,source)
end

function support_nodes(volume)
    types,tags,connectivity=API.mesh.get_elements(3,volume)
    result=Dict{UInt64,Tuple}()
    for (msh,ids,nodes) in zip(types,tags,connectivity)
        for edge in eachcol(reshape(API.mesh.get_element_edge_nodes(msh,volume,false),3,:))
            support=Cert.key(edge[1:2]);node=edge[3]
            @test !haskey(result,node) || result[node]==support
            result[node]=support
        end
        if msh in (12,13,14)
            for face in eachcol(reshape(API.mesh.get_element_face_nodes(msh,4,volume,false),9,:))
                support=Cert.key(face[1:4]);node=face[9]
                @test !haskey(result,node) || result[node]==support
                result[node]=support
            end
        end
        if msh==12
            for cell in eachcol(reshape(nodes,27,:))
                result[cell[27]]=Cert.key(cell[1:8])
            end
        end
    end
    return result
end

function check_parameters(dim,entity,boundary,count)
    tags,xyz,uv=API.mesh.get_nodes(dim,entity,boundary,true)
    @test length(tags)==count && allunique(tags)
    @test length(xyz)==3count && length(uv)==dim*count
    @test maximum(abs.(API.model.get_value(dim,entity,uv).-xyz))<=2e-11
    for node in tags
        p,parameters,d,owner=API.mesh.get_node(node)
        @test length(parameters)==d
        @test maximum(abs.(API.model.get_value(d,owner,parameters).-p))<=2e-11
    end
    return Set(tags)
end

function check_quadratic(carriers,f,expected)
    n=length(f.levels)-1
    volume_nodes=API.mesh.get_nodes(3,carriers.volume,true)[1]
    @test length(volume_nodes)==54n+(f.laterals ? 59 : 27)
    counts=zeros(Int,4)
    for node in volume_nodes
        p,uv,dim,owner=API.mesh.get_node(node)
        counts[dim+1]+=1
        @test dim in (1,2) ? length(uv)==dim : isempty(uv)
        dim in (1,2) && (@test maximum(abs.(API.model.get_value(dim,owner,uv).-p))<=2e-11)
    end
    @test counts==[8,8n+28,32n-2,14n+(f.laterals ? 25 : -7)]
    types,tags,_=API.mesh.get_elements(3,carriers.volume)
    families=Dict(Int(t)=>length(ids) for (t,ids) in zip(types,tags))
    # Free choices are constrained by the complete joined relation, rather
    # than a pointer-dependent family signature. Every actual elevated cell
    # and node is checked below; recombined terminal fans have fixed counts.
    if f.laterals
        correct=Dict(12=>4(n-1),11=>8,14=>20)
        filter!(pair->last(pair)>0,correct)
        @test families==correct
    end
    supports=support_nodes(carriers.volume)
    @test length(supports)==44n+(f.laterals ? 45 : 17)
    @test length(volume_nodes)-length(supports)==10(n+1)+(f.laterals ? 4 : 0)
    signature=Dict{Tuple,Tuple}()
    for (node,support) in supports
        p,uv,dim,owner=API.mesh.get_node(node)
        carrier=length(support)==2 ? get(expected.edges,support,(3,carriers.volume)) :
            length(support)==4 ? get(expected.quads,support,(3,carriers.volume)) : (3,carriers.volume)
        @test (dim,owner)==carrier
        corners=[API.mesh.get_node(primary)[1] for primary in support]
        mean=[sum(c[d] for c in corners)/length(corners) for d in 1:3]
        @test maximum(abs.(p.-mean))<=2e-11
        signature[Tuple(sort!(Tuple.(corners)))]=(Tuple(p),carrier)
    end
    for surface in carriers.laterals
        curve=abs(API.CURRENT[].meshing.extrude_sources[(2,surface)][2])
        long=curve in (f.curve_tags[i] for i in f.long_pair)
        own=check_parameters(2,surface,false,long ? 14n-7 : 2n-1)
        closure=check_parameters(2,surface,true,long ? 18n+9 : 6n+3)
        @test own⊆closure
    end
    for surface in (carriers.source,carriers.top)
        own=check_parameters(2,surface,false,7)
        closure=check_parameters(2,surface,true,27)
        @test own⊆closure
    end
    for curve in f.curve_tags
        long=curve in (f.curve_tags[i] for i in f.long_pair)
        check_parameters(1,curve,false,long ? 7 : 1)
        check_parameters(1,curve,true,long ? 9 : 3)
    end
    source_and_top=Set(Int.(f.curve_tags))
    union!(source_and_top,(abs(curve) for (_,curve) in
        API.model.get_boundary([(2,carriers.top)],false,false,false)))
    lateral_curves=Set(abs(curve) for surface in carriers.laterals for (_,curve) in
        API.model.get_boundary([(2,surface)],false,false,false))
    vertical_curves=setdiff(lateral_curves,source_and_top)
    @test length(vertical_curves)==4
    for curve in vertical_curves
        # Only this region's four native corner columns participate here;
        # another independent region can have coincident curve geometry.
        check_parameters(1,curve,false,2n-1)
        check_parameters(1,curve,true,2n+1)
    end
    for msh in types
        q,_=API.mesh.get_integration_points(msh,"Gauss3")
        jac,determinants,xyz=API.mesh.get_jacobians(msh,q,carriers.volume)
        @test all(isfinite,jac) && all(>(0),determinants) && all(isfinite,xyz)
        @test length(jac)==9length(determinants) && length(xyz)==3length(determinants)
    end
    return signature
end

artifact_name(profile,laterals,direction,order)=
    "four_quad_strip_$(profile)_R$(laterals)_D$(direction)_P$(order)"

function artifact_record(name,order)
    base=QuadTriNoNewTriangleCertificates.api_record(name,API)
    buffer=IOBuffer();println(buffer,"four_quad_strip_support_query_v1");println(buffer,base.sha)
    class=API.LAST_MESH_CLASS[]
    for (label,catalog) in (("edge",class.edge_entities),("tri",class.face_entities),("quad",class.quad_entities))
        println(buffer,label);write(buffer,UInt64(length(catalog)))
        for (support,owner) in sort!(collect(catalog);by=first)
            write(buffer,UInt64(length(support)))
            for node in support;write(buffer,UInt64(node));end
            write(buffer,Int64(owner[1]),Int64(owner[2]))
        end
    end
    for (_,surface) in sort(API.model.get_entities(2)),boundary in (false,true)
        tags,xyz,uv=API.mesh.get_nodes(2,surface,boundary,true)
        write(buffer,UInt64(surface),UInt8(boundary),UInt64(length(tags)))
        for node in tags;write(buffer,UInt64(node));end
        for array in (xyz,uv)
            write(buffer,UInt64(length(array)))
            for value in array;write(buffer,reinterpret(UInt64,Float64(value)));end
        end
    end
    return merge(base,(;kind="api_p$(order)_supports",sha=bytes2hex(sha256(take!(buffer)))))
end

function artifact_records()
    records=NamedTuple[]
    for profile in PROFILES,laterals in (false,true),direction in (-1,1)
        f=Cert.fixture("artifact";layers=profile,laterals,height=direction,pins=(1,2,3,4))
        with_model(f,carriers->begin
            linear=API.mesh.generate(3)
            expected=actual_carriers(carriers)
            @test Cert.certify(linear,f,expected.source).total==1
            push!(records,artifact_record(artifact_name(profile,laterals,direction,1),1))
            API.mesh.set_order(2)
            check_quadratic(carriers,f,expected)
            push!(records,artifact_record(artifact_name(profile,laterals,direction,2),2))
        end)
    end
    return records
end

function artifact_pins()
    path=joinpath(@__DIR__,"..","artifacts","quadtri_nonew_four_quad_strip_crc.txt")
    records=Dict{String,NamedTuple}()
    for line in readlines(path)
        isempty(line) || startswith(line,"#") || begin
            fields=split(line,'\t');length(fields)==6 || error("four-Quad strip artifact row width")
            name,kind,n_nodes,n_cells,families,sha=String.(fields)
            haskey(records,name) && error("quad strip duplicate artifact name")
            records[name]=(;name,kind,n_nodes=parse(Int,n_nodes),n_cells=parse(Int,n_cells),families,sha)
        end
    end
    length(records)==24 || error("four-Quad strip artifact must contain24 products")
    return records
end

function saved_oracle_records()
    path=joinpath(@__DIR__,"..","artifacts","quadtri_nonew_four_quad_strip_oracle.toml")
    artifact=TOML.parsefile(path)
    @test artifact["gmsh_version"]=="4.15.2" && artifact["fixture_count"]==24
    return artifact["fixtures"]
end

function saved_oracles(records)
    originals=filter(f->occursin(r"^next_four_quad_strip_L[13]_",f["name"]),records)
    @test length(originals)==8
    return Dict((f["intervals"]==1 ? "L1" : "L3",f["recombine_laterals"],
                Int(sign(f["translation_height"])))=>f for f in originals)
end

# The fixture is reconstructed from the captured input and its raw metadata,
# rather than inferring geometry from a newly generated volume or template.
function oracle_fixture(saved)
    raw=saved["input_geo"]
    declarations=collect(eachmatch(r"Point\((\d+)\)=\{([^}]+)\};",raw))
    @test length(declarations)==4
    point_tags=Tuple(parse(Int,m.captures[1]) for m in declarations)
    corners=Tuple(Tuple(parse.(Float64,split(m.captures[2],','))[1:3]) for m in declarations)
    axis=Int(saved["variant"]["axis"])+1
    axes=axis==1 ? (2,3,1) : axis==2 ? (1,3,2) : (1,2,3)
    plane=axis==1 ? :YZ : axis==2 ? :XZ : :XY
    xy=Tuple(Tuple(Cert.Q(corner[d]) for d in axes[1:2]) for corner in corners)
    winding=Int(sign(Cert.orient(xy[1],xy[2],xy[3])))
    curves=Tuple(Int.(saved["source_curve_tags"]))
    long_pair=Tuple(i for i in 1:4 if length(saved["source_chains"][i]["nodes"])==5)
    @test length(long_pair)==2
    base=Cert.fixture(saved["name"])
    return merge(base,(;source=replace(raw,"Extrude{"=>"sweep[]=Extrude{";count=1)*";\n",
        corners,point_tags,curve_tags=curves,axes,plane,winding,long_pair,
        surface=Int(saved["source_surface"]),laterals=Bool(saved["recombine_laterals"]),
        height=Float64(saved["translation_height"]),
        levels=Float64.(saved["logical"]["requested_normalized_levels"])))
end

function check_saved_source_chains(saved)
    for chain in saved["source_chains"]
        curve=Int(chain["curve"])
        native=API.CURRENT[].curve_params[curve]
        actual=Float64.(chain["stored_parameters"])
        @test length(native)==length(actual)+2
        @test first(native)==0. && last(native)==1.
        isempty(actual) || (@test maximum(abs.(native[2:end-1].-actual))<=2e-11)
        tags,xyz,uv=API.mesh.get_nodes(1,curve,false,true)
        @test length(tags)==length(actual) && length(uv)==length(actual)
        isempty(actual) || (@test maximum(abs.(sort(uv).-sort(actual)))<=2e-11)
    end
end

function check_saved_primary(saved,volume)
    tags=API.mesh.get_nodes(3,volume,true)[1]
    original=saved["p1"]["nodes"]
    @test length(tags)==length(original)
    matched=Set{Int}()
    for node in tags
        p,_,dim,owner=API.mesh.get_node(node)
        hits=findall(q->maximum(abs.(p.-q["coordinates"]))<=2e-11,original)
        @test length(hits)==1
        length(hits)==1 || continue
        index=only(hits);push!(matched,index)
        @test (dim,owner)==Tuple(original[index]["owner"])
    end
    @test length(matched)==length(original)
end

function check_saved_boundary_counts(saved)
    @test length(API.mesh.get_nodes()[1])==length(saved["p2"]["nodes"])
    # Stored face-center UV can be empty upstream; native computed UV is
    # verified by reevaluation rather than copied from that raw provenance.
    for entity in saved["p2"]["entities"]
        dim,tag=entity["dim"],entity["tag"]
        dim in (1,2) || continue
        @test length(API.mesh.get_nodes(dim,tag,false)[1])==length(entity["owned"]["tags"])
        @test length(API.mesh.get_nodes(dim,tag,true)[1])==length(entity["closure"]["tags"])
    end
end

function node_geometry(volume=nothing)
    tags=volume===nothing ? API.mesh.get_nodes()[1] : API.mesh.get_nodes(3,volume,true)[1]
    Set((Tuple(p),dim,owner,Tuple(uv)) for node in tags
        for (p,uv,dim,owner) in (API.mesh.get_node(node),))
end

function cell_geometry(volume)
    types,tags,connectivity=API.mesh.get_elements(3,volume)
    Set((Int(msh),Tuple(sort!([Tuple(API.mesh.get_node(v)[1]) for v in cell])))
        for (msh,ids,nodes) in zip(types,tags,connectivity)
        for cell in eachcol(reshape(nodes,Tessella.Elements.msh_spec(msh).nnodes,:)))
end

function lower_primary_tags()
    primary=Set{UInt64}()
    types,tags,connectivity=API.mesh.get_elements()
    for (msh,nodes) in zip(types,connectivity)
        msh_dimension(msh)<=1 || continue
        width=Tessella.Elements.msh_spec(msh).nnodes
        count=msh==15 ? 1 : 2
        for cell in eachcol(reshape(nodes,width,:)),node in cell[1:count]
            push!(primary,node)
        end
    end
    return primary
end

function check_lower_regeneration(entity,old,current,order,old_primary,new_primary)
    oldtags,oldxyz,oldparams=old
    tags,xyz,params=current
    oldpoints=reshape(oldxyz,3,:);points=reshape(xyz,3,:)
    olduv=reshape(oldparams,entity[1],length(oldtags));uv=reshape(params,entity[1],length(tags))
    @test length(tags)==length(oldtags) && allunique(tags)
    @test Set(Tuple(p) for p in eachcol(uv))==Set(Tuple(p) for p in eachcol(olduv))
    signature(a,b,ids,selected)=Set((Tuple(a[:,i]),Tuple(b[:,i])) for i in eachindex(ids) if ids[i] in selected)
    @test signature(points,uv,tags,new_primary)==signature(oldpoints,olduv,oldtags,old_primary)
    if order==2
        # Gmsh Generator.cpp runs SetOrder1 before MeshOnlyEmpty and creates
        # requested P2 supports on CAD afterward. Primary sampling is retained;
        # old linear-average P2 bits need not equal new CAD evaluation bits.
        # A pinned edited-midpoint control independently proves this rebuild.
        @test xyz==API.model.get_value(entity[1],entity[2],params)
    else
        @test signature(points,uv,tags,Set(tags))==signature(oldpoints,olduv,oldtags,Set(oldtags))
    end
end

# MSH2 and MSH4 may group the same cells differently. Primary node identity,
# complete typed connectivity, physical tag and elementary owner are retained.
function file_cells(mesh;owners=true)
    return sort!([(block.msh,Tuple(cell),block.tags[column],
        !owners ? Int32(0) : mesh.entity_data===nothing ? mesh.elementary_entities[index][column] :
            mesh.entity_data.block_entities[index][column])
        for (index,block) in enumerate(mesh.blocks)
        for (column,cell) in enumerate(eachcol(block.nodes))])
end

function check_files(mesh,directory,stem;names=mesh.physical_names)
    classified=mesh.entity_data!==nothing || mesh.elementary_entities!==nothing
    for version in (2.2,4.1),binary in (false,true)
        path=joinpath(directory,"$(stem)_$(version)_$(binary).msh")
        Tessella.Elements.write_mixed_msh(path,mesh;version,binary)
        back=Tessella.Elements.read_mixed_msh(path)
        @test back.physical_names==names
        @test back.coords==mesh.coords
        @test file_cells(back;owners=classified)==file_cells(mesh;owners=classified)
        if version==4.1 && mesh.entity_data!==nothing
            @test back.entity_data.node_entities==mesh.entity_data.node_entities
            @test back.entity_data.node_parametric==mesh.entity_data.node_parametric
        end
    end
end

function run_tests()
@testset "Four-Quad strip actual P2 carriers and public lifecycle" begin
    captured=saved_oracle_records()
    @testset "Twelve P1/P2 native golden products" begin
        pins=artifact_pins()
        oracle=saved_oracles(captured)
        for profile in PROFILES,laterals in (false,true),direction in (-1,1)
            f=Cert.fixture("golden";layers=profile,laterals,height=direction,pins=(1,2,3,4))
            with_model(f,carriers->begin
                linear=API.mesh.generate(3)
                @test length(API.mesh.get_nodes()[1])==10length(f.levels)+(f.laterals ? 4 : 0)
                source=Tessella.Model.mesh_model_surface(API.CURRENT[],carriers.source)
                certificate=Cert.certify(linear,f,source)
                @test certificate.total==1
                owned=zeros(Int,4)
                for node in API.mesh.get_nodes(3,carriers.volume,true)[1]
                    owned[API.mesh.get_node(node)[3]+1]+=1
                end
                n=length(f.levels)-1
                @test owned==[8,4n+8,6n-6,f.laterals ? 4 : 0]
                # Higher-dimensional generation retains the published volume
                # cell cache while lower carriers are exposed through nodes.
                @test isempty(API.mesh.get_elements(2)[1])
                saved=get(oracle,(profile===:one ? "L1" : profile===:three ? "L3" : "graded",laterals,direction),nothing)
                saved===nothing || check_saved_primary(saved,carriers.volume)
                @test artifact_record(artifact_name(profile,laterals,direction,1),1)==pins[artifact_name(profile,laterals,direction,1)]
                expected=actual_carriers(carriers)
                API.mesh.set_order(2)
                @test artifact_record(artifact_name(profile,laterals,direction,2),2)==pins[artifact_name(profile,laterals,direction,2)]
                saved===nothing || check_saved_boundary_counts(saved)
                initial=check_quadratic(carriers,f,expected)
                before=snapshot();API.mesh.set_order(2)
                @test (API.mesh.get_nodes(),API.mesh.get_elements())==before.payload
                API.mesh.set_order(1)
                @test length(API.mesh.get_nodes()[1])==10length(f.levels)+(f.laterals ? 4 : 0)
                API.mesh.set_order(2)
                @test check_quadratic(carriers,f,expected)==initial
            end)
        end
    end

    @testset "Sixteen independent primary variants preserve actual source samples and carriers" begin
        variants=filter(record->startswith(record["name"],"next_four_quad_strip_v2_"),captured)
        @test length(variants)==16
        for saved in variants
            f=oracle_fixture(saved)
            with_model(f,carriers->begin
                linear=API.mesh.generate(3)
                expected=actual_carriers(carriers)
                @test Cert.certify(linear,f,expected.source).total>0
                check_saved_primary(saved,carriers.volume)
                check_saved_source_chains(saved)
                API.mesh.set_order(2)
                check_saved_boundary_counts(saved)
                check_quadratic(carriers,f,expected)
            end)
        end
    end

    @testset "Opposite winding, axes, rounded straight-CAD sampling and sparse entity tags" begin
        for plane in (:XY,:YZ,:ZX),laterals in (false,true)
            f=Cert.fixture("owners";plane,laterals,height=-1.,layers=:graded,shape=:rounded,
                           winding=-1,pins=:opposite,curve_reverse=5,tags=:sparse,
                           progression=4.,long_pair=plane===:YZ ? (2,4) : (1,3))
            with_model(f,carriers->begin
                linear=API.mesh.generate(3);expected=actual_carriers(carriers)
                @test Cert.certify(linear,f,expected.source).total>0
                API.mesh.set_order(2);check_quadratic(carriers,f,expected)
            end)
        end
    end

    @testset "Surface generation and immediate quadratic volume generation" begin
        for laterals in (false,true)
            f=Cert.fixture("shell";laterals,layers=:three)
            n=length(f.levels)-1
            with_model(f,carriers->begin
                API.mesh.generate(2)
                @test length(API.mesh.get_nodes()[1])==10(n+1)
                types,tags,_=API.mesh.get_elements(2)
                @test sum(length,tags)==(laterals ? 10n+12 : 20n+12)
                @test isempty(API.mesh.get_elements(3)[1])
                API.mesh.set_order(2)
                @test length(API.mesh.get_nodes()[1])==40n+34
                owned=zeros(Int,4)
                for node in API.mesh.get_nodes()[1]
                    owned[API.mesh.get_node(node)[3]+1]+=1
                end
                @test owned==[8,8n+28,32n-2,0]
                for surface in (carriers.source,carriers.top)
                    check_parameters(2,surface,false,7)
                    check_parameters(2,surface,true,27)
                end
                for surface in carriers.laterals
                    curve=abs(API.CURRENT[].meshing.extrude_sources[(2,surface)][2])
                    long=curve in (f.curve_tags[i] for i in f.long_pair)
                    check_parameters(2,surface,false,long ? 14n-7 : 2n-1)
                    check_parameters(2,surface,true,long ? 18n+9 : 6n+3)
                end
            end)
            with_model(f,carriers->begin
                API.option("Mesh.ElementOrder",2)
                API.mesh.generate(3)
                expected=actual_carriers(carriers)
                check_quadratic(carriers,f,expected)
                name=artifact_name(:three,laterals,1,2)
                @test artifact_record(name,2)==artifact_pins()[name]
            end)
        end
    end

    @testset "Actual identity and support geometry survive lifecycle operations" begin
        for laterals in (false,true),direction in (-1,1)
            f=Cert.fixture("lifecycle";laterals,height=direction,layers=:graded,pins=(1,2,3,4))
            with_model(f,carriers->begin
                API.mesh.generate(3);expected=actual_carriers(carriers)
                API.mesh.set_order(2)
                initial=check_quadratic(carriers,f,expected)
                nodes=node_geometry();cells=cell_geometry(carriers.volume)
                before=snapshot()
                @test_throws ArgumentError API.mesh.renumber_nodes([1,2],[1,1])
                unchanged(before)
                API.mesh.reverse([(3,carriers.volume)])
                @test node_geometry()==nodes && cell_geometry(carriers.volume)==cells
                API.mesh.reverse([(3,carriers.volume)])
                primary_count=length(API.LAST_MESH_CLASS[].node_entities)
                API.mesh.renumber_nodes(collect(1:primary_count),collect(primary_count:-1:1))
                expected=actual_carriers(carriers)
                @test check_quadratic(carriers,f,expected)==initial
                for msh in API.mesh.get_element_types(3,carriers.volume)
                    tags=API.mesh.get_elements_by_type(msh,carriers.volume)[1]
                    API.mesh.reorder_elements(msh,carriers.volume,collect(length(tags)-1:-1:0))
                end
                @test node_geometry()==nodes && cell_geometry(carriers.volume)==cells
                types,tags,_=API.mesh.get_elements(3,carriers.volume)
                tag=first(first(tags));msh,connectivity,_,_=API.mesh.get_element(tag)
                removed=(Int(msh),Tuple(sort!([Tuple(API.mesh.get_node(v)[1]) for v in connectivity])))
                API.mesh.remove_elements(3,carriers.volume,[tag])
                @test cell_geometry(carriers.volume)==setdiff(cells,Set([removed]))
                API.mesh.clear();@test isempty(API.mesh.get_elements()[1])
            end)
        end
    end

    @testset "Coincident independent CAD carriers remain distinct" begin
        left=Cert.fixture("left";laterals=false,layers=:three)
        right=Cert.fixture("right";laterals=true,layers=:three,tags=:sparse,offset=(3.,0.,0.))
        pair=replace(left.source,"sweep[]="=>"left[]=")*replace(right.source,"sweep[]="=>"right[]=")
        API.initialize()
        try
            mktempdir() do directory
                path=joinpath(directory,"independent_patches.geo");write(path,pair)
                opened=API.open_geo!(path;mesh_dim=0)
                for (tag,p) in collect(API.CURRENT[].points)
                    p[1]>=3 || continue
                    API.model.set_coordinates(tag,p[1]-3,p[2],p[3])
                end
                outputs=[Int.(opened.lists[name]) for name in ("left","right")]
                carriers=[(;source=f.surface,top=out[1],volume=out[2],laterals=out[3:end])
                          for (f,out) in zip((left,right),outputs)]
                API.mesh.generate(3)
                expected=actual_carriers.(carriers)
                API.mesh.set_order(2)
                first_nodes=API.mesh.get_nodes(3,carriers[1].volume,true)[1]
                second_nodes=API.mesh.get_nodes(3,carriers[2].volume,true)[1]
                @test isempty(intersect(first_nodes,second_nodes))
                @test length(API.mesh.get_nodes()[1])==2*(54*3)+27+59
                for (c,f,e) in zip(carriers,(left,right),expected)
                    check_quadratic(c,f,e)
                end
                before=(node_geometry(carriers[2].volume),
                        cell_geometry(carriers[2].volume))
                API.mesh.clear([(3,carriers[1].volume)])
                @test isempty(API.mesh.get_elements(3,carriers[1].volume)[1])
                # Selective clear compacts the legacy dense cache. Compare the
                # surviving actual geometry, carriers and parameters by region.
                @test (node_geometry(carriers[2].volume),
                       cell_geometry(carriers[2].volume))==before
            end
        finally
            API.finalize()
        end
    end

    @testset "Refined child supports retain the actual CAD carriers" begin
        for laterals in (false,true),initial_order in (1,2)
            f=Cert.fixture("refine";laterals,layers=:one,pins=(1,2,3,4))
            with_model(f,carriers->begin
                API.mesh.generate(3)
                initial_order==2 && API.mesh.set_order(2)
                before_types,before_tags,_=API.mesh.get_elements(3,carriers.volume)
                parent_counts=Dict(Int(msh)=>length(ids) for (msh,ids) in zip(before_types,before_tags))
                # P1 and P2 represent the same parents. The pinned subdivision
                # contract has eight Tet/Hex/Pri children; one Pyramid produces
                # eight Tet and four Pyramid children.
                parent(msh)=get(parent_counts,initial_order==2 ? Dict(4=>11,5=>12,6=>13,7=>14)[msh] : msh,0)
                expected_children=Dict(4=>8*(parent(4)+parent(7)),5=>8parent(5),6=>8parent(6),7=>4parent(7))
                filter!(pair->last(pair)>0,expected_children)
                API.mesh.refine()
                types,tags,_=API.mesh.get_elements(3,carriers.volume)
                @test Dict(Int(msh)=>length(ids) for (msh,ids) in zip(types,tags))==
                    expected_children
                API.mesh.set_order(2)
                supports=support_nodes(carriers.volume)
                # Unit-square carriers are finite line segments and rectangular
                # surface patches. Entire primary supports must lie on a carrier;
                # no endpoint ownership or coordinate welding is used.
                cad=[(dim,entity,API.model.get_bounding_box(dim,entity)) for dim in (1,2)
                     for (_,entity) in API.model.get_entities(dim)]
                for (node,support) in supports
                    points=[API.mesh.get_node(v)[1] for v in support]
                    owner=(3,carriers.volume)
                    for (dim,entity,bounds) in cad
                        all(all(bounds[d]<=p[d]<=bounds[d+3] for d in 1:3) for p in points) || continue
                        owner=(dim,entity);break
                    end
                    @test API.mesh.get_node(node)[3:4]==owner
                end
                # Refining the width-five strip doubles each original edge:
                # its P2 cap/long patch has 17x5 nodes, with 15x3 interiors.
                # The width-two short patch has 5x5 nodes and 3x3 interiors.
                for surface in carriers.laterals
                    curve=abs(API.CURRENT[].meshing.extrude_sources[(2,surface)][2])
                    long=curve in (f.curve_tags[i] for i in f.long_pair)
                    own=check_parameters(2,surface,false,long ? 45 : 9)
                    closure=check_parameters(2,surface,true,long ? 85 : 25)
                    @test own⊆closure
                end
                for surface in (carriers.source,carriers.top)
                    check_parameters(2,surface,false,45)
                    check_parameters(2,surface,true,85)
                end
            end)
        end
    end

    @testset "Unrepresentable new P2 placement rejects atomically" begin
        for laterals in (false,true)
            base=Cert.fixture("P2_precision";laterals,layers=:one,offset=(0.,0.,2.0^50+2.))
            n=laterals ? 2 : 4
            f=merge(base,(;source=replace(base.source,"Layers{1}"=>"Layers{$n}"),levels=collect(range(0.,1.;length=n+1))))
            with_model(f,carriers->begin
                linear=API.mesh.generate(3)
                source=Tessella.Model.mesh_model_surface(API.CURRENT[],f.surface)
                @test Cert.certify(linear,f,source).total==1
                before=snapshot()
                @test_throws ArgumentError API.mesh.set_order(2)
                unchanged(before)
            end)
        end
    end

    @testset "Physical metadata survives actual order and file payloads" begin
        for laterals in (false,true)
            base=Cert.fixture("metadata";laterals,layers=:graded,tags=:sparse)
            f=merge(base,(;source=base.source*"""
                Physical Surface("four-strip source",501)={$(base.surface)};
                Physical Surface("four-strip cap",502)={sweep[0]};
                Physical Volume("four-strip body",503)={sweep[1]};
                """))
            with_model(f,carriers->begin
                API.mesh.generate(3)
                expected=Dict((2,501)=>"four-strip source",(2,502)=>"four-strip cap",(3,503)=>"four-strip body")
                @test API.model.get_physical_groups()==sort!(collect(keys(expected)))
                @test API.model.get_entities_for_physical_group(2,501)==[carriers.source]
                @test API.model.get_entities_for_physical_group(2,502)==[carriers.top]
                @test API.model.get_entities_for_physical_group(3,503)==[carriers.volume]
                # Complete owned metadata is the model projection contract;
                # the legacy API3 cache contains the actual volume cell set.
                linear=Tessella.Model.mesh_model_volume(API.CURRENT[],carriers.volume)
                projected=Tessella.Model.model_to_mixed(API.CURRENT[],linear,3,carriers.volume)
                @test projected.physical_names==expected
                for order in (2,1,2)
                    API.mesh.set_order(order)
                    @test API.model.get_physical_groups()==sort!(collect(keys(expected)))
                    for (entity,name) in expected
                        @test API.model.get_physical_name(entity...)==name
                    end
                    for node in API.mesh.get_nodes(3,carriers.volume)[1]
                        @test API.mesh.get_node(node)[3:4]==(3,carriers.volume)
                        @test API.model.get_physical_groups_for_entity(3,carriers.volume)==[503]
                    end
                end
                # MSH4 preserves node ownership and parameters; MSH2 carries
                # the actual cell ownership and first physical membership.
                @test all(all(==(503),b.tags) for b in projected.blocks if msh_dimension(b.msh)==3)
                mktempdir() do directory
                    check_files(projected,directory,"four_strip_classified";names=expected)
                    # The bare API3 cache declares its actual volume cells;
                    # classification is queried separately through the API.
                    # Serialize every actual P2 node and connectivity slot.
                    actual=API.mesh.get()
                    check_files(actual,directory,"four_strip_p2_cache")
                end
            end)
        end
    end

    @testset "Unsupported source changes and invalid selections preserve the live cache" begin
        for laterals in (false,true)
            f=Cert.fixture("cache_atomicity";laterals,layers=:one)
            with_model(f,carriers->begin
                API.mesh.generate(3);API.mesh.set_order(2)
                before=snapshot()
                @test_throws ArgumentError API.mesh.clear([(3,carriers.volume),(3,90001)])
                unchanged(before)
                for i in f.long_pair
                    API.mesh.set_transfinite_curve(f.curve_tags[i],6)
                end
                # Free M5 remains unsupported. Recombined arbitrary strips
                # are supported, so reject unmatched opposite source chains.
                laterals && API.mesh.set_transfinite_curve(f.curve_tags[last(f.long_pair)],7)
                before=snapshot()
                @test_throws r"QuadTriNoNewVerts" API.mesh.generate(3)
                unchanged(before)
            end)
        end
    end

    @testset "Invalid four-strip generation keeps geometry, options and allocator history" begin
        for laterals in (false,true)
            base=Cert.fixture("construction_atomicity";laterals,layers=:one)
            # These sources have valid CAD before meshing. Their unsupported
            # transform, incomplete normalized layers or missing source policy
            # must fail inside the owned generation transaction.
            invalid=(
                replace(base.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{0.125,0.0,1.0}"),
                replace(base.source,"Layers{1}"=>"Layers{{1},{0.5}}"),
                replace(base.source,"Recombine Surface{$(base.surface)};"=>""),
            )
            for source in invalid
                @test source!=base.source
                with_model(merge(base,(;source)),carriers->begin
                    API.option("Mesh.ElementOrder",2)
                    before=snapshot()
                    @test_throws r"QuadTriNoNewVerts" API.mesh.generate(3)
                    unchanged(before)
                end)
            end
            # The source and translated CAD endpoints are representable, but
            # intermediate requested planes round onto an adjacent plane.
            collapsed=Cert.fixture("collapsed_planes";laterals,layers=:three,
                height=.25,offset=(0.,0.,2.0^50+2.))
            with_model(collapsed,carriers->begin
                before=snapshot()
                @test_throws r"QuadTriNoNewVerts" API.mesh.generate(3)
                unchanged(before)
            end)
            with_model(base,carriers->begin
                API.mesh.generate(3);API.mesh.set_order(2)
                API.mesh.set_transfinite_surface(carriers.top)
                before=snapshot()
                @test_throws r"constrained-boundary|override" API.mesh.generate(3)
                unchanged(before)
            end)
        end
    end

    @testset "Actual lower supports survive an OnlyEmpty one-dimensional generation" begin
        for laterals in (false,true),initial_order in (1,2)
            f=Cert.fixture("lower_bridge";laterals,layers=:graded)
            with_model(f,carriers->begin
                API.mesh.generate(3)
                old_primary=Set(API.mesh.get_nodes()[1])
                initial_order==2 && API.mesh.set_order(2)
                before=Dict((dim,tag)=>API.mesh.get_nodes(dim,tag,false,true)
                    for dim in (0,1) for (_,tag) in API.model.get_entities(dim))
                high=(API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag())
                # Immediate set_order does not change the generation option.
                API.option("Mesh.ElementOrder",initial_order)
                API.option("Mesh.MeshOnlyEmpty",1)
                API.mesh.generate(1)
                new_primary=lower_primary_tags()
                for (entity,payload) in before
                    check_lower_regeneration(entity,payload,API.mesh.get_nodes(entity...,false,true),
                        initial_order,old_primary,new_primary)
                end
                types,tags,_=API.mesh.get_elements(1)
                @test !isempty(types) && sum(length,tags)==20+4(length(f.levels)-1)
                @test API.mesh.get_max_node_tag()>=high[1]
                @test API.mesh.get_max_element_tag()>=high[2]
                before=snapshot()
                @test_throws ArgumentError API.mesh.clear([(1,first(f.curve_tags)),(1,90001)])
                unchanged(before)
            end)
        end
    end
end
end

get(ENV,"TESSELLA_FOUR_STRIP_DEFINE_ONLY","")=="1" || run_tests()

end
