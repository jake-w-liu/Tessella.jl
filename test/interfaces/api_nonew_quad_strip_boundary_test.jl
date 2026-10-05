if !isdefined(@__MODULE__,:QuadTriNoNewQuadStripCertificates)
    include("../geometry/quadtri_nonew_quad_strip_certificates.jl")
end
if !isdefined(@__MODULE__,:QuadTriNoNewTriangleCertificates)
    include("../geometry/quadtri_nonew_triangle_certificates.jl")
end

module NoNewQuadStripBoundaryTests

using Test,Tessella,SHA,TOML
using Tessella.Elements: msh_dimension
using ..QuadTriNoNewQuadStripCertificates
using ..QuadTriNoNewTriangleCertificates
const API=Tessella.API
const Cert=QuadTriNoNewQuadStripCertificates
const PROFILES=(:one,:three,:graded)

function with_model(f,action)
    API.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"quad_strip_boundary.geo");write(path,f.source)
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
        history=(API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[]))
end

function unchanged(before)
    @test API.CURRENT[]===before.model && repr(before.model)==before.value
    @test before.model.curve_params==before.params
    @test API.LAST_MESH[]===before.cache && API.LAST_MESH_CLASS[]===before.class
    @test API.LAST_MESH_HIGH_ORDER[]===before.overlay
    @test API.LAST_MESH_EDGES[]===before.edges && API.LAST_MESH_FACES[]===before.faces
    @test (API.mesh.get_nodes(),API.mesh.get_elements())==before.payload
    @test (API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[])==before.history
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
    @test length(volume_nodes)==30n+(f.laterals ? 31 : 15)
    counts=zeros(Int,4)
    for node in volume_nodes
        p,uv,dim,owner=API.mesh.get_node(node)
        counts[dim+1]+=1
        @test dim in (1,2) ? length(uv)==dim : isempty(uv)
        dim in (1,2) && (@test maximum(abs.(API.model.get_value(dim,owner,uv).-p))<=2e-11)
    end
    @test counts==[8,8n+12,16n-2,6n+(f.laterals ? 13 : -3)]
    types,tags,_=API.mesh.get_elements(3,carriers.volume)
    families=Dict(Int(t)=>length(ids) for (t,ids) in zip(types,tags))
    # Free choices are constrained by the complete joined relation, rather
    # than a pointer-dependent family signature. Every actual elevated cell
    # and node is checked below; recombined terminal fans have fixed counts.
    if f.laterals
        correct=Dict(12=>2(n-1),11=>4,14=>10)
        filter!(pair->last(pair)>0,correct)
        @test families==correct
    end
    supports=support_nodes(carriers.volume)
    @test length(supports)==24n+(f.laterals ? 23 : 9)
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
        own=check_parameters(2,surface,false,long ? 6n-3 : 2n-1)
        closure=check_parameters(2,surface,true,long ? 10n+5 : 6n+3)
        @test own⊆closure
    end
    for surface in (carriers.source,carriers.top)
        own=check_parameters(2,surface,false,3)
        closure=check_parameters(2,surface,true,15)
        @test own⊆closure
    end
    for curve in f.curve_tags
        long=curve in (f.curve_tags[i] for i in f.long_pair)
        check_parameters(1,curve,false,long ? 3 : 1)
        check_parameters(1,curve,true,long ? 5 : 3)
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
    "quad_strip_$(profile)_R$(laterals)_D$(direction)_P$(order)"

function artifact_record(name,order)
    base=QuadTriNoNewTriangleCertificates.api_record(name,API)
    buffer=IOBuffer();println(buffer,"quad_strip_support_query_v1");println(buffer,base.sha)
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
            API.mesh.generate(3)
            push!(records,artifact_record(artifact_name(profile,laterals,direction,1),1))
            API.mesh.set_order(2)
            push!(records,artifact_record(artifact_name(profile,laterals,direction,2),2))
        end)
    end
    return records
end

function artifact_pins()
    path=joinpath(@__DIR__,"..","artifacts","quadtri_nonew_quad_strip_crc.txt")
    records=Dict{String,NamedTuple}()
    for line in readlines(path)
        isempty(line) || startswith(line,"#") || begin
            fields=split(line,'\t');length(fields)==6 || error("quad strip artifact row width")
            name,kind,n_nodes,n_cells,families,sha=String.(fields)
            haskey(records,name) && error("quad strip duplicate artifact name")
            records[name]=(;name,kind,n_nodes=parse(Int,n_nodes),n_cells=parse(Int,n_cells),families,sha)
        end
    end
    length(records)==24 || error("quad strip artifact must contain24 products")
    return records
end

function saved_oracles()
    path=joinpath(@__DIR__,"..","artifacts","quadtri_nonew_quad_strip_oracle.toml")
    artifact=TOML.parsefile(path)
    @test artifact["gmsh_version"]=="4.15.2" && artifact["fixture_count"]==16
    return Dict((f["profile"],f["laterals"],f["direction"])=>f
                for f in artifact["fixtures"] if f["capture_kind"]=="original_unit")
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

function run_tests()
@testset "Quad-strip actual P2 carriers and public lifecycle" begin
    @testset "Twelve P1/P2 native golden products" begin
        pins=artifact_pins()
        oracle=saved_oracles()
        for profile in PROFILES,laterals in (false,true),direction in (-1,1)
            f=Cert.fixture("golden";layers=profile,laterals,height=direction,pins=(1,2,3,4))
            with_model(f,carriers->begin
                API.mesh.generate(3)
                @test length(API.mesh.get_nodes()[1])==6length(f.levels)+(f.laterals ? 2 : 0)
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
                @test length(API.mesh.get_nodes()[1])==6length(f.levels)+(f.laterals ? 2 : 0)
                API.mesh.set_order(2)
                @test check_quadratic(carriers,f,expected)==initial
            end)
        end
    end

    @testset "Opposite winding, axes, curved-free CAD sampling and sparse entity tags" begin
        for plane in (:XY,:YZ,:ZX),laterals in (false,true)
            f=Cert.fixture("owners";plane,laterals,height=-1.,layers=:graded,shape=:rounded,
                           winding=-1,pins=:opposite,curve_reverse=5,tags=:sparse,
                           progression=4.,long_pair=plane===:YZ ? (2,4) : (1,3))
            with_model(f,carriers->begin
                API.mesh.generate(3);expected=actual_carriers(carriers)
                API.mesh.set_order(2);check_quadratic(carriers,f,expected)
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
                @test length(API.mesh.get_nodes()[1])==2*(30*3)+15+31
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
                for surface in carriers.laterals
                    curve=abs(API.CURRENT[].meshing.extrude_sources[(2,surface)][2])
                    long=curve in (f.curve_tags[i] for i in f.long_pair)
                    own=check_parameters(2,surface,false,long ? 21 : 9)
                    closure=check_parameters(2,surface,true,long ? 45 : 25)
                    @test own⊆closure
                end
                for surface in (carriers.source,carriers.top)
                    check_parameters(2,surface,false,21)
                    check_parameters(2,surface,true,45)
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
                Physical Surface("strip source",501)={$(base.surface)};
                Physical Surface("strip cap",502)={sweep[0]};
                Physical Volume("strip body",503)={sweep[1]};
                """))
            with_model(f,carriers->begin
                API.mesh.generate(3)
                expected=Dict((2,501)=>"strip source",(2,502)=>"strip cap",(3,503)=>"strip body")
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
                # MSH4 owns cell carriers in entity_data; its reader need not
                # duplicate that information into the optional MSH2 field.
                signature(mesh)=sort!([(b.msh,Tuple(cell),b.tags[i],mesh.entity_data.block_entities[j][i])
                    for (j,b) in enumerate(mesh.blocks) for (i,cell) in enumerate(eachcol(b.nodes))])
                @test all(all(==(503),b.tags) for b in projected.blocks if msh_dimension(b.msh)==3)
                mktempdir() do directory
                    for binary in (false,true)
                        path=joinpath(directory,"strip.msh")
                        Tessella.Elements.write_mixed_msh(path,projected;version=4.1,binary)
                        back=Tessella.Elements.read_mixed_msh(path)
                        @test back.physical_names==expected
                        @test back.coords==projected.coords
                        @test signature(back)==signature(projected)
                        @test back.entity_data.node_entities==projected.entity_data.node_entities
                        @test back.entity_data.node_parametric==projected.entity_data.node_parametric
                    end
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
                # Arbitrary free and recombined strips are supported;
                # reject unmatched opposite source chains.
                API.mesh.set_transfinite_curve(f.curve_tags[last(f.long_pair)],7)
                before=snapshot()
                @test_throws r"QuadTriNoNewVerts" API.mesh.generate(3)
                unchanged(before)
            end)
        end
    end
end
end

get(ENV,"TESSELLA_STRIP_DEFINE_ONLY","")=="1" || run_tests()

end
