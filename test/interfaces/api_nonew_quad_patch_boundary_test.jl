if !isdefined(@__MODULE__,:QuadTriNoNewQuadPatchCertificates)
    include("../geometry/quadtri_nonew_quad_patch_certificates.jl")
end
if !isdefined(@__MODULE__,:QuadTriNoNewTriangleCertificates)
    include("../geometry/quadtri_nonew_triangle_certificates.jl")
end

module NoNewQuadPatchBoundaryTests

using Test,Tessella,SHA
using Tessella.Elements: msh_dimension
using ..QuadTriNoNewQuadPatchCertificates
using ..QuadTriNoNewTriangleCertificates
const API=Tessella.API
const Cert=QuadTriNoNewQuadPatchCertificates
const PROFILES=(:one,:three,:graded)

function with_model(f,action)
    API.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"quad_patch_boundary.geo");write(path,f.source)
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
        if msh in (12,14)
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
    @test length(volume_nodes)==50n+25
    counts=zeros(Int,4)
    for node in volume_nodes
        p,uv,dim,owner=API.mesh.get_node(node)
        counts[dim+1]+=1
        @test dim in (1,2) ? length(uv)==dim : isempty(uv)
        dim in (1,2) && (@test maximum(abs.(API.model.get_value(dim,owner,uv).-p))<=2e-11)
    end
    @test counts==[8,8n+20,24n+6,18n-9]
    types,tags,_=API.mesh.get_elements(3,carriers.volume)
    families=Dict(Int(t)=>length(ids) for (t,ids) in zip(types,tags))
    correct=f.laterals ? Dict(12=>4(n-1),14=>12) : Dict(11=>24n-8,14=>4)
    filter!(pair->last(pair)>0,correct)
    @test families==correct
    supports=support_nodes(carriers.volume)
    @test length(supports)==41n+16
    @test count(s->length(s)==2,values(supports))==(f.laterals ? 21n+24 : 41n+12)
    @test count(s->length(s)==4,values(supports))==(f.laterals ? 16n-4 : 4)
    @test count(s->length(s)==8,values(supports))==(f.laterals ? 4(n-1) : 0)
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
        own=check_parameters(2,surface,false,6n-3)
        closure=check_parameters(2,surface,true,10n+5)
        @test own⊆closure
    end
    for surface in (carriers.source,carriers.top)
        own=check_parameters(2,surface,false,9)
        closure=check_parameters(2,surface,true,25)
        @test own⊆closure
    end
    for curve in f.curve_tags
        check_parameters(1,curve,false,3)
        check_parameters(1,curve,true,5)
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
    "quad_patch_$(profile)_R$(laterals)_D$(direction)_P$(order)"

function artifact_record(name,order)
    base=QuadTriNoNewTriangleCertificates.api_record(name,API)
    buffer=IOBuffer();println(buffer,"quad_patch_support_query_v1");println(buffer,base.sha)
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
    path=joinpath(@__DIR__,"..","artifacts","quadtri_nonew_quad_patch_crc.txt")
    records=Dict{String,NamedTuple}()
    for line in readlines(path)
        isempty(line) || startswith(line,"#") || begin
            fields=split(line,'\t');length(fields)==6 || error("quad patch artifact row width")
            name,kind,n_nodes,n_cells,families,sha=String.(fields)
            haskey(records,name) && error("quad patch duplicate artifact name")
            records[name]=(;name,kind,n_nodes=parse(Int,n_nodes),n_cells=parse(Int,n_cells),families,sha)
        end
    end
    length(records)==24 || error("quad patch artifact must contain24 products")
    return records
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

@testset "Quad-patch actual P2 carriers and public lifecycle" begin
    @testset "Twelve P1/P2 native golden products" begin
        pins=artifact_pins()
        for profile in PROFILES,laterals in (false,true),direction in (-1,1)
            f=Cert.fixture("golden";layers=profile,laterals,height=direction,pins=(1,2,3,4))
            with_model(f,carriers->begin
                API.mesh.generate(3)
                @test length(API.mesh.get_nodes()[1])==9length(f.levels)
                @test artifact_record(artifact_name(profile,laterals,direction,1),1)==pins[artifact_name(profile,laterals,direction,1)]
                expected=actual_carriers(carriers)
                API.mesh.set_order(2)
                @test artifact_record(artifact_name(profile,laterals,direction,2),2)==pins[artifact_name(profile,laterals,direction,2)]
                initial=check_quadratic(carriers,f,expected)
                before=snapshot();API.mesh.set_order(2)
                @test (API.mesh.get_nodes(),API.mesh.get_elements())==before.payload
                API.mesh.set_order(1)
                @test length(API.mesh.get_nodes()[1])==9length(f.levels)
                API.mesh.set_order(2)
                @test check_quadratic(carriers,f,expected)==initial
            end)
        end
    end

    @testset "Opposite winding, axes, curved-free CAD sampling and sparse entity tags" begin
        for plane in (:XY,:YZ,:ZX),laterals in (false,true)
            f=Cert.fixture("owners";plane,laterals,height=-1.,layers=:graded,shape=:rounded,
                           winding=-1,pins=:opposite,curve_reverse=5,tags=:sparse)
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
                @test length(API.mesh.get_nodes()[1])==350
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
                API.mesh.refine()
                types,tags,_=API.mesh.get_elements(3,carriers.volume)
                @test Dict(Int(msh)=>length(ids) for (msh,ids) in zip(types,tags))==
                    (laterals ? Dict(4=>96,7=>48) : Dict(4=>160,7=>16))
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
                    own=check_parameters(2,surface,false,21)
                    closure=check_parameters(2,surface,true,45)
                    @test own⊆closure
                end
                for surface in (carriers.source,carriers.top)
                    check_parameters(2,surface,false,49)
                    check_parameters(2,surface,true,81)
                end
            end)
        end
    end

    @testset "Unrepresentable new P2 placement rejects atomically" begin
        for laterals in (false,true)
            f=Cert.fixture("P2_precision";laterals,layers=:four,offset=(0.,0.,2.0^50+2.))
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
end

end
