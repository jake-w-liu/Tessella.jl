if !isdefined(@__MODULE__,:QuadTriNoNewTriangleCertificates)
    include("../geometry/quadtri_nonew_triangle_certificates.jl")
end

module NoNewTwoTriBoundaryTests

using Test, Tessella, SHA
using Tessella.MeshTypes: nnodes
using ..QuadTriNoNewTriangleCertificates

const API = Tessella.API
const PROFILES = ((:one, "Layers{1}", (0., 1.)),
                  (:three, "Layers{3}", (0., 1/3, 2/3, 1.)),
                  (:graded, "Layers{{2,1},{0.25,1}}", (0., .125, .25, 1.)))

function fixture(layers, laterals, direction)
    policy = laterals ? " RecombLaterals" : ""
    return """
    Mesh.ElementOrder=1;Geometry.OldRuledSurface=0;
    Point(1)={0,0,0};Point(2)={1,0,0};Point(3)={1,1,0};Point(4)={0,1,0};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
    Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
    Transfinite Curve{1,2,3,4}=2;
    Transfinite Surface{1}={1,2,3,4} Left;
    out[]=Extrude{0,0,$direction}{Surface{1};$layers;Recombine;
        QuadTriNoNewVerts$policy;};
    """
end

function with_model(source, action)
    API.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"two_tri_boundary.geo")
            write(path,source)
            geometry=API.open_geo!(path;mesh_dim=0)
            out=Int.(geometry.lists["out"])
            action((;source=1,top=out[1],volume=out[2],laterals=out[3:end]))
        end
    finally
        API.finalize()
    end
end

function snapshot()
    return (;model=API.CURRENT[],value=repr(API.CURRENT[]),
        cache=API.LAST_MESH[],class=API.LAST_MESH_CLASS[],
        overlay=API.LAST_MESH_HIGH_ORDER[],mids=copy(API.LAST_MESH_HIGH_ORDER_MIDS[]),
        nodes=API.mesh.get_nodes(),cells=API.mesh.get_elements(),
        history=(API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[]))
end

function unchanged(before)
    @test API.CURRENT[]===before.model && repr(API.CURRENT[])==before.value
    @test API.LAST_MESH[]===before.cache && API.LAST_MESH_CLASS[]===before.class
    @test API.LAST_MESH_HIGH_ORDER[]===before.overlay
    @test API.LAST_MESH_HIGH_ORDER_MIDS[]==before.mids
    @test API.mesh.get_nodes()==before.nodes && API.mesh.get_elements()==before.cells
    @test (API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[])==before.history
end

function node_geometry()
    return Set(begin
        p,uv,dim,entity=API.mesh.get_node(node)
        (Tuple(p),dim,entity,Tuple(uv))
    end for node in API.mesh.get_nodes()[1])
end

function cell_geometry(volume)
    types,tags,connectivity=API.mesh.get_elements(3,volume)
    result=Set{Tuple{Int,Tuple}}()
    for (msh,ids,nodes) in zip(types,tags,connectivity)
        width=Tessella.Elements.msh_spec(msh).nnodes
        for column in eachindex(ids)
            points=[Tuple(API.mesh.get_node(nodes[(column-1)*width+row])[1])
                    for row in 1:width]
            push!(result,(Int(msh),Tuple(sort!(points))))
        end
    end
    return result
end

function support_nodes(msh,volume)
    supports=Dict{UInt64,Tuple}()
    edges=reshape(API.mesh.get_element_edge_nodes(msh,volume,false),3,:)
    for edge in eachcol(edges)
        key=Tuple(sort!(collect(edge[1:2])))
        @test !haskey(supports,edge[3]) || supports[edge[3]]==key
        supports[edge[3]]=key
    end
    if msh==13
        faces=reshape(API.mesh.get_element_face_nodes(msh,4,volume,false),9,:)
        for face in eachcol(faces)
            key=Tuple(sort!(collect(face[1:4])))
            @test !haskey(supports,face[9]) || supports[face[9]]==key
            supports[face[9]]=key
        end
    end
    return supports
end

function actual_carriers()
    # This bounded fixture has only axis-aligned straight curves and planar
    # rectangular surfaces. Use their actual CAD tags and bounding boxes,
    # together with the element's explicit primary support, to find carriers.
    # No node identities are merged or inferred from endpoint owners.
    return [(dim,entity,API.model.get_bounding_box(dim,entity))
            for dim in (1,2) for (_,entity) in API.model.get_entities(dim)]
end

function support_carrier(support,carriers,volume)
    points=[API.mesh.get_node(node)[1] for node in support]
    for (dim,entity,bounds) in carriers
        if all(all(bounds[i]<=p[i]<=bounds[i+3] for i in 1:3) for p in points)
            return (dim,entity)
        end
    end
    return (3,volume)
end

function check_supports(msh,volume)
    carriers=actual_carriers()
    supports=support_nodes(msh,volume)
    signature=Dict{Tuple,Tuple{Int,Int}}()
    for (node,support) in supports
        p,_,dim,owner=API.mesh.get_node(node)
        expected=support_carrier(support,carriers,volume)
        @test (dim,owner)==expected
        points=[API.mesh.get_node(primary)[1] for primary in support]
        mean=[sum(point[i] for point in points)/length(points) for i in 1:3]
        @test maximum(abs.(p.-mean))<=2e-11
        key=Tuple(sort!([Tuple(point) for point in points]))
        signature[key]=expected
    end
    return supports,signature
end

function check_parameters(dim,entity,boundary,expected)
    nodes,coordinates,uv=API.mesh.get_nodes(dim,entity,boundary,true)
    @test length(nodes)==expected && allunique(nodes)
    @test length(coordinates)==3expected && length(uv)==dim*expected
    @test maximum(abs.(API.model.get_value(dim,entity,uv).-coordinates))<=2e-11
    return nodes
end

function check_quadratic(carriers,levels,laterals,direction)
    n=length(levels)-1
    msh=laterals ? 13 : 11
    types,tags,_=API.mesh.get_elements(3,carriers.volume)
    @test types==[msh]
    @test sum(length,tags)==(laterals ? 2n : 6n)
    allnodes=API.mesh.get_nodes()[1]
    @test length(allnodes)==18n+9 && allunique(allnodes)
    dimensions=zeros(Int,4)
    for node in allnodes
        p,uv,dim,entity=API.mesh.get_node(node)
        dimensions[dim+1]+=1
        if dim in (1,2)
            @test length(uv)==dim
            @test maximum(abs.(API.model.get_value(dim,entity,uv).-p))<=2e-11
        else
            @test isempty(uv)
        end
    end
    @test dimensions==[8,8n+4,8n-2,2n-1]
    @test Set(Int.(last.(API.model.get_entities(2))))==
          Set((carriers.source,carriers.top,carriers.laterals...))
    @test length(carriers.laterals)==4
    for surface in carriers.laterals
        own=check_parameters(2,surface,false,2n-1)
        closure=check_parameters(2,surface,true,6n+3)
        @test Set(own) ⊆ Set(closure)
        for node in own
            @test API.mesh.get_node(node)[3:4]==(2,surface)
        end
    end
    for (surface,z) in ((carriers.source,0.),(carriers.top,Float64(direction)))
        own=check_parameters(2,surface,false,1)
        @test Tuple(API.mesh.get_node(only(own))[1])==(.5,.5,z)
        check_parameters(2,surface,true,9)
    end
    volume_nodes,xyz,_=API.mesh.get_nodes(3,carriers.volume,false,true)
    @test length(volume_nodes)==2n-1
    expected_heights=sort!([direction*z for z in
        (collect(levels[2:end-1])...,((levels[i]+levels[i+1])/2 for i in 1:n)...)])
    actual=reshape(xyz,3,:)
    @test all(actual[1,:].==.5) && all(actual[2,:].==.5)
    @test maximum(abs.(sort!(collect(actual[3,:])).-expected_heights))<=2e-11
    supports,signature=check_supports(msh,carriers.volume)
    @test length(supports)==14n+5
    @test count(s->length(s)==2,values(supports))==(laterals ? 9n+5 : 14n+5)
    @test count(s->length(s)==4,values(supports))==(laterals ? 5n : 0)
    return signature
end

artifact_name(profile,laterals,direction,order)=
    "two_tri_$(profile)_R$(laterals)_D$(direction)_P$(order)"

function artifact_record(name,order)
    # Keep the existing strict public node/cell/owner/parameter serializer
    # unchanged. Extend its digest with actual integer support catalogs and
    # computed own/closure surface queries, which have distinct UV contracts.
    base=QuadTriNoNewTriangleCertificates.api_record(name,API)
    buffer=IOBuffer()
    println(buffer,"two_tri_support_query_v1");println(buffer,base.sha)
    class=API.LAST_MESH_CLASS[]
    for (label,catalog) in (("edge",class.edge_entities),
                            ("tri",class.face_entities),("quad",class.quad_entities))
        println(buffer,label);write(buffer,UInt64(length(catalog)))
        for (support,owner) in sort!(collect(catalog);by=first)
            write(buffer,UInt64(length(support)))
            for node in support;write(buffer,UInt64(node));end
            write(buffer,Int64(owner[1]),Int64(owner[2]))
        end
    end
    for (_,surface) in sort(API.model.get_entities(2)),boundary in (false,true)
        tags,coordinates,parameters=API.mesh.get_nodes(2,surface,boundary,true)
        write(buffer,UInt64(surface),UInt8(boundary),UInt64(length(tags)))
        for node in tags;write(buffer,UInt64(node));end
        for array in (coordinates,parameters)
            write(buffer,UInt64(length(array)))
            for value in array;write(buffer,reinterpret(UInt64,Float64(value)));end
        end
    end
    return merge(base,(;kind="api_p$(order)_supports",sha=bytes2hex(sha256(take!(buffer)))))
end

function artifact_pins()
    path=joinpath(@__DIR__,"..","artifacts","quadtri_nonew_two_tri_crc.txt")
    records=Dict{String,NamedTuple}()
    for line in readlines(path)
        isempty(line) || startswith(line,"#") || begin
            fields=split(line,'\t')
            length(fields)==6 || error("two-Tri CRC artifact row width")
            name,kind,n_nodes,n_cells,families,sha=String.(fields)
            haskey(records,name) && error("duplicate two-Tri CRC artifact name")
            records[name]=(;name,kind,n_nodes=parse(Int,n_nodes),n_cells=parse(Int,n_cells),families,sha)
        end
    end
    length(records)==24 || error("two-Tri CRC artifact must contain 24 products")
    return records
end

@testset "Two-Tri NoNew actual P2 carriers and public lifecycle" begin
    pins=artifact_pins()
    @testset "Twelve normal-translation ownership fixtures" begin
        for (name,layers,levels) in PROFILES,laterals in (false,true),direction in (-1,1)
            @testset "$name RecombLaterals=$laterals direction=$direction" begin
                with_model(fixture(layers,laterals,direction),carriers->begin
                    API.mesh.generate(3)
                    @test length(API.mesh.get_nodes()[1])==4length(levels)
                    @test API.mesh.get_element_types(3,carriers.volume)==[laterals ? 6 : 4]
                    record_name=artifact_name(name,laterals,direction,1)
                    @test artifact_record(record_name,1)==pins[record_name]
                    API.mesh.set_order(2)
                    record_name=artifact_name(name,laterals,direction,2)
                    @test artifact_record(record_name,2)==pins[record_name]
                    initial=check_quadratic(carriers,levels,laterals,direction)
                    before=snapshot()
                    API.mesh.set_order(2)
                    @test API.mesh.get_nodes()==before.nodes && API.mesh.get_elements()==before.cells
                    @test check_quadratic(carriers,levels,laterals,direction)==initial
                    API.mesh.set_order(1)
                    @test length(API.mesh.get_nodes()[1])==4length(levels)
                    @test API.mesh.get_element_types(3,carriers.volume)==[laterals ? 6 : 4]
                    API.mesh.set_order(2)
                    @test check_quadratic(carriers,levels,laterals,direction)==initial
                end)
            end
        end
    end

    @testset "Explicit identity maps survive permutation and exact selection" begin
        for laterals in (false,true),direction in (-1,1)
            _,layers,levels=PROFILES[3]
            with_model(fixture(layers,laterals,direction),carriers->begin
                API.mesh.generate(3);API.mesh.set_order(2)
                initial=check_quadratic(carriers,levels,laterals,direction)
                initial_nodes=node_geometry();initial_cells=cell_geometry(carriers.volume)
                primary_count=laterals ? nnodes(API.LAST_MESH[]) : 4length(levels)
                before=snapshot()
                @test_throws ArgumentError API.mesh.renumber_nodes(collect(1:primary_count),fill(1,primary_count))
                unchanged(before)
                for route in (:reverse,:renumber_nodes,:renumber_elements,:reorder)
                    if route===:reverse
                        API.mesh.reverse([(3,carriers.volume)])
                    elseif route===:renumber_nodes
                        API.mesh.renumber_nodes(collect(1:primary_count),collect(primary_count:-1:1))
                    elseif route===:renumber_elements
                        tags=sort(vcat(API.mesh.get_elements(3,carriers.volume)[2]...))
                        API.mesh.renumber_elements(tags,reverse(tags))
                    else
                        msh=laterals ? 13 : 11
                        count=sum(length,API.mesh.get_elements(3,carriers.volume)[2])
                        API.mesh.reorder_elements(msh,carriers.volume,collect(count-1:-1:0))
                    end
                    @test node_geometry()==initial_nodes
                    @test cell_geometry(carriers.volume)==initial_cells
                    @test check_quadratic(carriers,levels,laterals,direction)==initial
                end
                types,tags,_=API.mesh.get_elements(3,carriers.volume)
                removed_type,removed_nodes,_,_=API.mesh.get_element(first(only(tags)))
                removed=(Int(removed_type),Tuple(sort!([
                    Tuple(API.mesh.get_node(node)[1]) for node in removed_nodes])))
                API.mesh.remove_elements(3,carriers.volume,[first(only(tags))])
                @test sum(length,API.mesh.get_elements(3,carriers.volume)[2])==length(only(tags))-1
                @test cell_geometry(carriers.volume)==setdiff(initial_cells,Set([removed]))
                _,selected=check_supports(only(types),carriers.volume)
                @test all(get(initial,key,nothing)==owner for (key,owner) in selected)
            end)
        end
    end

    @testset "Refinement preserves actual parent carriers on child supports" begin
        for laterals in (false,true),direction in (-1,1),initial_order in (1,2)
            with_model(fixture("Layers{1}",laterals,direction),carriers->begin
                API.mesh.generate(3)
                initial_order==2 && API.mesh.set_order(2)
                API.mesh.refine()
                # Exact Mixed refinement returns P1; the supported legacy
                # simplex refinement retains its documented P2 query overlay.
                @test API.mesh.get_element_types(3,carriers.volume)==
                    [laterals ? 6 : initial_order==2 ? 11 : 4]
                API.mesh.set_order(2)
                msh=laterals ? 13 : 11
                @test sum(length,API.mesh.get_elements(3,carriers.volume)[2])==(laterals ? 16 : 48)
                _,supports=check_supports(msh,carriers.volume)
                @test count(owner->owner[1]==1,values(supports))>0
                @test count(owner->owner[1]==2,values(supports))>0
                @test count(owner->owner[1]==3,values(supports))>0
                for surface in carriers.laterals
                    own,xyz,uv=API.mesh.get_nodes(2,surface,false,true)
                    closure,closed_xyz,closed_uv=API.mesh.get_nodes(2,surface,true,true)
                    @test length(own)==9 && length(closure)==25
                    @test length(uv)==18 && length(closed_uv)==50
                    @test Set(own) ⊆ Set(closure)
                    @test maximum(abs.(API.model.get_value(2,surface,uv).-xyz))<=2e-11
                    @test maximum(abs.(API.model.get_value(2,surface,closed_uv).-closed_xyz))<=2e-11
                end
                API.mesh.reverse()
                @test last(check_supports(msh,carriers.volume))==supports
            end)
        end
    end

    @testset "One existing triangular free source retains refined classification" begin
        source="""
        Mesh.TransfiniteTri=1;
        Point(1)={0,0,0};Point(2)={1,0,0};Point(3)={0,1,0};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,1};
        Curve Loop(1)={1,2,3};Plane Surface(1)={1};
        Transfinite Curve{1,2,3}=2;Transfinite Surface{1};
        out[]=Extrude{0,0,1}{Surface{1};Layers{1};Recombine;QuadTriNoNewVerts;};
        """
        with_model(source,carriers->begin
            API.mesh.generate(3);API.mesh.set_order(2);API.mesh.refine()
            @test API.LAST_MESH_CLASS[]!==nothing
            @test API.mesh.get_element_types(3,carriers.volume)==[11]
            @test sum(length,API.mesh.get_elements(3,carriers.volume)[2])==24
            for surface in carriers.laterals
                check_parameters(2,surface,false,9)
                check_parameters(2,surface,true,25)
            end
            curve_nodes,xyz,uv=API.mesh.get_nodes(1,1,false,true)
            @test length(curve_nodes)==3 && length(uv)==3
            @test Set(Tuple(p) for p in eachcol(reshape(xyz,3,:)))==
                Set(((.25,0.,0.),(.5,0.,0.),(.75,0.,0.)))
            @test maximum(abs.(API.model.get_value(1,1,uv).-xyz))<=2e-11
        end)
    end

    @testset "Unsupported source grids reject without publishing partial state" begin
        for laterals in (false,true)
            with_model(fixture("Layers{3}",laterals,1),carriers->begin
                API.mesh.generate(3);API.mesh.set_order(2)
                for curve in 1:4;API.mesh.set_transfinite_curve(curve,3);end
                before=snapshot()
                err=try API.mesh.generate(3);nothing catch caught;caught end
                @test err isa ArgumentError
                @test occursin("QuadTriNoNewVerts",sprint(showerror,err))
                unchanged(before)
            end)
        end
    end
end

end # module
