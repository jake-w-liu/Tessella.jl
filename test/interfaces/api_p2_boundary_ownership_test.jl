using Test,Tessella
using Tessella.MeshTypes: nnodes

if !isdefined(@__MODULE__,:QuadTriNoNewTriangleCertificates)
    include("../geometry/quadtri_nonew_triangle_certificates.jl")
end
const _P2BOQ=QuadTriNoNewTriangleCertificates

function _p2bo_with_model(source,action)
    api=Tessella.API
    api.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"p2_boundary.geo");write(path,source)
            geometry=api.open_geo!(path;mesh_dim=0)
            action(api,geometry)
        end
    finally
        api.finalize()
    end
end

function _p2bo_node_geometry(api,nodes)
    return Set(begin
        p,uv,dim,entity=api.mesh.get_node(node)
        (Tuple(p),dim,entity,Tuple(uv))
    end for node in nodes)
end

function _p2bo_lateral_queries(api,surfaces)
    for surface in surfaces
        own,coordinates,uv=api.mesh.get_nodes(2,surface,false,true)
        @test length(own)==5 && length(uv)==10
        @test maximum(abs.(api.model.get_value(2,surface,uv).-coordinates))<=2e-11
        for node in own
            p,parameters,dim,entity=api.mesh.get_node(node)
            @test (dim,entity)==(2,surface) && length(parameters)==2
            @test maximum(abs.(api.model.get_value(2,surface,parameters).-p))<=2e-11
        end
        closure,coordinates,uv=api.mesh.get_nodes(2,surface,true,true)
        @test length(closure)==21 && length(uv)==42 && allunique(closure)
        @test Set(own) ⊆ Set(closure)
        @test maximum(abs.(api.model.get_value(2,surface,uv).-coordinates))<=2e-11
    end
end

function _p2bo_snapshot(api)
    return (;model=api.CURRENT[],value=repr(api.CURRENT[]),
        cache=api.LAST_MESH[],class=api.LAST_MESH_CLASS[],
        overlay=api.LAST_MESH_HIGH_ORDER[],mids=copy(api.LAST_MESH_HIGH_ORDER_MIDS[]),
        nodes=api.mesh.get_nodes(),cells=api.mesh.get_elements(),
        history=(api.NODE_TAG_MAX[],api.ELEMENT_TAG_MAX[]))
end

function _p2bo_atomic(api,before)
    @test api.CURRENT[]===before.model && repr(api.CURRENT[])==before.value
    @test api.LAST_MESH[]===before.cache && api.LAST_MESH_CLASS[]===before.class
    @test api.LAST_MESH_HIGH_ORDER[]===before.overlay
    @test api.LAST_MESH_HIGH_ORDER_MIDS[]==before.mids
    @test api.mesh.get_nodes()==before.nodes && api.mesh.get_elements()==before.cells
    @test (api.NODE_TAG_MAX[],api.ELEMENT_TAG_MAX[])==before.history
end

function _p2bo_mixed_catalog_fixture(api,coincident_edges)
    for (surface,offset,points) in (
            (1,0,((0.,0.),(1.,0.),(0.,1.))),
            (2,100,((3.,0.),(4.,0.),(3.,1.))),
            (3,200,((10.,0.),(11.,0.),(11.,1.),(10.,1.))))
        for (i,p) in enumerate(points)
            api.model.add_point(p[1],p[2],0;tag=offset+i)
        end
        for i in eachindex(points)
            api.model.add_line(offset+i,offset+mod1(i+1,length(points));tag=offset+i)
            api.mesh.set_transfinite_curve(offset+i,2)
        end
        api.model.add_curve_loop(collect(offset+1:offset+length(points));tag=surface)
        api.model.add_plane_surface([surface];tag=surface)
        api.mesh.set_transfinite_surface(surface)
    end
    api.option("Mesh.TransfiniteTri",1)
    api.mesh.set_recombine(2,3)
    for tag in (coincident_edges ? (101:103) : (101:101))
        p=api.CURRENT[].points[tag]
        api.model.set_coordinates(tag,p[1]-3,p[2],p[3])
    end
    api.mesh.generate(2)
    @test api.LAST_MESH[] isa Tessella.Elements.MixedMesh
    @test length(api.LAST_MESH_CLASS[].edge_entities)==10
end

function _p2bo_support_coordinates(api)
    mesh=api.LAST_MESH[]
    result=Dict{Tuple{Int,Int32},Set{NTuple{3,Float64}}}()
    for (edge,owner) in api.LAST_MESH_CLASS[].edge_entities
        @test owner[1]==1
        points=Set(Tuple(mesh.coords[:,node]) for node in edge)
        xyz=reshape(api.model.get_value(1,Int(owner[2]),[0.,1.]),3,:)
        @test points==Set(Tuple(point) for point in eachcol(xyz))
        result[owner]=points
    end
    return result
end

@testset "Actual simplex P2 boundary support ownership" begin
    @testset "Mixed compaction remaps exact supports before publication" begin
        api=Tessella.API
        try
            api.initialize()
            _p2bo_mixed_catalog_fixture(api,false)
            before=_p2bo_support_coordinates(api)
            count=length(api.mesh.get_nodes()[1])
            api.mesh.remove_duplicate_nodes()
            @test length(api.mesh.get_nodes()[1])==count-1
            @test _p2bo_support_coordinates(api)==before
            @test length(api.LAST_MESH_CLASS[].face_entities)==2
            api.mesh.set_order(2)
            @test _p2bo_support_coordinates(api)==before
            class=api.LAST_MESH_CLASS[]
            for block in api.LAST_MESH[].blocks
                block.msh==9 || continue
                for column in axes(block.nodes,2),
                    (slot,i,j) in ((4,1,2),(5,2,3),(6,3,1))
                    support=minmax(block.nodes[i,column],block.nodes[j,column])
                    @test class.node_entities[block.nodes[slot,column]]==class.edge_entities[support]
                end
            end
        finally
            api.finalize()
        end
        try
            api.initialize()
            _p2bo_mixed_catalog_fixture(api,true)
            before=_p2bo_snapshot(api)
            err=try api.mesh.remove_duplicate_nodes();nothing catch e;e end
            @test err isa ArgumentError
            @test occursin("ambiguous entity ownership",sprint(showerror,err))
            _p2bo_atomic(api,before)
        finally
            api.finalize()
        end
    end
    @testset "Hex27 face centers keep actual Surface carriers after Point merge" begin
        api=Tessella.API
        try
            api.initialize()
            api.model.add_box(0,0,0,1,1,1)
            api.model.add_box(1,1,1,1,1,1)
            for (_,curve) in api.model.get_entities(1)
                api.mesh.set_transfinite_curve(curve,2)
            end
            for (_,surface) in api.model.get_entities(2)
                api.mesh.set_transfinite_surface(surface)
                api.mesh.set_recombine(2,surface)
            end
            for (_,volume) in api.model.get_entities(3)
                api.mesh.set_transfinite_volume(volume)
            end
            api.mesh.generate(3)
            @test length(api.mesh.get_nodes()[1])==16
            @test length(api.LAST_MESH_CLASS[].quad_entities)==12
            api.mesh.remove_duplicate_nodes()
            @test length(api.mesh.get_nodes()[1])==15
            api.mesh.set_order(2)
            @test length(api.mesh.get_nodes()[1])==53
            for (_,surface) in api.model.get_entities(2)
                bounds=api.model.get_bounding_box(2,surface)
                expected=ntuple(i->(bounds[i]+bounds[i+3])/2,3)
                own,xyz,uv=api.mesh.get_nodes(2,surface,false,true)
                @test length(own)==1 && Tuple(xyz)==expected
                @test length(uv)==2
                point,parameters,dim,entity=api.mesh.get_node(only(own))
                @test (dim,entity)==(2,surface) && Tuple(point)==expected
                @test api.model.get_value(2,surface,parameters)==point
            end
            for volume in 1:2
                own,xyz,_=api.mesh.get_nodes(3,volume)
                @test length(own)==1
                @test Tuple(xyz)==(volume-.5,volume-.5,volume-.5)
            end
        finally
            api.finalize()
        end
    end
    @testset "Tri6 curve supports and includeBoundary" begin
        source="""
        Point(1)={0,0,0};Point(2)={1,0,0};Point(3)={1,1,0};Point(4)={0,1,0};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Transfinite Curve{1,2,3,4}=2;Transfinite Surface{1};
        """
        _p2bo_with_model(source,(api,_)->begin
            api.mesh.generate(2);api.mesh.set_order(2)
            @test api.LAST_MESH_HIGH_ORDER[] isa Tessella.HighOrder.P2TriMesh
            @test length(api.mesh.get_nodes()[1])==9
            mids=((.5,0.,0.),(1.,.5,0.),(.5,1.,0.),(0.,.5,0.))
            for (curve,expected) in enumerate(mids)
                own,xyz,uv=api.mesh.get_nodes(1,curve,false,true)
                @test length(own)==1 && Tuple(xyz)==expected && length(uv)==1
                p,parameters,dim,entity=api.mesh.get_node(only(own))
                @test (dim,entity)==(1,curve) && Tuple(p)==expected
                @test api.model.get_value(1,curve,parameters)≈p
                closure,xyz,uv=api.mesh.get_nodes(1,curve,true,true)
                @test length(closure)==3 && length(uv)==3 && allunique(closure)
                @test maximum(abs.(api.model.get_value(1,curve,uv).-xyz))<=2e-11
            end
            own,xyz,uv=api.mesh.get_nodes(2,1,false,true)
            @test length(own)==1 && Tuple(xyz)==(.5,.5,0.) && length(uv)==2
            closure,xyz,uv=api.mesh.get_nodes(2,1,true,true)
            @test length(closure)==9 && length(uv)==18 && allunique(closure)
            @test maximum(abs.(api.model.get_value(2,1,uv).-xyz))<=2e-11
            before=_p2bo_node_geometry(api,api.mesh.get_nodes()[1])
            api.mesh.reverse()
            @test _p2bo_node_geometry(api,api.mesh.get_nodes()[1])==before
            api.mesh.renumber_nodes(collect(1:4),collect(4:-1:1))
            @test _p2bo_node_geometry(api,api.mesh.get_nodes()[1])==before
            for curve in 1:4
                @test length(api.mesh.get_nodes(1,curve,true)[1])==3
            end
        end)
    end

    @testset "Coincident Tet10 regions retain separate carriers through mutations" begin
        # Both regions are fully simplex, so this exercises the legacy overlay
        # rather than the native MixedMesh high-order representation.
        source=replace(_P2BOQ.paired_source(),
            "QuadTriNoNewVerts RecombLaterals"=>"QuadTriNoNewVerts")
        _p2bo_with_model(source,(api,geometry)->begin
            for (tag,p) in collect(api.CURRENT[].points)
                p[1]>=3 || continue
                api.model.set_coordinates(tag,p[1]-3,p[2],p[3])
            end
            cache=api.mesh.generate(3)
            @test cache isa Tessella.MeshTypes.Mesh && nnodes(cache)==24
            api.mesh.set_order(2)
            @test api.LAST_MESH_HIGH_ORDER[] isa Tessella.HighOrder.P2Mesh
            left=Int(geometry.lists["left"][2]);right=Int(geometry.lists["right"][2])
            surfaces=(Int.(geometry.lists["left"][3:end]),Int.(geometry.lists["right"][3:end]))
            function inspect_regions()
                a=api.mesh.get_nodes(3,left,true)[1]
                b=api.mesh.get_nodes(3,right,true)[1]
                @test length(a)==42 && length(b)==42 && isempty(intersect(a,b))
                @test isempty(api.mesh.get_nodes(3,left,false)[1])
                @test isempty(api.mesh.get_nodes(3,right,false)[1])
                for group in surfaces;_p2bo_lateral_queries(api,group);end
                return (_p2bo_node_geometry(api,a),_p2bo_node_geometry(api,b))
            end
            initial=inspect_regions()
            @test length(api.mesh.get_nodes()[1])==84
            @test Set(first.(initial[1]))==Set(first.(initial[2]))
            before=_p2bo_snapshot(api)
            @test_throws ArgumentError api.mesh.renumber_nodes(collect(1:24),fill(1,24))
            _p2bo_atomic(api,before)
            other=first(api.mesh.get_elements(3,right)[2][1])
            @test_throws ArgumentError api.mesh.remove_elements(3,left,[other])
            _p2bo_atomic(api,before)
            for route in (:reverse,:renumber_nodes,:renumber_elements,:reorder)
                if route===:reverse
                    api.mesh.reverse()
                elseif route===:renumber_nodes
                    api.mesh.renumber_nodes(collect(1:24),collect(24:-1:1))
                elseif route===:renumber_elements
                    tags=sort(vcat(api.mesh.get_elements(3)[2]...))
                    api.mesh.renumber_elements(tags,reverse(tags))
                else
                    api.mesh.reorder_elements(11,right,collect(8:-1:0))
                end
                @test inspect_regions()==initial
            end
            api.mesh.remove_elements(3,left)
            @test isempty(api.mesh.get_elements(3,left)[1])
            @test sum(length,api.mesh.get_elements(3,right)[2])==9
            surviving=api.mesh.get_nodes(3,right,true)[1]
            @test length(surviving)==42
            @test _p2bo_node_geometry(api,surviving)==initial[2]
            _p2bo_lateral_queries(api,surfaces[2])
        end)
    end

    @testset "Refined multi-region carriers include child-face supports" begin
        source=replace(_P2BOQ.paired_source(),
            "QuadTriNoNewVerts RecombLaterals"=>"QuadTriNoNewVerts")
        _p2bo_with_model(source,(api,geometry)->begin
            api.mesh.generate(3);api.mesh.set_order(2);api.mesh.refine()
            @test sum(length,api.mesh.get_elements(3)[2])==144
            for name in ("left","right"),surface in Int.(geometry.lists[name][3:end])
                own,xyz,uv=api.mesh.get_nodes(2,surface,false,true)
                @test length(own)==33 && length(uv)==66 && allunique(own)
                @test maximum(abs.(api.model.get_value(2,surface,uv).-xyz))<=2e-11
                closure,xyz,uv=api.mesh.get_nodes(2,surface,true,true)
                @test length(closure)==65 && length(uv)==130 && allunique(closure)
                @test Set(own) ⊆ Set(closure)
                @test maximum(abs.(api.model.get_value(2,surface,uv).-xyz))<=2e-11
            end
            entities=Int.([geometry.lists["left"][2],geometry.lists["right"][2]])
            nodes=[api.mesh.get_nodes(3,entity,true)[1] for entity in entities]
            @test isempty(intersect(nodes...))
            before=[_p2bo_node_geometry(api,group) for group in nodes]
            api.mesh.reverse()
            @test [_p2bo_node_geometry(api,api.mesh.get_nodes(3,entity,true)[1])
                for entity in entities]==before
        end)
    end
end
