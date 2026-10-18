# PER02 slice: dimension-0 periodic queries and independent multi-direction
# periodic pairs. Run standalone:
#   julia --project=. --startup-file=no -O0 test/geometry/periodic_dim0_multidir_test.jl
# (`--compile=min` cannot run this path: the recombination libm shims use
# `llvmcall`, which requires the compiler.)
module PeriodicDim0MultidirChecks

using Test
using Tessella
using Tessella.MeshTypes: validate
const A=Tessella.API

# affine maps a MASTER entity onto its slave (Gmsh row-major 4x4).
const TX=Float64[1,0,0,1, 0,1,0,0, 0,0,1,0, 0,0,0,1]
const TY=Float64[1,0,0,0, 0,1,0,1, 0,0,1,0, 0,0,0,1]
const TZ=Float64[1,0,0,0, 0,1,0,0, 0,0,1,1, 0,0,0,1]

function _geo(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

const _SQUARE="""
    Point(1) = {0, 0, 0, 0.5}; Point(2) = {1, 0, 0, 0.5};
    Point(3) = {1, 1, 0, 0.5}; Point(4) = {0, 1, 0, 0.5};
    Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
    Curve Loop(1) = {1, 2, 3, 4};
    Plane Surface(1) = {1};
    """

const _DOUBLE_PERIODIC_SQUARE=_SQUARE*"""
    Periodic Curve {2} = {4} Translate {1, 0, 0};
    Periodic Curve {3} = {1} Translate {0, 1, 0};
    """

function _api_square()
    A.model.add_point(0,0,0;tag=1,meshSize=0.5)
    A.model.add_point(1,0,0;tag=2,meshSize=0.5)
    A.model.add_point(1,1,0;tag=3,meshSize=0.5)
    A.model.add_point(0,1,0;tag=4,meshSize=0.5)
    A.model.add_line(1,2;tag=1);A.model.add_line(2,3;tag=2)
    A.model.add_line(3,4;tag=3);A.model.add_line(4,1;tag=4)
    A.model.add_curve_loop([1,2,3,4];tag=1)
    A.model.add_plane_surface([1];tag=1)
    A.mesh.set_periodic(1,[2],[4],TX)
    A.mesh.set_periodic(1,[3],[1],TY)
    return nothing
end

function _api_cube()
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
            (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,p) in enumerate(points)
        A.model.add_point(p...;tag,meshSize=1.)
    end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
           (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(edges)
        A.model.add_line(ends...;tag)
        A.mesh.set_transfinite_curve(tag,2)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,lines) in enumerate(loops)
        A.model.add_curve_loop(collect(lines);tag)
        A.model.add_plane_surface([tag];tag)
        tag in (1,3,6) && A.mesh.set_transfinite_surface(tag)
        A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1)
    A.model.add_volume([1];tag=1)
    # Three independent directions: x=1 face from x=0, y=1 from y=0,
    # z=1 from z=0 — one master per slave face, shared corner vertices.
    A.mesh.set_periodic(2,[4],[6],TX)
    A.mesh.set_periodic(2,[5],[3],TY)
    A.mesh.set_periodic(2,[2],[1],TZ)
    return nothing
end

@testset "two periodic directions share a corner without duplicate nodes" begin
    built=_geo(_DOUBLE_PERIODIC_SQUARE;mesh_dim=2)
    @test built.msg_error_count==0
    @test validate(built.mesh).ok
    @test [(Int(c.slave_entity),Int(c.master_entity))
           for c in model_periodic_constraints(built.model)]==[(2,4),(3,1)]
    mesh=built.mesh
    for (slave,master,shift) in ((2,4,(1.0,0.0,0.0)),(3,1,(0.0,1.0,0.0)))
        mapping=model_periodic_nodes(built.model,mesh,1,slave)
        @test mapping.master_entity==master
        @test !isempty(mapping.slave_nodes)
        for (sn,mn) in zip(mapping.slave_nodes,mapping.master_nodes)
            @test Tuple(mesh.coords[:,sn])==Tuple(mesh.coords[:,mn]).+shift
        end
    end
    # Vertex relations resolve through the dim-1 endpoint graph — the same
    # spanning-forest pairs the $Periodic dim-0 links serialize. On this
    # square the shared corner (point 3) takes its y-direction parent 2; the
    # x-direction edge (3<-4) stays a curve-level record only.
    for (slave_point,master_point,shift) in
            ((2,1,(1.0,0.0,0.0)),(4,1,(0.0,1.0,0.0)),(3,2,(0.0,1.0,0.0)))
        mapping=model_periodic_nodes(built.model,mesh,0,slave_point)
        @test mapping.master_entity==master_point
        @test mapping.affine!==nothing
        @test length(mapping.slave_nodes)==length(mapping.master_nodes)==1
        @test Tuple(mesh.coords[:,only(mapping.slave_nodes)])==
              built.model.points[slave_point]
        @test Tuple(mesh.coords[:,only(mapping.master_nodes)])==
              built.model.points[master_point]
    end
    @test_throws ArgumentError model_periodic_nodes(built.model,mesh,0,1)
    mixed=model_to_mixed(built.model,geo_entity_mesh(built,2,1),1)
    dim0=[link for link in mixed.periodic_links if link.dim==0]
    dim1=[link for link in mixed.periodic_links if link.dim==1]
    @test length(dim0)==3 && length(dim1)==2
    @test Dict{Int,Int}(Int(link.slave_entity)=>Int(link.master_entity)
                        for link in dim0)==Dict(2=>1,4=>1,3=>2)
    @test allunique(vcat((Vector(link.slave_nodes) for link in dim0)...))
    @test all(link->length(link.slave_nodes)==1,dim0)
end

@testset "API periodic vertex queries and same-slave rejection" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        _api_square()
        # Pre-mesh: vertex masters resolve through the curve endpoint graph
        # (upstream reports the stored `GVertex` master and affine the same
        # way); node lists stay empty until a mesh exists.
        @test A.mesh.get_periodic(0,[1,2,3,4])==Int32[1,1,2,1]
        pre=A.mesh.get_periodic_nodes(0,3)
        @test pre.master_entity==2 && isempty(pre.slave_nodes) &&
              isempty(pre.master_nodes) && pre.affine==TY
        # generate(1) exercises the periodic curve copy plus the endpoint
        # vertex correspondence — generate(2) over declared boundary curve
        # relations is a separate pre-existing authority-path gap (the same
        # `.geo` surface succeeds through `mesh_model_surface`).
        generated=A.mesh.generate(1)
        @test validate(generated).ok
        @test A.mesh.get_periodic(0,[1,2,3,4])==Int32[1,1,2,1]
        post=A.mesh.get_periodic_nodes(0,3)
        @test post.master_entity==2
        @test length(post.slave_nodes)==length(post.master_nodes)==1
        @test first(A.mesh.get_node(only(post.slave_nodes)))==[1.0,1.0,0.0]
        @test first(A.mesh.get_node(only(post.master_nodes)))==[1.0,0.0,0.0]
        @test post.affine==TY
        lone=A.mesh.get_periodic_nodes(0,1)
        @test lone.master_entity==1 && isempty(lone.slave_nodes) &&
              isempty(lone.master_nodes) && isempty(lone.affine)
        curve_link=A.mesh.get_periodic_nodes(1,2)
        @test curve_link.master_entity==4 && !isempty(curve_link.slave_nodes)
        @test_throws ArgumentError A.mesh.get_periodic_nodes(0,99)
        @test_throws ArgumentError A.mesh.get_periodic_nodes(0,1.5)
        # A single slave may not gain a second master through the strict API —
        # the unsupported multi-master direction rejects loudly rather than
        # silently replacing or dropping the stored relation.
        @test_throws ArgumentError A.mesh.set_periodic(1,[2],[3],TX)
        @test A.mesh.get_periodic(0,[1,2,3,4])==Int32[1,1,2,1]
    finally
        A.finalize()
    end
end

@testset "three-direction periodic box volumes through the hybrid boundary" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        _api_cube()
        A.mesh.generate(3)
        @test validate(A.mesh.get()).ok
        # Each pair answers independently: four corner correspondences each.
        for (slave,master,affine,shift) in
                ((4,6,TX,(1.,0.,0.)),(5,3,TY,(0.,1.,0.)),(2,1,TZ,(0.,0.,1.)))
            link=A.mesh.get_periodic_nodes(2,slave)
            @test link.master_entity==master
            @test link.affine==affine
            @test length(link.slave_nodes)==length(link.master_nodes)==4
            @test allunique(link.slave_nodes) && allunique(link.master_nodes)
            @test all(zip(link.slave_nodes,link.master_nodes)) do (s,m)
                first(A.mesh.get_node(s))≈first(A.mesh.get_node(m)).+collect(shift)
            end
        end
        # Corner vertices resolve single masters through the derived curve
        # graph — no duplicate vertex nodes anywhere.
        @test A.mesh.get_periodic(0,collect(1:8))==Int32[1,1,2,1,1,2,3,4]
        @test validate(A.mesh.get()).ok
        families,labels,vertices=A.mesh.get_elements(3,1)
        @test Set(families)==Set(Int32[4,7])
        ids,xyz,_=A.mesh.get_nodes()
        points=Dict(uid=>xyz[3row-2:3row] for (row,uid) in enumerate(ids))
        function tet_volume(nodes)
            a,b,c,d=(points[uid] for uid in nodes)
            return abs(((b[1]-a[1])*((c[2]-a[2])*(d[3]-a[3])-(c[3]-a[3])*(d[2]-a[2]))+
                        (b[2]-a[2])*((c[3]-a[3])*(d[1]-a[1])-(c[1]-a[1])*(d[3]-a[3]))+
                        (b[3]-a[3])*((c[1]-a[1])*(d[2]-a[2])-(c[2]-a[2])*(d[1]-a[1])))/6)
        end
        cells=Dict(Int(kind)=>reshape(nodes,
                Tessella.Elements.msh_num_nodes(kind),:)
                   for (kind,ids2,nodes) in zip(families,labels,vertices))
        volumes=[tet_volume(cell) for cell in eachcol(cells[4])]
        append!(volumes,[tet_volume(cell[[1,2,3,5]])+
                         tet_volume(cell[[1,3,4,5]])
                         for cell in eachcol(cells[7])])
        @test all(>(0),volumes)
        @test sum(volumes)≈1.0 atol=2e-10
    finally
        A.finalize()
    end
end

end
