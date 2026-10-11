module RigidPeriodicHybridFunctionalChecks
using Test,LinearAlgebra,Tessella
const A=Tessella.API
const E=Tessella.Elements
const RZ90=[0.,-1.,0.,1., 1.,0.,0.,0., 0.,0.,1.,1., 0.,0.,0.,1.]
const RZ180=[-1.,0.,0.,1., 0.,-1.,0.,1., 0.,0.,1.,1., 0.,0.,0.,1.]

function cube(;affine)
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
        tag!=2 && A.mesh.set_transfinite_surface(tag)
        A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1)
    A.model.add_volume([1];tag=1)
    A.mesh.set_periodic(2,[2],[1],affine)
end

snapshot()=(A.mesh.get_nodes(),A.mesh.get_elements(),
            A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())

function map_point(affine,point)
    return [sum(affine[4(axis-1)+component]*getfield(point isa Tuple ? point : Tuple(point),component)
                for component in 1:3)+affine[4axis] for axis in 1:3]
end

function correspondence(affine,high_order,expected)
    link=A.mesh.get_periodic_nodes(2,2,high_order)
    @test link.master_entity==1 && link.affine==affine
    @test length(link.slave_nodes)==length(link.master_nodes)==expected
    @test allunique(link.slave_nodes) && allunique(link.master_nodes)
    @test all(first(A.mesh.get_node(slave))≈map_point(affine,first(A.mesh.get_node(master)))
              for (slave,master) in zip(link.slave_nodes,link.master_nodes))
    return Dict(zip(link.master_nodes,link.slave_nodes))
end

function verify_cube(affine)
    families,labels,vertices=A.mesh.get_elements(3,1)
    @test Set(families)==Set(Int32[4,7])
    cells=Dict(Int(kind)=>reshape(nodes,E.msh_num_nodes(kind),length(ids))
               for (kind,ids,nodes) in zip(families,labels,vertices))
    @test size(cells[4],2)==24 && size(cells[7],2)==6
    ids,xyz,_=A.mesh.get_nodes()
    points=Dict(uid=>xyz[3row-2:3row] for (row,uid) in enumerate(ids))
    function tet_volume(nodes)
        a,b,c,d=(points[uid] for uid in nodes)
        det(hcat(b-a,c-a,d-a))/6
    end
    volumes=[tet_volume(cell) for cell in eachcol(cells[4])]
    append!(volumes,[tet_volume(cell[[1,2,3,5]])+tet_volume(cell[[1,3,4,5]])
                     for cell in eachcol(cells[7])])
    @test all(>(0),volumes)
    @test sum(volumes)≈1.0 atol=2e-10
    incidence=Dict{Tuple,Int}()
    for (kind,nodes) in cells
        faces=kind==4 ? ((1,2,3),(1,2,4),(1,3,4),(2,3,4)) :
            ((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5))
        for cell in eachcol(nodes),face in faces
            key=Tuple(sort(cell[collect(face)]))
            incidence[key]=get(incidence,key,0)+1
        end
    end
    quads=Tuple[]
    for surface in 1:6
        kinds,tags,nodes=A.mesh.get_elements(2,surface)
        @test kinds==Int32[3] && length(only(tags))==1
        push!(quads,Tuple(sort(only(nodes))))
    end
    @test all(count->count in (1,2),values(incidence))
    @test Set(face for (face,count) in incidence if count==1)==Set(quads)
    @test Set(Tuple(sort(cell[1:4])) for cell in eachcol(cells[7]))==Set(quads)
    pairs=correspondence(affine,false,4)
    _,master_tags,master_nodes=A.mesh.get_elements(2,1)
    _,slave_tags,slave_nodes=A.mesh.get_elements(2,2)
    @test map(uid->pairs[uid],only(master_nodes))==only(slave_nodes)
    @test isempty(intersect(only(master_tags),only(slave_tags)))
    return cells
end

@testset "quarter-turn copied Quad boundary admits positive hybrid cube" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        cube(;affine=RZ90)
        A.mesh.generate(3)
        cells=verify_cube(RZ90)
        # The slave face's rotated copy binds the actual master rows: slave
        # corner uids map onto rotated master corner coordinates exactly.
        link=A.mesh.get_periodic_nodes(2,2,false)
        @test all(zip(link.slave_nodes,link.master_nodes)) do (slave,master)
            sm,pm=first(A.mesh.get_node(slave)),first(A.mesh.get_node(master))
            sm==map_point(RZ90,pm)
        end
        # Surface-derived cross-side curve masters are enforced by the owned
        # lower producer: every slave edge node is the affine image of its
        # master edge node, so side faces share one set of vertices.
        @test Tessella.validate(A.mesh.get()).ok
        mktempdir() do directory
            path=joinpath(directory,"rigid-cube.msh")
            A.write(path)
            restored=E.read_mixed_msh(path)
            @test Tessella.validate(restored).ok
            link=only(row for row in restored.periodic_links if row.dim==2 && row.slave_entity==2)
            @test link.master_entity==1 && collect(link.affine)==RZ90
            data=restored.entity_data
            ids=first(A.mesh.get_nodes())
            @test Set(data.external_node_tags)==Set(ids)
        end
        A.mesh.set_order(2)
        @test Set(first(A.mesh.get_elements(3,1)))==Set(Int32[11,14])
        correspondence(RZ90,false,4)
        correspondence(RZ90,true,9)
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "half-turn copied Quad boundary admits positive hybrid cube" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        cube(;affine=RZ180)
        A.mesh.generate(3)
        verify_cube(RZ180)
        link=A.mesh.get_periodic_nodes(2,2,false)
        @test all(zip(link.slave_nodes,link.master_nodes)) do (slave,master)
            sm,pm=first(A.mesh.get_node(slave)),first(A.mesh.get_node(master))
            sm==map_point(RZ180,pm)
        end
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "non-rigid periodic affines reject atomically" begin
    shear=[1.,0.,0.25,0.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.]
    reflection=[-1.,0.,0.,1.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.]
    for affine in (shear,reflection)
        A.initialize()
        try
            A.option("General.Terminal",0)
            cube(;affine)
            A.mesh.generate(2)
            before=snapshot()
            @test_throws ArgumentError A.mesh.generate(3)
            @test snapshot()==before
        finally
            A.finalize()
        end
    end
end

@testset "directly declared cross-surface curve relation keeps its guard" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
                (0,0,1),(1,0,1),(1,1,1),(0,1,1))
        for (tag,p) in enumerate(points)
            A.model.add_point(p...;tag,meshSize=1.)
        end
        edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
               (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
        for (tag,ends) in enumerate(edges)
            A.model.add_line(ends...;tag)
        end
        loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
               (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
        for (tag,lines) in enumerate(loops)
            A.model.add_curve_loop(collect(lines);tag)
            A.model.add_plane_surface([tag];tag)
        end
        # A directly declared (non-derived) relation pairing Curve 5 — a
        # boundary of Surfaces 2 and 3 — to Curve 1 off Surface 2 stays an
        # explicit admission failure; only Surface-induced cross-side
        # relations route through the owned lower producer.
        A.mesh.set_periodic(1,[5],[1],
            [1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.])
        error=try
            A.mesh.generate(2);nothing
        catch failure
            failure isa InterruptException && rethrow();failure
        end
        @test error isa ArgumentError
        @test occursin("has only one entity",sprint(showerror,error))
    finally
        A.finalize()
    end
end
end
