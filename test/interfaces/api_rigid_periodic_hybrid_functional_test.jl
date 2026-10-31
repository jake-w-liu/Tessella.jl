module RigidPeriodicHybridFunctionalChecks
using Test,LinearAlgebra,Tessella
const A=Tessella.API
const E=Tessella.Elements
const RZ90=[0.,-1.,0.,1., 1.,0.,0.,0., 0.,0.,1.,1., 0.,0.,0.,1.]
const RZ180=[-1.,0.,0.,1., 0.,-1.,0.,1., 0.,0.,1.,1., 0.,0.,0.,1.]
const TX=[1.,0.,0.,1., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.]
const TY=[1.,0.,0.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.]
const TX2=[1.,0.,0.,2., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.]
const SCALE2=[2.,0.,0.,12., 0.,2.,0.,0., 0.,0.,2.,0., 0.,0.,0.,1.]

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

# Same hybrid cube as `cube`, but the Faces listed in `slaves` stay
# unstructured so independent periodic copies produce their Quads —
# pairs may be declared afterwards per direction.
function multi_cube(;slaves=(2,4,5))
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
    # face tags: 1=z0, 2=z1, 3=y0, 4=x1, 5=y1, 6=x0
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,lines) in enumerate(loops)
        A.model.add_curve_loop(collect(lines);tag)
        A.model.add_plane_surface([tag];tag)
        tag in slaves || A.mesh.set_transfinite_surface(tag)
        A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1)
    A.model.add_volume([1];tag=1)
end

# Two adjacent unit cubes sharing the wall Face 4 (x=1): Region 1 spans
# [0,1]^3 with Face tags 1..6 like `multi_cube`, Region 2 spans x in [1,2]
# with exterior Faces 101..105; `slaves` marks Faces left unstructured for
# periodic copies and `tris` marks Faces left unrecombined (a periodic
# slave's cells inherit their master's block type, so a Tri wall needs a
# Tri master).
function adjacent!(;slaves=Int[],tris=Int[])
    corners=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,point) in enumerate(corners);A.model.add_point(point...;tag,meshSize=1.);end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,pair) in enumerate(edges);A.model.add_line(pair...;tag);A.mesh.set_transfinite_curve(tag,2);end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),(2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,loop) in enumerate(loops)
        A.model.add_curve_loop(tag==1 ? collect(loop) : [-line for line in reverse(loop)];tag)
        A.model.add_plane_surface([tag];tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1);A.model.add_volume([1];tag=1)
    for (tag,point) in zip(101:104,((2,0,0),(2,1,0),(2,0,1),(2,1,1)))
        A.model.add_point(point...;tag,meshSize=1.)
    end
    more=((2,101),(101,102),(102,3),(6,103),(103,104),(104,7),(101,103),(102,104))
    for (tag,pair) in zip(101:108,more);A.model.add_line(pair...;tag);A.mesh.set_transfinite_curve(tag,2);end
    right=((101,102,103,-2),(104,105,106,-6),(101,107,-104,-10),(102,108,-105,-107),(103,11,-106,-108))
    for (index,loop) in enumerate(right)
        tag=100+index
        A.model.add_curve_loop(index==1 ? collect(loop) : [-line for line in reverse(loop)];tag)
        A.model.add_plane_surface([tag];tag)
    end
    A.model.add_surface_loop([101,102,103,104,105,-4];tag=101);A.model.add_volume([101];tag=2)
    for surface in vcat(collect(1:6),collect(101:105))
        surface in slaves || A.mesh.set_transfinite_surface(surface)
        (surface==4 || surface in tris) || A.mesh.set_recombine(2,surface)
    end
end
# Region 1 faces: 1(z0),2(z1),3(y0),4(wall x1),5(y1),6(x0)
# Region 2 faces: 101(z0),102(z1),103(y0),104(x2),105(y1),4(wall,-4)

# Two DISJOINT unit boxes: Region 1 = [0,1]^3, Region 2 = [10,12]x[0,2]x[0,2]
# — a scale-2 periodic pair joins Faces on disconnected Regions.
function two_boxes_scale!()
    pts1=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,p) in enumerate(pts1); A.model.add_point(p...;tag,meshSize=1.); end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(edges); A.model.add_line(ends...;tag); end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),(2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,lines) in enumerate(loops)
        A.model.add_curve_loop(collect(lines);tag)
        A.model.add_plane_surface([tag];tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1)
    A.model.add_volume([1];tag=1)
    pts2=((10,0,0),(12,0,0),(12,2,0),(10,2,0),(10,0,2),(12,0,2),(12,2,2),(10,2,2))
    for (tag,p) in enumerate(pts2); A.model.add_point(p...;tag=100+tag,meshSize=1.); end
    e2=((101,102),(102,103),(103,104),(104,101),(105,106),(106,107),
        (107,108),(108,105),(101,105),(102,106),(103,107),(104,108))
    for (tag,ends) in enumerate(e2); A.model.add_line(ends...;tag=100+tag); end
    l2=((101,102,103,104),(105,106,107,108),(101,110,-105,-109),
        (102,111,-106,-110),(103,112,-107,-111),(104,109,-108,-112))
    for (tag,lines) in enumerate(l2)
        A.model.add_curve_loop(collect(lines);tag=100+tag)
        A.model.add_plane_surface([100+tag];tag=100+tag)
    end
    A.model.add_surface_loop(collect(101:106);tag=101)
    A.model.add_volume([101];tag=2)
end
# box2 faces: 101=z0,102=z1,103=y0,104=x12,105=y2,106=x10

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

# Generalized link check for non-`cube` fixtures: the emitted dim-2 link
# carries the declared master tag and affine, and every emitted slave row
# sits exactly on its master's affine image.
function slave_link(slave,affine,master)
    link=A.mesh.get_periodic_nodes(2,slave,false)
    @test link.master_entity==master && link.affine==affine
    @test length(link.slave_nodes)==length(link.master_nodes)>0
    @test allunique(link.slave_nodes) && allunique(link.master_nodes)
    @test all(zip(link.slave_nodes,link.master_nodes)) do (s,mo)
        first(A.mesh.get_node(s))≈map_point(affine,first(A.mesh.get_node(mo)))
    end
    return link
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

@testset "non-rigid periodic affines copy under the declared map" begin
    # Gmsh 4.15.2 accepts and generates the same maps: the applied 3×4
    # transform copies verbatim under per-`GFace` `setMeshMaster` — the
    # declaration's corner/edge correspondence is the reachability gate, not
    # a rigid-body class. Both fixtures keep the master's Quads so the slave
    # copies exercise the hybrid boundary certificate.
    shear=[1.,0.,0.25,0.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.]
    reflection=[-1.,0.,0.,1.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.]
    for affine in (shear,reflection)
        A.initialize()
        try
            A.option("General.Terminal",0)
            cube(;affine)
            A.mesh.generate(3)
            verify_cube(affine)
            @test Tessella.validate(A.mesh.get()).ok
        finally
            A.finalize()
        end
    end
end

@testset "non-corresponding affine still rejects at declaration" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # (x,y,0) -> (x,y+0.25x,1) does not map the unit-square master onto
        # the slave Face — upstream rejects the same non-corresponding map
        # when the vertex correspondence cannot be resolved.
        offplane=[1.,0.,0.,0.,0.25,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.]
        @test_throws ArgumentError cube(;affine=offplane)
    finally
        A.finalize()
    end
end

@testset "independent periodic directions share a copied corner" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        multi_cube()
        A.mesh.set_periodic(2,[2],[1],RZ90)
        A.mesh.set_periodic(2,[4],[6],TX)
        A.mesh.set_periodic(2,[5],[3],TY)
        A.mesh.generate(3)
        families,_,_=A.mesh.get_elements(3,1)
        @test Set(families)==Set(Int32[4,7])
        links=(slave_link(2,RZ90,1),slave_link(4,TX,6),slave_link(5,TY,3))
        # The (1,1,1) corner is a slave vertex of all three slave Faces; each
        # relation binds that same actual vertex to a different master under
        # its own affine — joint consistency keys bindings by (uid,affine),
        # so distinct transforms at one corner are not a conflict.
        ids,xyz,_=A.mesh.get_nodes()
        corner=only([uid for (row,uid) in enumerate(ids)
                     if xyz[3row-2:3row]==[1.0,1.0,1.0]])
        bound=Set{eltype(ids)}()
        for link in links
            index=findfirst(==(corner),link.slave_nodes)
            @test index!==nothing
            push!(bound,link.master_nodes[index])
        end
        @test length(bound)==3
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "cross-Region pair joins Faces on different Regions" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        adjacent!(;slaves=[104])
        A.mesh.set_periodic(2,[104],[6],TX2)
        A.mesh.generate(3)
        for region in (1,2)
            families,_,_=A.mesh.get_elements(3,region)
            @test Set(families)==Set(Int32[4,7])
        end
        slave_link(104,TX2,6)
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "shared Tri wall admits a periodic slave copy" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # Tri master keeps the shared wall a Tri interface — upstream
        # generates the same fused two-Region case (recombined shared walls
        # stay the qualified-reject case covered below).
        adjacent!(;slaves=[4],tris=[6])
        A.mesh.set_periodic(2,[4],[6],TX)
        A.mesh.generate(3)
        slave_link(4,TX,6)
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "shared wall admits a periodic master" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        adjacent!(;slaves=[104])
        A.mesh.set_periodic(2,[104],[4],TX)
        A.mesh.generate(3)
        slave_link(104,TX,4)
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "recombined master on a shared wall keeps the shared-Quad reject" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # The slave's copy inherits the recombined master's Quads onto the
        # shared wall — upstream fails the same non-manifold quad-boundary
        # case, so the two-Region Pyramid contract stays a loud reject.
        adjacent!(;slaves=[4])
        A.mesh.set_periodic(2,[4],[6],TX)
        error=try
            A.mesh.generate(3);nothing
        catch failure
            failure isa InterruptException && rethrow();failure
        end
        @test error isa ArgumentError
        @test occursin("shared Quad",sprint(showerror,error))
    finally
        A.finalize()
    end
end

@testset "uniform-scale pair across disjoint Regions" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        two_boxes_scale!()
        A.mesh.set_periodic(2,[104],[6],SCALE2)
        A.mesh.generate(3)
        slave_link(104,SCALE2,6)
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "chained pair reuses a slave copy product as the next master" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # Face 6 is master of 104->6 and slave of 6->101: its produced copy
        # rows source the second copy across disjoint Regions. Gmsh 4.15.2
        # accepts the same chained declaration and publishes both links.
        # half maps Face 101 (z=0 of box2, x in [10,12], y in [0,2]) onto
        # Face 6 (x=0 of box1): (x,y,0) -> (0,(x-10)/2,y/2).
        half=[0.,0.,1.,0.,0.5,0.,0.,-5.,0.,0.5,0.,0.,0.,0.,0.,1.]
        two_boxes_scale!()
        A.mesh.set_periodic(2,[104],[6],SCALE2)
        A.mesh.set_periodic(2,[6],[101],half)
        A.mesh.generate(3)
        slave_link(104,SCALE2,6)
        slave_link(6,half,101)
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "adjacent-Face pairs sharing boundary points keep the declaration guard" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        multi_cube()
        # Faces 6 (x=0) and 3 (y=0) share a boundary edge — the rotation maps
        # each shared Point to itself, which upstream resolves but the
        # declaration's disjoint-tag contract still guards; chained pairs
        # over adjacent Faces remain a classified residual.
        rot=[0.,-1.,0.,0.,1.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.]
        @test_throws ArgumentError A.mesh.set_periodic(2,[6],[3],rot)
    finally
        A.finalize()
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
