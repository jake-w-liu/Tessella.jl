module APIOCCCurvedSizedFillChecks
using Test, Tessella
const A=Tessella.API

# Curved OCC faces (Sphere/Cylinder/Cone/Torus) used to be filled by a fixed
# (u,v) lattice seeded from the boundary spacing: the interior ignored the
# background size (a selected size field was rejected outright), sphere poles
# collected 8-degree slivers and the cell count was off by -34%..+27% against
# Gmsh 4.15.2. They are now refined in the parameter plane under the surface
# metric by the same Frontal-Delaunay kernel as planes. Reference cell counts
# below are Gmsh 4.15.2 (`gmsh.model.mesh.generate(2)`, default options).

# Closed-shell statistics of the current dim-2 mesh: cell count, smallest
# triangle angle, mean edge length over `size_of(midpoint)`, whether every
# edge is shared by exactly two consistently oriented triangles, and the
# largest `distance(node)` residual.
function shell_statistics(size_of,distance)
    tags,xyz,_=A.mesh.get_nodes()
    row=Dict(tag=>index for (index,tag) in enumerate(tags))
    point(tag)=(xyz[3row[tag]-2],xyz[3row[tag]-1],xyz[3row[tag]])
    directed=Dict{Tuple{UInt64,UInt64},Int}()
    cells=0;smallest=180.0
    for (_,surface) in A.model.get_entities(2)
        types,_,nodes=A.mesh.get_elements(2,surface)
        for (type,connectivity) in zip(types,nodes)
            type==2 || continue
            for k in 1:3:length(connectivity)
                a,b,c=connectivity[k],connectivity[k+1],connectivity[k+2]
                pa,pb,pc=point(a),point(b),point(c)
                for (u,v,w) in ((pa,pb,pc),(pb,pc,pa),(pc,pa,pb))
                    e1=v.-u;e2=w.-u
                    cosine=sum(e1.*e2)/sqrt(sum(abs2,e1)*sum(abs2,e2))
                    smallest=min(smallest,acosd(clamp(cosine,-1,1)))
                end
                for edge in ((a,b),(b,c),(c,a))
                    directed[edge]=get(directed,edge,0)+1
                end
                cells+=1
            end
        end
    end
    closed=all(==(1),values(directed)) &&
        all(edge->haskey(directed,(edge[2],edge[1])),keys(directed))
    ratio=0.0;edges=0
    for (a,b) in keys(directed)
        a<b || continue
        pa,pb=point(a),point(b)
        middle=(pa.+pb)./2
        ratio+=sqrt(sum(abs2,pb.-pa))/size_of(middle...)
        edges+=1
    end
    residual=maximum(distance(point(tag)...) for tag in tags)
    return (cells=cells,smallest=smallest,ratio=ratio/edges,closed=closed,
            residual=residual,point=point,directed=directed,tags=tags)
end

function generate(build;h=nothing,field=nothing)
    A.initialize()
    try
        build();A.model.occ.synchronize()
        if h!==nothing
            A.option("Mesh.MeshSizeMin",h);A.option("Mesh.MeshSizeMax",h)
        end
        if field!==nothing
            tag=A.mesh.field.add("MathEval")
            A.mesh.field.set_string(tag,"F",field)
            A.mesh.field.set_as_background_mesh(tag)
            A.option("Mesh.MeshSizeExtendFromBoundary",0)
            A.option("Mesh.MeshSizeFromPoints",0)
            A.option("Mesh.MeshSizeFromCurvature",0)
        end
        A.mesh.generate(2)
        return true
    catch
        A.finalize()
        rethrow()
    end
end

@testset "sized Frontal-Delaunay fill of curved OCC faces" begin
    @testset "sphere follows the uniform size, poles included" begin
        generate(()->A.model.occ.add_sphere(0,0,0,1);h=0.3)
        s=try
            shell_statistics((x,y,z)->0.3,(x,y,z)->abs(hypot(x,y,z)-1))
        finally
            A.finalize()
        end
        @test s.closed
        @test s.residual<=1e-12
        @test 340<=s.cells<=420          # Gmsh: 380; the lattice gave 440
        # The lattice gave 8 deg at the poles; Gmsh's own minimum at this
        # size is 25.9 deg, so the bound separates slivers, not one degree.
        @test s.smallest>20
        @test 0.9<=s.ratio<=1.06
        # Each pole is an ordinary vertex of the closed mesh.
        for tag in s.tags
            abs(abs(s.point(tag)[3])-1)<=1e-12 || continue
            valence=count(edge->edge[1]==tag,keys(s.directed))
            @test 5<=valence<=7
        end
    end

    @testset "sphere pole fans at ordinary sizes (h=$h)" for h in (0.2,0.15,0.1)
        # The wrapped-cell rule added for sizes larger than the model once
        # also acted on the open sector of a growing pole fan and left two
        # slivers at the south pole: minimum angle 6.8, 6.3 and 9.6 deg at
        # these sizes (Gmsh 4.15.2: 29.0, 27.8, 26.6). It now runs only as a
        # second pass when wrapped cells are left, so a resolved sphere
        # never sees it.
        generate(()->A.model.occ.add_sphere(0,0,0,1);h=h)
        s=try
            shell_statistics((x,y,z)->h,(x,y,z)->abs(hypot(x,y,z)-1))
        finally
            A.finalize()
        end
        @test s.closed
        @test s.smallest>20
        for tag in s.tags
            abs(abs(s.point(tag)[3])-1)<=1e-12 || continue
            @test 4<=count(edge->edge[1]==tag,keys(s.directed))<=7
        end
    end

    @testset "torus (doubly periodic chart)" begin
        generate(()->A.model.occ.add_torus(0,0,0,1.5,0.5);h=0.3)
        s=try
            shell_statistics((x,y,z)->0.3,
                (x,y,z)->abs(hypot(hypot(x,y)-1.5,z)-0.5))
        finally
            A.finalize()
        end
        @test s.closed
        @test s.residual<=1e-12
        @test 740<=s.cells<=910          # Gmsh: 824; the lattice gave 560
        @test s.smallest>30
        @test 0.9<=s.ratio<=1.06
    end

    @testset "cone with apex and its base disk" begin
        generate(()->A.model.occ.add_cone(0,0,0,0,0,1,1,0);h=0.25)
        s=try
            shell_statistics((x,y,z)->0.25,(x,y,z)->
                min(abs(z),abs(hypot(x,y)-(1-z))/sqrt(2)))
        finally
            A.finalize()
        end
        @test s.closed
        @test s.residual<=1e-12
        @test 300<=s.cells<=385          # Gmsh: 342
        @test s.smallest>28
    end

    @testset "background size field on a sphere" begin
        # Previously: ArgumentError "a background size field on curved
        # Surface[1] is not implemented".
        size_of=(x,y,z)->0.15+0.35*abs(z)
        generate(()->A.model.occ.add_sphere(0,0,0,1);
                 field="0.15+0.35*Fabs(z)")
        s=try
            shell_statistics(size_of,(x,y,z)->abs(hypot(x,y,z)-1))
        finally
            A.finalize()
        end
        @test s.closed
        @test s.residual<=1e-12
        @test 385<=s.cells<=500          # Gmsh: 442
        @test 0.85<=s.ratio<=1.1
        @test s.smallest>22
    end
end
end
