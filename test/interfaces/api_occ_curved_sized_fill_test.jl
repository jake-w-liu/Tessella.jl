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
    # One mesh per behaviour, each as small as still shows it: the suite runs
    # interpreted, and under coverage every kernel line is interpreted too.
    @testset "sphere follows the uniform size, poles included" begin
        # h=0.2 is also the pole-fan regression: the wrapped-cell rule added
        # for sizes larger than the model once acted on the open sector of a
        # growing pole fan and left two slivers at the south pole (minimum
        # angle 6.8 deg here; Gmsh 4.15.2: 29.0). It now runs only as a
        # second pass when wrapped cells are left. The fan itself then still
        # kept a cell whose far edge passes beside the pole (23.8 deg) until
        # such edges were exchanged; see the last testset.
        generate(()->A.model.occ.add_sphere(0,0,0,1);h=0.2)
        s=try
            shell_statistics((x,y,z)->0.2,(x,y,z)->abs(hypot(x,y,z)-1))
        finally
            A.finalize()
        end
        @test s.closed
        @test s.residual<=1e-12
        @test 710<=s.cells<=870          # Gmsh: 790
        @test s.smallest>30              # Gmsh: 29.0; was 23.8, slivers 6.8
        @test 0.9<=s.ratio<=1.06
        # Each pole is an ordinary vertex of the closed mesh.
        for tag in s.tags
            abs(abs(s.point(tag)[3])-1)<=1e-12 || continue
            @test 4<=count(edge->edge[1]==tag,keys(s.directed))<=7
        end
    end

    @testset "torus (doubly periodic chart)" begin
        generate(()->A.model.occ.add_torus(0,0,0,1,0.4);h=0.35)
        s=try
            shell_statistics((x,y,z)->0.35,
                (x,y,z)->abs(hypot(hypot(x,y)-1,z)-0.4))
        finally
            A.finalize()
        end
        @test s.closed
        @test s.residual<=1e-12
        @test 310<=s.cells<=400          # Gmsh: 354
        @test s.smallest>25
        @test 0.9<=s.ratio<=1.06
    end

    @testset "cone with apex and its base disk" begin
        generate(()->A.model.occ.add_cone(0,0,0,0,0,1,1,0);h=0.35)
        s=try
            shell_statistics((x,y,z)->0.35,(x,y,z)->
                min(abs(z),abs(hypot(x,y)-(1-z))/sqrt(2)))
        finally
            A.finalize()
        end
        @test s.closed
        @test s.residual<=1e-12
        @test 165<=s.cells<=215          # Gmsh: 190
        @test s.smallest>28
    end

    @testset "background size field on a sphere" begin
        # Previously: ArgumentError "a background size field on curved
        # Surface[1] is not implemented".
        size_of=(x,y,z)->0.25+0.4*abs(z)
        generate(()->A.model.occ.add_sphere(0,0,0,1);
                 field="0.25+0.4*Fabs(z)")
        s=try
            shell_statistics(size_of,(x,y,z)->abs(hypot(x,y,z)-1))
        finally
            A.finalize()
        end
        @test s.closed
        @test s.residual<=1e-12
        @test 180<=s.cells<=240          # Gmsh: 208
        @test 0.85<=s.ratio<=1.1
        @test s.smallest>20
    end

    @testset "edge passing beside a pole is exchanged" begin
        # A frontal candidate beyond a pole is outside the chart and rejected
        # (Gmsh does the same), which leaves a fan cell whose far edge passes
        # beside the pole: minimum angles of 3 to 18 deg on about one sphere
        # in nine over centre, radius and size, in Gmsh 4.15.2 as well. Here
        # the cell (pole,a,b) spans more than a half turn, so it lies folded
        # over the pole; d is the vertex across a-b.
        M=Tessella.Model
        A.initialize()
        geometry=try
            A.model.occ.add_sphere(0,0,0,1);A.model.occ.synchronize()
            A.CURRENT[].surface_geometry[1]
        finally
            A.finalize()
        end
        north=pi/2
        uv=[1.8 0.2 3.4 1.8; north north-0.25 north-0.25 north-0.6]
        world=Matrix{Float64}(undef,3,4)
        for k in 1:4
            world[:,k].=M._occ_surface_point(geometry,uv[1,k],uv[2,k])
        end
        owned=[false,true,true,true]
        none=Set{NTuple{2,NTuple{2,Float64}}}()
        cells()=Int32[1 3; 2 2; 3 4]
        tris=cells()
        @test M._occ_curved_surface_fan_flips!(uv,tris,world,owned,geometry,none)==1
        @test tris==Int32[1 1; 2 4; 4 3]
        # A constrained edge stays, and so does an edge that faces an
        # ordinary vertex.
        tris=cells()
        kept=Set([((uv[1,2],uv[2,2]),(uv[1,3],uv[2,3]))])
        @test M._occ_curved_surface_fan_flips!(uv,tris,world,owned,geometry,kept)==0
        @test tris==cells()
        @test M._occ_curved_surface_fan_flips!(uv,tris,world,trues(4),geometry,none)==0
        @test tris==cells()
    end

    @testset "background size of a curved face is a plain number" begin
        # The size is read inside the metric of every cavity test. Its general
        # evaluation goes through untyped field and callback slots, which
        # boxed the value at every step and every quantity computed from it
        # (about 150 of the 240 heap objects per triangle of a sphere mesh).
        # Without a field, a callback or curvature sizing it is a constant
        # and takes a typed path; both paths must agree.
        M=Tessella.Model
        A.initialize()
        try
            A.model.occ.add_sphere(0,0,0,1);A.model.occ.synchronize()
            A.option("Mesh.MeshSizeMin",0.0);A.option("Mesh.MeshSizeMax",0.2)
            m=A.CURRENT[]
            options=M._model_algorithm_size_options(m,nothing,"test")
            size=M._OCCSurfaceParamBackgroundSize(m,1,m.surface_geometry[1],
                options,nothing,"test")
            for (u,v) in ((0.3,0.2),(4.0,-1.1),(6.0,1.5))
                @test size(u,v)===M._occ_surface_param_size(size,u,v)
            end
            @test size(0.3,0.2) isa Float64
            if Base.JLOptions().compile_enabled==1
                size(1.0,0.5)
                @test (@allocated size(1.0,0.5))==0
            end
        finally
            A.finalize()
        end
    end
end
end
