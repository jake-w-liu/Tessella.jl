module Mesh2DRefineIncrementalChecks
using Test, Tessella
const K=Tessella.Mesh2D

# Two passes of the Ruppert loop used to cost time proportional to the whole
# triangulation on every step: each subsegment split reclassified every cell
# as inside or outside by flood fill (67% of a refinement to 150,000
# triangles, and quadratic: 747,000 triangles took 70 s), and each step
# scanned every subsegment for an encroached one (21%). The split now flags
# only the new vertex's star, and the search examines only the subsegments
# whose adjacent cells changed. Both must agree with the full computation at
# every step.
@testset "refinement keeps interior flags and encroachment candidates incrementally" begin
    # A square with a square hole and one embedded segment, so boundary
    # constraints (which change side) and internal ones (which do not) both
    # get split.
    xs=[0.0,1.0,1.0,0.0,0.4,0.6,0.6,0.4,0.1,0.9]
    ys=[0.0,0.0,1.0,1.0,0.4,0.4,0.6,0.6,0.8,0.9]
    T=K.constrained_delaunay(xs,ys,
        [(1,2),(2,3),(3,4),(4,1),(5,8),(8,7),(7,6),(6,5)];
        internal_segments=[(9,10)])
    interior=Vector{Bool}(collect(K.classify_interior(T)))
    pointids=Dict{NTuple{2,Float64},Int32}()
    for i in 1:T.nreal
        pointids[(T.x[i],T.y[i])]=Int32(i)
    end
    agrees()=begin
        full=K._classify_interior(T)
        all(t->!T.alive[t] || interior[t]==full[t],eachindex(T.alive))
    end
    B=1/(2*sind(25.0));area=0.004
    T.bad_scan_from=Int32(1);K._reset_bad_pending!(T)
    candidates=K._enc_candidates(T)
    cells=Int32[]
    steps=0;splits=0;internal_splits=0
    while steps<3000
        steps+=1
        enc=K._find_encroached!(T,candidates)
        @test enc==K._find_encroached(T)
        if enc!==nothing
            internal_splits+=K._skey(enc...) in T.internal
            K._split_subsegment!(T,enc[1],enc[2],interior,pointids,candidates,
                nothing,cells)
            splits+=1
            @test agrees()
            continue
        end
        bad=K._find_bad(T,interior,B,area,nothing,nothing,nothing)
        bad==0 && break
        a,b,c=K._vert(T,bad,1),K._vert(T,bad,2),K._vert(T,bad,3)
        cc=K._circumcenter2(K._pt(T,a),K._pt(T,b),K._pt(T,c))
        sp=K._encroached_by_point(T,cc)
        if sp!==nothing
            K._split_subsegment!(T,sp[1],sp[2],interior,pointids,candidates,
                nothing,cells)
            splits+=1
            @test agrees()
            continue
        end
        pa,pb,pc=K._pt(T,a),K._pt(T,b),K._pt(T,c)
        K._insert_steiner!(T,cc,interior,pointids;
            fallback=(pa[1]/3+pb[1]/3+pc[1]/3,pa[2]/3+pb[2]/3+pc[2]/3),
            encpending=candidates,newt=cells)
    end
    @test steps<3000
    @test splits>=8 && internal_splits>0
    @test agrees()
    @test K._find_encroached(T)===nothing
    # The hole stays a hole.
    mesh=K.to_mesh(T;interior)
    covered=sum(abs((mesh.coords[1,mesh.tris[2,t]]-mesh.coords[1,mesh.tris[1,t]])*
                    (mesh.coords[2,mesh.tris[3,t]]-mesh.coords[2,mesh.tris[1,t]])-
                    (mesh.coords[1,mesh.tris[3,t]]-mesh.coords[1,mesh.tris[1,t]])*
                    (mesh.coords[2,mesh.tris[2,t]]-mesh.coords[2,mesh.tris[1,t]]))/2
                for t in axes(mesh.tris,2))
    @test covered≈0.96 atol=1e-12

    # The whole refinement, compiled: under 1 kB allocated per triangle (it
    # was 2.2 kB at this size and grew with it, the flood fill's bit vectors).
    # The suite runs interpreted, where every call allocates.
    if Base.JLOptions().compile_enabled==1
        square=(xs[1:8],ys[1:8],[(1,2),(2,3),(3,4),(4,1),(5,8),(8,7),(7,6),(6,5)])
        Tessella.mesh_planar(square...;max_area=1e-3)
        bytes=@allocated(result=Tessella.mesh_planar(square...;max_area=1e-4))
        @test bytes<=1000*size(result.tris,2)
    end
end
end
