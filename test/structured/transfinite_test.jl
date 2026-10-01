using Test
using Tessella
using Tessella.MeshTypes: Mesh, boundary_edges, mesh_crc, nnodes, node, nsegs,
                          ntris, triangle_area, validate

if !isdefined(Tessella,:Transfinite)
    Base.include(Tessella, joinpath(
        @__DIR__, "..", "..", "src", "structured", "Transfinite.jl"))
end
using Tessella.Transfinite: mesh_transfinite_patch

function _transfinite_rectangle(L::Int,H::Int;
                               xmin=0.0,xmax=Float64(L),
                               ymin=0.0,ymax=Float64(H))
    L>0&&H>0 || throw(ArgumentError("positive cell counts required"))
    bottom=[(xmin+(xmax-xmin)*i/L,ymin,0.0) for i in 0:L]
    right=[(xmax,ymin+(ymax-ymin)*j/H,0.0) for j in 0:H]
    top=[(xmax-(xmax-xmin)*i/L,ymax,0.0) for i in 0:L]
    left=[(xmin,ymax-(ymax-ymin)*j/H,0.0) for j in 0:H]
    return bottom,right,top,left
end

function _transfinite_surface_area(mesh)
    sum(triangle_area(node(mesh,mesh.tris[1,t]),node(mesh,mesh.tris[2,t]),
                      node(mesh,mesh.tris[3,t])) for t in 1:ntris(mesh);init=0.0)
end

function _transfinite_polygon_area(sides)
    ring=NTuple{3,Float64}[]
    for side in sides
        append!(ring,@view side[1:end-1])
    end
    area=0.0
    for i in eachindex(ring)
        p=ring[i];q=ring[mod1(i+1,length(ring))]
        area+=p[1]*q[2]-p[2]*q[1]
    end
    return abs(area)/2
end

function _transfinite_canonical_triangles(mesh)
    result=NTuple{3,Int32}[]
    for t in axes(mesh.tris,2)
        values=sort(mesh.tris[:,t])
        push!(result,(values[1],values[2],values[3]))
    end
    sort!(result)
end

function _expected_triangles(L,H,arrangement)
    id(i,j)=Int32(i+1+j*(L+1))
    result=NTuple{3,Int32}[]
    for i in 0:L-1,j in 0:H-1
        v1=id(i,j);v2=id(i+1,j);v3=id(i+1,j+1);v4=id(i,j+1)
        right=arrangement===:right ||
              (arrangement===:alternate_right&&isodd(i+j)) ||
              (arrangement===:alternate_left&&iseven(i+j))
        first,second=right ? ((v1,v2,v3),(v3,v4,v1)) :
                             ((v1,v2,v4),(v4,v2,v3))
        for triangle in (first,second)
            values=sort(collect(triangle))
            push!(result,(values[1],values[2],values[3]))
        end
    end
    sort!(result)
end

function _transfinite_edge_set(matrix)
    result=Set{NTuple{2,Int32}}()
    for i in axes(matrix,2)
        a=matrix[1,i];b=matrix[2,i]
        push!(result,a<b ? (a,b) : (b,a))
    end
    result
end

struct _CountOnlySide <: AbstractVector{NTuple{3,Float64}}
    count::Int
end
Base.size(side::_CountOnlySide)=(side.count,)
Base.getindex(::_CountOnlySide,::Int)=error("resource counts were not checked first")

@noinline function _transfinite_allocated(sides)
    GC.gc()
    return @allocated mesh_transfinite_patch(sides...;arrangement=:alternate_left)
end

@testset "four-sided planar transfinite patches" begin
    @testset "Gmsh arrangements, deterministic CRCs, and physical tags" begin
        sides=_transfinite_rectangle(3,2)
        expected_crc=Dict(
            :left=>"bc5c8cf1ba3bcbd22315c3095990658c9b95efe33ba8f7dcee99b9cc6dadf9e4",
            :right=>"b968747720c166a0a8f25a787528277a5265f8c187db70f7c73903354f05e151",
            :alternate_left=>"4265a3d7e815fb0510d77e6b6bda1b96f776c18bfa8b75fc03fa0a6a0df13f2d",
            :alternate_right=>"8c7c0070649022cbc5dcd81ae8abf0898d091e352ca8fb86b7e8a20e262dc8de")
        for arrangement in (:left,:right,:alternate_left,:alternate_right)
            mesh=mesh_transfinite_patch(sides...;arrangement=arrangement,
                                        face_tag=21,side_tags=(11,12,13,14))
            @test validate(mesh).ok
            @test (nnodes(mesh),nsegs(mesh),ntris(mesh))==(12,10,12)
            @test mesh_crc(mesh).bbox==((0.0,0.0,0.0),(3.0,2.0,0.0))
            @test mesh_crc(mesh).sha==expected_crc[arrangement]
            @test _transfinite_canonical_triangles(mesh)==_expected_triangles(3,2,arrangement)
            @test mesh.tri_tag==fill(Int32(21),12)
            @test mesh.seg_tag==Int32[11,11,11,12,12,13,13,13,14,14]
            boundary,maxincidence=boundary_edges(mesh.tris)
            @test maxincidence==2
            @test Set(boundary)==_transfinite_edge_set(mesh.segs)
            @test _transfinite_surface_area(mesh)==6.0
            @test mesh_crc(mesh)==mesh_crc(mesh_transfinite_patch(
                sides...;arrangement=arrangement,face_tag=21,
                side_tags=(11,12,13,14)))
        end
    end

    @testset "average-chord Coons interpolation and boundary conservation" begin
        bottom=[(0.,0.,0.),(0.7,-0.25,0.),(2.1,-0.55,0.),(3.2,-0.2,0.),(4.,0.,0.)]
        right=[(4.,0.,0.),(4.35,0.8,0.),(4.2,2.1,0.),(4.,3.,0.)]
        top=[(4.,3.,0.),(3.1,3.55,0.),(2.0,3.35,0.),(0.8,3.15,0.),(0.,3.,0.)]
        left=[(0.,3.,0.),(-0.3,2.15,0.),(-0.5,0.9,0.),(0.,0.,0.)]
        sides=(bottom,right,top,left)
        mesh=mesh_transfinite_patch(sides...;arrangement=:alternate_left)
        @test (nnodes(mesh),nsegs(mesh),ntris(mesh))==(20,14,24)
        @test validate(mesh).ok
        @test _transfinite_surface_area(mesh)≈_transfinite_polygon_area(sides) atol=64eps(Float64)

        # Boundary nodes are copied bit-for-bit, in deterministic row-major order.
        @test [node(mesh,i) for i in 1:5]==bottom
        @test [node(mesh,5+5j) for j in 0:3]==right
        @test [node(mesh,i+15) for i in 1:5]==reverse(top)
        @test [node(mesh,1+5j) for j in 0:3]==reverse(left)

        # Independent direct evaluation of the documented average-chord formula.
        top_grid=reverse(top);left_grid=reverse(left)
        chord(a,b)=hypot(a[1]-b[1],a[2]-b[2],a[3]-b[3])
        ustep=[0.5(chord(bottom[i+1],bottom[i])+chord(top_grid[i+1],top_grid[i]))
               for i in 1:4]
        vstep=[0.5(chord(right[j+1],right[j])+chord(left_grid[j+1],left_grid[j]))
               for j in 1:3]
        u=sum(ustep[1:2])/sum(ustep);v=vstep[1]/sum(vstep)
        c1=bottom[1];c2=bottom[end];c3=top_grid[end];c4=top_grid[1]
        expected=ntuple(3) do d
            (1-u)*left_grid[2][d]+u*right[2][d]+
            (1-v)*bottom[3][d]+v*top_grid[3][d]-
            ((1-u)*(1-v)*c1[d]+u*(1-v)*c2[d]+u*v*c3[d]+(1-u)*v*c4[d])
        end
        @test all(isapprox(node(mesh,8)[d],expected[d];
                           atol=16eps(Float64),rtol=16eps(Float64)) for d in 1:3)
    end

    @testset "tilted and clockwise patches" begin
        base=_transfinite_rectangle(4,3;xmax=2.0,ymax=1.5)
        ex=(inv(sqrt(2.0)),inv(sqrt(2.0)),0.0)
        ey=(-inv(sqrt(6.0)),inv(sqrt(6.0)),2inv(sqrt(6.0)))
        origin=(1.0,-2.0,3.0)
        transform(p)=(origin[1]+p[1]*ex[1]+p[2]*ey[1],
                      origin[2]+p[1]*ex[2]+p[2]*ey[2],
                      origin[3]+p[1]*ex[3]+p[2]*ey[3])
        tilted=map(side->transform.(side),base)
        mesh=mesh_transfinite_patch(tilted...;arrangement=:right)
        @test validate(mesh).ok
        @test _transfinite_surface_area(mesh)≈3.0 atol=256eps(Float64)

        clockwise=([(0.,0.,0.),(0.,1.,0.),(0.,2.,0.)],
                   [(0.,2.,0.),(1.,2.,0.),(2.,2.,0.),(3.,2.,0.)],
                   [(3.,2.,0.),(3.,1.,0.),(3.,0.,0.)],
                   [(3.,0.,0.),(2.,0.,0.),(1.,0.,0.),(0.,0.,0.)])
        reversed_mesh=mesh_transfinite_patch(clockwise...;
                                              arrangement=:alternate_right)
        @test validate(reversed_mesh).ok
        @test _transfinite_surface_area(reversed_mesh)==6.0

        long_patch=mesh_transfinite_patch(_transfinite_rectangle(2048,1)...)
        @test (nnodes(long_patch),ntris(long_patch))==(4098,4096)
        @test validate(long_patch).ok
    end

    @testset "validated blockers and resource limits" begin
        sides=_transfinite_rectangle(3,2)
        @test_throws ArgumentError mesh_transfinite_patch(
            [(0.,0.,0.)],sides[2],sides[3],sides[4])
        @test_throws ArgumentError mesh_transfinite_patch(
            sides[1],sides[2],[(3.,2.,0.),(0.,2.,0.)],sides[4])
        mismatched=copy(sides[2]);insert!(mismatched,2,(3.,0.5,0.))
        @test_throws ArgumentError mesh_transfinite_patch(
            sides[1],mismatched,sides[3],sides[4])
        broken=copy(sides[2]);broken[1]=(3.,nextfloat(0.0),0.)
        @test_throws ArgumentError mesh_transfinite_patch(
            sides[1],broken,sides[3],sides[4])
        nonfinite=copy(sides[1]);nonfinite[2]=(NaN,0.,0.)
        @test_throws ArgumentError mesh_transfinite_patch(
            nonfinite,sides[2],sides[3],sides[4])
        extra=Any[(point...,point==(1.,0.,0.) ? NaN : 0.0)
                  for point in sides[1]]
        @test_throws ArgumentError mesh_transfinite_patch(
            extra,sides[2],sides[3],sides[4])
        duplicate=copy(sides[1]);duplicate[2]=duplicate[1]
        @test_throws ArgumentError mesh_transfinite_patch(
            duplicate,sides[2],sides[3],sides[4])
        nonplanar=copy(sides[1]);nonplanar[2]=(1.,0.,1e-6)
        @test_throws ArgumentError mesh_transfinite_patch(
            nonplanar,sides[2],sides[3],sides[4])

        bowtie=([(0.,0.,0.),(1.,1.,0.)],
                [(1.,1.,0.),(0.,1.,0.)],
                [(0.,1.,0.),(1.,0.,0.)],
                [(1.,0.,0.),(0.,0.,0.)])
        @test_throws ArgumentError mesh_transfinite_patch(bowtie...)
        backtrack=([(0.,0.,0.),(1.,0.,0.),(0.5,0.,0.),(2.,0.,0.)],
                   [(2.,0.,0.),(2.,1.,0.)],
                   [(2.,1.,0.),(1.5,1.,0.),(1.,1.,0.),(0.,1.,0.)],
                   [(0.,1.,0.),(0.,0.,0.)])
        @test_throws ArgumentError mesh_transfinite_patch(backtrack...)
        folded=([(0.,0.,0.),(1.138676204796158,-0.8998793245683743,0.),
                 (1.,0.,0.)],
                [(1.,0.,0.),(0.9753676924206318,-0.6302141619964883,0.),
                 (1.,1.,0.)],
                [(1.,1.,0.),(-0.14986440109197496,1.9513598146977102,0.),
                 (0.,1.,0.)],
                [(0.,1.,0.),(-0.9877714621697844,-0.56016161339637,0.),
                 (0.,0.,0.)])
        @test_throws ArgumentError mesh_transfinite_patch(folded...)

        @test_throws ArgumentError mesh_transfinite_patch(sides...;arrangement=:alternate)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;arrangement="Left")
        @test_throws ArgumentError mesh_transfinite_patch(sides...;face_tag=true)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;face_tag=-1)
        @test_throws ArgumentError mesh_transfinite_patch(
            sides...;face_tag=big(typemax(Int32))+1)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;side_tags=(1,2,3))
        @test_throws ArgumentError mesh_transfinite_patch(sides...;side_tags=(1,2,true,4))
        @test_throws ArgumentError mesh_transfinite_patch(sides...;side_tags=(1,2,-1,4))
        @test_throws ArgumentError mesh_transfinite_patch(sides...;max_nodes=true)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;max_triangles=false)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;max_nodes=-1)
        @test_throws ArgumentError mesh_transfinite_patch(
            sides...;max_nodes=big(typemax(Int32))+1)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;max_nodes=12.0)
        boolean_side = Any[(true,0.,0.),(1.,0.,0.),(2.,0.,0.),(3.,0.,0.)]
        @test_throws ArgumentError mesh_transfinite_patch(
            boolean_side,sides[2],sides[3],sides[4])
        @test_throws ArgumentError mesh_transfinite_patch(sides...;max_nodes=11)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;max_triangles=11)
        bounded=mesh_transfinite_patch(sides...;max_nodes=12,max_triangles=12)
        @test (nnodes(bounded),ntris(bounded))==(12,12)

        huge=_CountOnlySide(typemax(Int))
        @test_throws ArgumentError mesh_transfinite_patch(huge,huge,huge,huge)
        wide=_CountOnlySide(100_000)
        @test_throws ArgumentError mesh_transfinite_patch(wide,wide,wide,wide)

        maximum=Float64(floatmax(Float64))
        overflowing=([(maximum,0.,0.),(-maximum,0.,0.)],
                     [(-maximum,0.,0.),(-maximum,1.,0.)],
                     [(-maximum,1.,0.),(maximum,1.,0.)],
                     [(maximum,1.,0.),(maximum,0.,0.)])
        @test_throws ArgumentError mesh_transfinite_patch(overflowing...)
    end

    @testset "allocation growth remains linear in output size" begin
        small_sides=_transfinite_rectangle(64,64)
        large_sides=_transfinite_rectangle(128,64)
        mesh_transfinite_patch(small_sides...;arrangement=:alternate_left)
        mesh_transfinite_patch(large_sides...;arrangement=:alternate_left)
        small=_transfinite_allocated(small_sides)
        large=_transfinite_allocated(large_sides)
        @test small>0
        @test large>small
        @test large<=2.25small+262_144
        @info "transfinite allocation ratchet" small_bytes=small large_bytes=large
    end

    @testset "warped transfinite patches (3-D Coons interpolation)" begin
        # z = 0.5*u*v bilinear warp: 25 nodes / 32 triangles / 16 segments.
        function _warped_quad(L::Int,H::Int;amplitude=0.5)
            bottom=[(Float64(i)/L,0.0,0.0) for i in 0:L]
            right=[(1.0,Float64(j)/H,amplitude*Float64(j)/H) for j in 0:H]
            top=[(1.0-Float64(i)/L,1.0,amplitude*(1.0-Float64(i)/L))
                 for i in 0:L]
            left=[(0.0,1.0-Float64(j)/H,0.0) for j in 0:H]
            return bottom,right,top,left
        end
        sides=_warped_quad(4,4)
        @test_throws ArgumentError mesh_transfinite_patch(sides...)
        mesh=mesh_transfinite_patch(sides...;arrangement=:left,
                                    allow_warped=true,face_tag=21,
                                    side_tags=(11,12,13,14))
        @test validate(mesh).ok
        @test (nnodes(mesh),nsegs(mesh),ntris(mesh))==(25,16,32)
        @test all(mesh.tri_tag.==21)
        @test (unique(mesh.seg_tag)|>sort)==[11,12,13,14]
        for i in 1:nnodes(mesh)
            x,y,z=node(mesh,i)
            @test z==0.5x*y
        end
        @test _transfinite_canonical_triangles(mesh)==
              sort!(_expected_triangles(4,4,:left))

        # Every arrangement preserves the warped interior exactly.
        for arrangement in (:right,:alternate_left,:alternate_right)
            candidate=mesh_transfinite_patch(
                sides...;arrangement=arrangement,allow_warped=true)
            @test validate(candidate).ok
            @test candidate.coords==mesh.coords
            @test _transfinite_canonical_triangles(candidate)==
                  sort!(_expected_triangles(4,4,arrangement))
        end

        # A coplanar boundary under `allow_warped` is bit-identical to the
        # planar path.
        planar=_transfinite_rectangle(4,4)
        plain=mesh_transfinite_patch(planar...;arrangement=:alternate_left)
        flagged=mesh_transfinite_patch(planar...;arrangement=:alternate_left,
                                      allow_warped=true)
        @test flagged.coords==plain.coords
        @test flagged.tris==plain.tris
        @test mesh_crc(flagged)==mesh_crc(plain)

        # The warped path allocates linearly, same growth envelope as planar.
        @noinline _warped_alloc(s)=@allocated mesh_transfinite_patch(
            s...;arrangement=:alternate_left,allow_warped=true)
        flat64=_transfinite_rectangle(64,64)
        warp_sides=_warped_quad(64,64)
        _warped_alloc(flat64);_warped_alloc(warp_sides)
        flat_alloc=_warped_alloc(flat64)
        warp_alloc=_warped_alloc(warp_sides)
        @test warp_alloc<=1.5flat_alloc+262_144

        # A genuinely non-coplanar ring that crosses itself in mid-air is
        # rejected by the 3-D boundary audit.
        c1=(0.0,0.0,0.0);c2=(1.0,1.0,1.0);c3=(0.0,1.0,0.0);c4=(1.0,0.0,1.0)
        lin(a,b,n)=[ntuple(d->a[d]+(b[d]-a[d])*i/(n-1),3) for i in 0:n-1]
        crossed=(lin(c1,c2,5),lin(c2,c3,5),lin(c3,c4,5),lin(c4,c1,5))
        @test_throws ArgumentError mesh_transfinite_patch(
            crossed...;allow_warped=true)

        # A nonadjacent vertex-vertex contact (the ring touches itself) is
        # also rejected: side 1's interior vertex lies on side 3's chain.
        t1=[(0.0,0.0,0.0),(1.0,0.0,1.0),(2.0,0.0,0.0)]
        t2=[(2.0,0.0,0.0),(2.0,1.0,0.0)]
        t3=[(2.0,1.0,0.0),(1.0,0.0,1.0),(0.0,1.0,0.0)]
        t4=[(0.0,1.0,0.0),(0.0,0.0,0.0)]
        @test_throws ArgumentError mesh_transfinite_patch(
            t1,t2,t3,t4;allow_warped=true)
    end

    @testset ".geo ruled-surface transfinite patches" begin
        warped_geo=mktemp() do path,io
            write(io,"""
Point(1) = {0,0,0};
Point(2) = {1,0,0};
Point(3) = {1,1,0.5};
Point(4) = {0,1,0};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Surface(1) = {1};
Transfinite Curve {1,2,3,4} = 5;
Transfinite Surface {1};
""")
            close(io)
            execute_geo(path;mesh_dim=2)
        end
        mesh=warped_geo.mesh
        @test validate(mesh).ok
        @test (size(mesh.coords,2),size(mesh.tris,2))==(25,32)
        for i in axes(mesh.coords,2)
            x,y,z=mesh.coords[:,i]
            @test z==0.5x*y
        end
        # Deterministic output across runs.
        rerun=mktemp() do path,io
            write(io,"""
Point(1) = {0,0,0}; Point(2) = {1,0,0}; Point(3) = {1,1,0.5}; Point(4) = {0,1,0};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Surface(1) = {1};
Transfinite Curve {1,2,3,4} = 5;
Transfinite Surface {1};
""")
            close(io)
            execute_geo(path;mesh_dim=2)
        end
        @test mesh_crc(rerun.mesh)==mesh_crc(mesh)

        # A `Plane Surface` whose boundary vertices are not coplanar meshes
        # anyway, matching upstream `planeSurface`: the declared plane comes
        # from the first non-collinear on-curve boundary samples — here the
        # plane through the samples of edges 1 and 2, z = 0.5y — boundary
        # nodes keep their true positions while the interior interpolates
        # exactly on the plane.
        tilted=mktemp() do path,io
            write(io,"""
Point(1) = {0,0,0}; Point(2) = {1,0,0}; Point(3) = {1,1,0.5}; Point(4) = {0,1,0};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Transfinite Curve {1,2,3,4} = 5;
Transfinite Surface {1};
""")
            close(io)
            execute_geo(path;mesh_dim=2)
        end
        @test validate(tilted.mesh).ok
        @test (size(tilted.mesh.coords,2),size(tilted.mesh.tris,2))==(25,32)
        for i in axes(tilted.mesh.coords,2)
            x,y,z=tilted.mesh.coords[:,i]
            interior=1e-9<x<1-1e-9 && 1e-9<y<1-1e-9
            interior && @test z ≈ 0.5y atol=1e-12
        end
        # The lifted corner keeps its true off-plane coordinate.
        @test any(i->tilted.mesh.coords[:,i] ≈ [1.0,1.0,0.5],
                  axes(tilted.mesh.coords,2))

        # A `Plane Surface` whose boundary curve leaves the declared plane
        # takes the same projected route: the arc's interior nodes keep
        # their true off-plane coordinates while the patch interior stays
        # exactly on the corner-defined y=0 plane.
        curved=mktemp() do path,io
            write(io,"""
Point(1) = {0,0,0}; Point(2) = {1,0,0};
Point(5) = {0,0,1}; Point(6) = {1,0,1}; Point(9) = {0.5,0,1};
Line(1) = {1,2}; Line(5) = {2,6}; Circle(6) = {6,9,5}; Line(7) = {5,1};
Curve Loop(1) = {1,5,6,7};
Plane Surface(1) = {1};
Transfinite Curve{1,5,6,7} = 5;
Transfinite Surface{1};
""")
            close(io)
            execute_geo(path;mesh_dim=2)
        end
        @test validate(curved.mesh).ok
        @test (size(curved.mesh.coords,2),size(curved.mesh.tris,2))==(25,32)
        # The arc bulges off the plane: its interior nodes keep y > 0, while
        # every patch-interior node lies exactly on y=0.
        arc_bulge=count(i->curved.mesh.coords[2,i]>1e-9,
                        axes(curved.mesh.coords,2))
        @test arc_bulge==3
        for i in axes(curved.mesh.coords,2)
            x,y,z=curved.mesh.coords[:,i]
            (abs(y)>1e-9 || x<1e-9 || x>1-1e-9 || z<1e-9) && continue
            @test y==0.0
        end
        # A coplanar `Surface` (ruled kind) now meshes identically to a plane
        # patch rather than rejecting the kind outright.
        flat=mktemp() do path,io
            write(io,"""
Point(1) = {0,0,0}; Point(2) = {1,0,0}; Point(3) = {1,1,0}; Point(4) = {0,1,0};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Surface(1) = {1};
Transfinite Curve {1,2,3,4} = 5;
Transfinite Surface {1};
""")
            close(io)
            execute_geo(path;mesh_dim=2)
        end
        @test (size(flat.mesh.coords,2),size(flat.mesh.tris,2))==(25,32)

        # `Surface ... In Sphere` keeps its own parameterization and does not
        # take the warped Coons path.
        sphere_err=try
            mktemp() do path,io
                write(io,"""
Point(1) = {0,0,0}; Point(2) = {1,0,0}; Point(3) = {1,1,0.5}; Point(4) = {0,1,0};
Point(5) = {0.5,0.5,-2};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Surface(1) = {1} In Sphere{5};
Transfinite Curve {1,2,3,4} = 5;
Transfinite Surface {1};
""")
                close(io)
                execute_geo(path;mesh_dim=2)
            end
            nothing
        catch caught
            caught
        end
        @test sphere_err isa ArgumentError
    end
end
