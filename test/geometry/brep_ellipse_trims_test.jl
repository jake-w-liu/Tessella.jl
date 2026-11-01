module EllipseTrimChecks
using Test,Tessella,LinearAlgebra
using Tessella.BRep: parse_step_entities
using Tessella.MeshTypes: nnodes,ntris,ntets,tet_volume,validate
const E=Tessella.Elements
const BRep=Tessella.BRep
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function _signed_surface_volume(mesh)
    volume=0.0
    for k in 1:size(mesh.tris,2)
        a,b,c=(mesh.coords[:,mesh.tris[i,k]] for i in 1:3)
        volume+=dot(a,cross(b,c))/6
    end
    return volume
end

_tet_volumes(mesh)=[det(hcat((mesh.coords[:,cell[i]]-mesh.coords[:,cell[1]]
                    for i in 2:4)...))/6 for cell in eachcol(mesh.tets)]

# Every undirected surface edge must be bounded by exactly two triangles,
# once in each direction — the shared-rim welds included.
function _edge_incidence(mesh)
    counts=Dict{NTuple{2,Int32},Int}()
    for t in eachcol(mesh.tris),(a,b) in
            ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
        key=a<b ? (a,b) : (b,a)
        counts[key]=get(counts,key,0)+1
    end
    return counts
end

function _with_source(body,fixture,replacements...)
    source=read(joinpath(FIXTURES,fixture),String)
    for (old,new) in replacements
        occursin(old,source) || throw(ArgumentError("fixture rewrite anchor missing: $old"))
        source=replace(source,old=>new)
    end
    return mktemp() do path,io
        write(io,source);close(io);body(path)
    end
end

function _entity_kinds(fixture)
    ents=parse_step_entities(read(joinpath(FIXTURES,fixture),String))
    kinds=Dict{String,Int}()
    for e in values(ents); kinds[e.kind]=get(kinds,e.kind,0)+1; end
    return ents,kinds
end

# n-gon polygonal area under an ellipse arc sampled at m parametric
# segments in [ta,tb] — the independent analytic area certificate for the
# bounded ellipse grid.
function _ellpoly_area(r1,r2,ta,tb,m)
    pts=[(r1*cos(ta+(tb-ta)*k/m),r2*sin(ta+(tb-ta)*k/m)) for k in 0:m]
    a=0.0
    for k in 1:m
        a+=pts[k][1]*pts[k+1][2]-pts[k+1][1]*pts[k][2]
    end
    return 0.5*abs(a)
end

@testset "ellipse/partial/rational entity audit" begin
    _,kinds=_entity_kinds("brep_occ_elliptical_pipe.step")
    @test kinds["MANIFOLD_SOLID_BREP"]==1 && kinds["CLOSED_SHELL"]==1
    @test kinds["ADVANCED_FACE"]==3
    @test kinds["ELLIPSE"]==2 && kinds["PLANE"]==2
    @test kinds["SURFACE_OF_LINEAR_EXTRUSION"]==1
    @test kinds["EDGE_CURVE"]==3 && kinds["LINE"]==1

    _,kinds=_entity_kinds("brep_occ_ellipse_arc_plate.step")
    @test kinds["ADVANCED_FACE"]==6
    @test kinds["ELLIPSE"]==2 && kinds["PLANE"]==5
    @test kinds["SURFACE_OF_LINEAR_EXTRUSION"]==1

    _,kinds=_entity_kinds("brep_occ_ellipse_hole_plate.step")
    @test kinds["ADVANCED_FACE"]==7
    @test kinds["ELLIPSE"]==2 && kinds["FACE_BOUND"]==2
    @test kinds["SURFACE_OF_LINEAR_EXTRUSION"]==1

    _,kinds=_entity_kinds("brep_occ_spline_partial_plate.step")
    @test kinds["ADVANCED_FACE"]==6 && kinds["B_SPLINE_CURVE_WITH_KNOTS"]==2
    @test kinds["SURFACE_OF_LINEAR_EXTRUSION"]==1

    _,kinds=_entity_kinds("brep_occ_prism_multi_rim.step")
    @test kinds["ADVANCED_FACE"]==6 && kinds["B_SPLINE_CURVE_WITH_KNOTS"]==2
    @test kinds["SURFACE_OF_LINEAR_EXTRUSION"]==1
end

@testset "elliptical pipe: full ELLIPSE trims weld a swept wall" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_elliptical_pipe.step");
                        fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test size(surface.tets,2)==0
    @test nnodes(surface)==96 && ntris(surface)==188
    @test all(values(_edge_incidence(surface)).==(2))
    @test vec(minimum(surface.coords,dims=2))≈[-2.,-1.,0.] atol=1e-12
    @test vec(maximum(surface.coords,dims=2))≈[2.,1.,1.] atol=1e-12
    # Every node sits on the r1=2, r2=1 ellipse about (0,0) at z in {0,1}.
    for i in axes(surface.coords,2)
        p=surface.coords[:,i]
        @test isapprox((p[1]/2)^2+p[2]^2,1.0;atol=1e-9)
        @test p[3]==0.0 || p[3]==1.0
    end
    # The 48-gon ellipse area, sampled parametrically, is the exact volume.
    @test _signed_surface_volume(surface)≈_ellpoly_area(2.0,1.0,0.0,2π,48) atol=1e-9

    filled=import_step(joinpath(FIXTURES,"brep_occ_elliptical_pipe.step"))
    @test filled isa E.Mesh && validate(filled).ok && ntets(filled)>0
    volumes=_tet_volumes(filled)
    @test all(>(0),volumes)
    @test sum(volumes)≈_ellpoly_area(2.0,1.0,0.0,2π,48) atol=1e-9
end

@testset "ellipse-arc plate: open ELLIPSE trims bound a swept wall" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_ellipse_arc_plate.step");
                        fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test nnodes(surface)==54 && ntris(surface)==104
    @test all(values(_edge_incidence(surface)).==(2))
    @test vec(minimum(surface.coords,dims=2))≈[0.,0.,0.] atol=1e-12
    @test vec(maximum(surface.coords,dims=2))≈[4.5,2.,1.] atol=1e-12
    # The half ellipse centre (3,1), r1=1.5 along x, r2=1 along y, spans
    # theta in [-pi/2,pi/2]: 24 sampled segments on the canonical grid.
    wallnodes=[surface.coords[:,i] for i in axes(surface.coords,2) if
        isapprox(((surface.coords[1,i]-3)/1.5)^2+(surface.coords[2,i]-1)^2,
                 1.0;atol=1e-9)]
    @test length(wallnodes)==25*2
    @test all(p->p[1]>=3-1e-12,wallnodes)
    @test all(p->p[3]==0.0 || p[3]==1.0,wallnodes)
    # Plate area = 3x2 rectangle plus the 24-gon half-ellipse.
    @test _signed_surface_volume(surface)≈
        6.0+_ellpoly_area(1.5,1.0,-π/2,π/2,24) atol=1e-9

    filled=import_step(joinpath(FIXTURES,"brep_occ_ellipse_arc_plate.step"))
    @test filled isa E.Mesh && validate(filled).ok && ntets(filled)>0
    @test all(>(0),_tet_volumes(filled))
end

@testset "ellipse-hole plate: closed ELLIPSE bounds a swept hole wall" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_ellipse_hole_plate.step");
                        fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test nnodes(surface)==104 && ntris(surface)==208
    @test all(values(_edge_incidence(surface)).==(2))
    @test vec(minimum(surface.coords,dims=2))≈[0.,0.,0.] atol=1e-12
    @test vec(maximum(surface.coords,dims=2))≈[5.,4.,1.] atol=1e-12
    # Hole wall nodes lie on the r1=1, r2=0.5 ellipse about (2.5,2).
    wallnodes=[surface.coords[:,i] for i in axes(surface.coords,2) if
        isapprox((surface.coords[1,i]-2.5)^2+((surface.coords[2,i]-2)/0.5)^2,
                 1.0;atol=1e-9)]
    @test length(wallnodes)==48*2
    # 5x4 plate minus the 48-gon elliptical hole.
    @test _signed_surface_volume(surface)≈
        20.0-_ellpoly_area(1.0,0.5,0.0,2π,48) atol=1e-9

    filled=import_step(joinpath(FIXTURES,"brep_occ_ellipse_hole_plate.step"))
    @test filled isa E.Mesh && validate(filled).ok && ntets(filled)>0
    @test all(>(0),_tet_volumes(filled))
end

@testset "partial-domain B-spline trims sample only the active interval" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_spline_partial_plate.step");
                        fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test all(values(_edge_incidence(surface)).==(2))
    @test vec(minimum(surface.coords,dims=2))≈[0.,0.,0.] atol=1e-12
    @test vec(maximum(surface.coords,dims=2))≈[4.,2.,1.] atol=1e-12
    # The sampled chain covers only [0.25,0.75] of the quadratic Bézier:
    # 8 segments on the overlapped single span, endpoints on the chord.
    bulge=[surface.coords[:,i] for i in axes(surface.coords,2) if
        surface.coords[1,i]>3.75+1e-9]
    @test length(bulge)==7*2   # interior grid points at both heights
    @test all(p->p[3]==0.0 || p[3]==1.0,bulge)
    # The partial window of x(u)=3+4u-4u^2 stays inside [3.75,4].
    @test all(p->3.75-1e-9<=p[1]<=4.0+1e-9,bulge)
    # Volume: trapezoid minus the covered corner — bounded independently by
    # the polygon area of the exact Bézier window sampled at the same grid.
    quadarea=0.5*(3.75*0.5+3.75*1.5+3.75*2)  # shoelace of the chord quad
    @test 5.0<_signed_surface_volume(surface)<quadarea+1.0

    # A partial use whose vertex is off the curve still rejects loudly:
    # sliding the shared vertex along the a->p1 chord keeps the wall plane
    # valid but lands the vertex beside the spline — inversion fails and
    # the edge defers to the curve-kind diagnostic.
    _with_source("brep_occ_spline_partial_plate.step",
        "#2 = CARTESIAN_POINT('',(3.75,0.5,0.));"=>
        "#2 = CARTESIAN_POINT('',(3.8,0.5066666666666667,0.));") do path
        err=try import_step(path;fill=false);nothing catch e e end
        @test err isa ArgumentError
        @test occursin("B_SPLINE_CURVE_WITH_KNOTS",sprint(showerror,err))
    end
end

@testset "multi-chain swept rims weld N congruent chains" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_prism_multi_rim.step");
                        fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test nnodes(surface)==38 && ntris(surface)==72
    @test all(values(_edge_incidence(surface)).==(2))
    @test vec(minimum(surface.coords,dims=2))≈[0.,0.,0.] atol=1e-12
    # The degree-2 profile approximates its 4.6-x control points; the
    # sampled chain's extreme is the grid point x(0.75)=x(1.25)=4.1625.
    @test vec(maximum(surface.coords,dims=2))≈[4.1625,2.,1.] atol=1e-12
    # The two rim pieces cover [0,1] and [1,2] of the shared basis curve;
    # every sampled point lies on the extrusion profile at both heights.
    top=[surface.coords[:,i] for i in axes(surface.coords,2) if
         surface.coords[3,i]≈1.0 && surface.coords[1,i]>3.0+1e-9]
    @test length(top)==15   # 17-node rim chain minus both x=3 endpoints
    filled=import_step(joinpath(FIXTURES,"brep_occ_prism_multi_rim.step"))
    @test filled isa E.Mesh && validate(filled).ok && ntets(filled)>0
    @test all(>(0),_tet_volumes(filled))
end

@testset "rational weights: exact circle through RATIONAL_B_SPLINE_CURVE" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_rational_bulge_plate.step");
                        fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test nnodes(surface)==38 && ntris(surface)==72
    @test all(values(_edge_incidence(surface)).==(2))
    # The semicircular bulge centre (3,1) r=1 is EXACT via the rational
    # quadratic — every sampled rim node satisfies the circle equation.
    rim=[surface.coords[:,i] for i in axes(surface.coords,2) if
         isapprox((surface.coords[1,i]-3)^2+(surface.coords[2,i]-1)^2,
                  1.0;atol=1e-9)]
    @test length(rim)==17*2   # 16 sampled segments on the two-span domain
    @test all(p->p[1]>=3-1e-12,rim)
    @test all(p->p[3]==0.0 || p[3]==1.0,rim)
    # Bulge polygon area — independent closed form: the sampled chain is
    # the rational-quadratic circle map evaluated on the u-grid 0:1/8:2
    # (uniform in u, not in angle — every node still lands on the circle).
    wq=sqrt(0.5)
    rq(p0,p1,p2,u)=((1-u)^2 .*p0 .+2wq*u*(1-u).*p1 .+u^2 .*p2) ./
                   ((1-u)^2+2wq*u*(1-u)+u^2)
    rim=[(3.,0.)]
    for u in (1:8)./8; push!(rim,rq((3.,0.),(4.,0.),(4.,1.),u)); end
    for u in (1:8)./8; push!(rim,rq((4.,1.),(4.,2.),(3.,2.),u)); end
    a=0.0
    for k in 1:length(rim)-1
        a+=rim[k][1]*rim[k+1][2]-rim[k+1][1]*rim[k][2]
    end
    a+=rim[end][1]*rim[1][2]-rim[1][1]*rim[end][2]
    @test _signed_surface_volume(surface)≈6.0+abs(a)/2 atol=1e-9

    # A non-unit weight carries the true circle; rewiring the weight list
    # to all-ones must change the certified shell — the weights are live.
    _with_source("brep_occ_rational_bulge_plate.step",
        "(1.,0.7071067811865476,1.,0.7071067811865476,1.)"=>
        "(1.,1.,1.,1.,1.)") do path
        s=import_step(path;fill=false)
        @test s isa E.Mesh && validate(s).ok
        @test _signed_surface_volume(s)!=_signed_surface_volume(surface)
    end
end

@testset "rational NURBS surface: weights lift the trimmed cap" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_nurbs_rational_bump.step");
                        fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test all(values(_edge_incidence(surface)).==(2))
    @test vec(minimum(surface.coords,dims=2))≈[0.,0.,0.] atol=1e-12
    ents,_=_entity_kinds("brep_occ_nurbs_rational_bump.step")
    sent=only([e for e in values(ents) if e.kind=="COMPLEX" &&
               BRep._step_nurbs_surface(e,ents)!==nothing])
    s=BRep._step_nurbs_surface(sent,ents)
    @assert s isa Tessella.NURBSSurface
    @test any(!=(1.0),s.weights)
    # Analytic apex on the parsed weighted patch: centre weight 0.5 pulls
    # z(0.5,0.5) to exactly 15/14 — a unit-weight patch would apex at 1.125.
    @test isapprox(collect(nurbs_eval(s,0.5,0.5)),[1.,0.5,15/14];atol=1e-12)
    # The metric interior lattice undershoots the apex but stays lifted and
    # strictly below the unit-weight profile at the same seeding.
    zmax=maximum(surface.coords[3,:])
    @test 1.06<zmax<15/14
    # Every lifted node satisfies the RATIONAL patch equation — verified
    # by inverting (x,y) on the parsed weighted surface, not a poly fit.
    u_lo,u_hi=s.knots_u[3],s.knots_u[end-2]
    v_lo,v_hi=s.knots_v[3],s.knots_v[end-2]
    for i in axes(surface.coords,2)
        p=surface.coords[:,i]
        p[3]>1+1e-9 || continue
        u,v=0.5,0.5
        for _ in 1:40
            S=nurbs_eval(s,u,v)
            h=1e-7
            Su=nurbs_eval(s,clamp(u+h,u_lo,u_hi),v)
            Sv=nurbs_eval(s,u,clamp(v+h,v_lo,v_hi))
            ru=(Su[1]-S[1])/h,(Su[2]-S[2])/h
            rv=(Sv[1]-S[1])/h,(Sv[2]-S[2])/h
            det=ru[1]*rv[2]-ru[2]*rv[1]
            abs(det)>0 || break
            r=(p[1]-S[1],p[2]-S[2])
            du=(r[1]*rv[2]-r[2]*rv[1])/det
            dv=(ru[1]*r[2]-ru[2]*r[1])/det
            u=clamp(u+du,u_lo,u_hi);v=clamp(v+dv,v_lo,v_hi)
        end
        q=nurbs_eval(s,u,v)
        @test hypot(q[1]-p[1],q[2]-p[2],q[3]-p[3])<1e-9
    end

    # The weights are live: rewiring the list to all ones must change the
    # certified shell.
    _with_source("brep_occ_nurbs_rational_bump.step",
        "(1.,1.,1.),(1.,0.5,1.),(1.,1.,1.)"=>
        "(1.,1.,1.),(1.,1.,1.),(1.,1.,1.)") do path
        s2=import_step(path;fill=false)
        @test s2 isa E.Mesh && validate(s2).ok
        @test _signed_surface_volume(s2)!=_signed_surface_volume(surface)
    end

    filled=import_step(joinpath(FIXTURES,"brep_occ_nurbs_rational_bump.step"))
    @test filled isa E.Mesh && validate(filled).ok && ntets(filled)>0
    @test all(>(0),_tet_volumes(filled))
end

end # module
