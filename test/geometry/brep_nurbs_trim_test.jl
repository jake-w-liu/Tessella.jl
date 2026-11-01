module NurbsTrimChecks
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

# Every undirected surface edge must be bounded by exactly two triangles —
# an unwelded periodic seam column or an unwelded trim shows up as a
# boundary (count-1) or over-used (count-3+) edge here.
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

# Independent on-surface residual: project the node's plan position onto
# the degree-2 octagonal rim B-spline parsed straight from the STEP text,
# using only the public curve evaluator — the wall is ruled along +z, so a
# wall node must project to zero residual regardless of height.
function _rim_curve(fixture)
    ents,_=_entity_kinds(fixture)
    rims=[e for e in values(ents) if e.kind=="B_SPLINE_CURVE_WITH_KNOTS"]
    @assert !isempty(rims)
    curve=BRep._step_nurbs_curve(first(rims),ents)
    @assert curve isa Tessella.NURBSCurve
    return curve
end

function _curve_residual(curve,p)
    u_lo=curve.knots[curve.degree+1]; u_hi=curve.knots[end-curve.degree]
    u_hi>u_lo || error("degenerate rim domain")
    θ=atan(p[2],p[1]);θ<0 && (θ+=2π)
    # Angle-seeded Newton (the rim is near-arc-length uniform) plus a small
    # coarse safety net; keep the converged root nearest the target.
    seeds=Float64[u_lo+(u_hi-u_lo)*θ/(2π)]
    for k in 1:8
        push!(seeds,u_lo+(u_hi-u_lo)*k/9)
    end
    best=Inf
    for s in seeds
        u=s
        for _ in 1:40
            h=max(1e-7*(u_hi-u_lo),1e-9)
            ua=clamp(u-h,u_lo,u_hi); ub=clamp(u+h,u_lo,u_hi)
            c0=nurbs_eval(curve,u);ca=nurbs_eval(curve,ua);cb=nurbs_eval(curve,ub)
            du=((cb[1]-ca[1])/(ub-ua),(cb[2]-ca[2])/(ub-ua),(cb[3]-ca[3])/(ub-ua))
            r=(p[1]-c0[1],p[2]-c0[2],p[3]-c0[3])
            den=du[1]^2+du[2]^2+du[3]^2
            den>0 || break
            step=(r[1]*du[1]+r[2]*du[2]+r[3]*du[3])/den
            un=clamp(u+step,u_lo,u_hi)
            done=abs(un-u)<=1e-13*(u_hi-u_lo);u=un
            done && break
        end
        c=nurbs_eval(curve,u)
        d=hypot(c[1]-p[1],c[2]-p[2],c[3]-p[3])
        d<best && (best=d)
        best<1e-12 && break
    end
    return best
end

@testset "trimmed-NURBS entity audit" begin
    _,kinds=_entity_kinds("brep_occ_nurbs_bump_plate.step")
    @test kinds["MANIFOLD_SOLID_BREP"]==1 && kinds["CLOSED_SHELL"]==1
    @test kinds["ADVANCED_FACE"]==6
    @test kinds["PLANE"]==5 && kinds["B_SPLINE_SURFACE_WITH_KNOTS"]==1
    @test kinds["B_SPLINE_CURVE_WITH_KNOTS"]==4 && kinds["LINE"]==8
    @test kinds["EDGE_CURVE"]==12
    # The AP214 representation context keeps the shape transferable.
    @test kinds["SHAPE_DEFINITION_REPRESENTATION"]==1
    @test kinds["ADVANCED_BREP_SHAPE_REPRESENTATION"]==1

    _,kinds=_entity_kinds("brep_occ_nurbs_pipe.step")
    @test kinds["ADVANCED_FACE"]==3
    @test kinds["PLANE"]==2 && kinds["B_SPLINE_SURFACE_WITH_KNOTS"]==1
    @test kinds["B_SPLINE_CURVE_WITH_KNOTS"]==2 && kinds["LINE"]==1
    @test kinds["EDGE_CURVE"]==3

    _,kinds=_entity_kinds("brep_occ_nurbs_pillow.step")
    @test kinds["ADVANCED_FACE"]==2
    @test kinds["B_SPLINE_SURFACE_WITH_KNOTS"]==2
    @test get(kinds,"PLANE",0)==0
    @test kinds["B_SPLINE_CURVE_WITH_KNOTS"]==4 && kinds["EDGE_CURVE"]==4
end

@testset "bump plate: sampled B-spline trims lift onto a quadratic patch" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_nurbs_bump_plate.step");
                        fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test size(surface.tets,2)==0
    @test nnodes(surface)==85 && ntris(surface)==166
    @test vec(minimum(surface.coords,dims=2))≈[0.,0.,0.] atol=1e-12
    @test vec(maximum(surface.coords,dims=2))≈[2.,1.,1.125] atol=1e-12
    @test all(values(_edge_incidence(surface)).==(2))
    @test _signed_surface_volume(surface)≈2.107791 atol=1e-6
    # The patch is a tensor Bézier: z = 1 + x(1-x/2)·y(1-y). Every lifted
    # node above the rim plane lies on the patch interior — including the
    # exact apex sample at (1,0.5,1.125).
    apexes=[i for i in axes(surface.coords,2) if surface.coords[3,i]≈1.125]
    @test length(apexes)>=1
    for i in axes(surface.coords,2)
        p=surface.coords[:,i]
        p[3]>1+1e-9 || continue
        residual=p[3]-1-(p[1]-p[1]^2/2)*p[2]*(1-p[2])
        @test abs(residual)<1e-9
    end

    filled=import_step(joinpath(FIXTURES,"brep_occ_nurbs_bump_plate.step"))
    @test filled isa E.Mesh && validate(filled).ok && ntets(filled)>0
    volumes=_tet_volumes(filled)
    @test all(>(0),volumes)
    @test sum(volumes)≈2.107791 atol=1e-6
end

@testset "periodic B-spline pipe: unwrapped uv strip welds the seam" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_nurbs_pipe.step");
                        fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test nnodes(surface)==112 && ntris(surface)==220
    @test vec(minimum(surface.coords,dims=2))[3]≈0.0 atol=1e-12
    @test vec(maximum(surface.coords,dims=2))[3]≈1.0 atol=1e-12
    @test all(values(_edge_incidence(surface)).==(2))
    @test _signed_surface_volume(surface)≈2.668243 atol=1e-6
    # The unwrapped seam welds positionally: no duplicated node anywhere,
    # and the seam line at theta=0 carries the two rim junctions.
    @test length(Set(Tuple(surface.coords[:,i]) for i in
                     axes(surface.coords,2)))==nnodes(surface)
    seam=[i for i in axes(surface.coords,2) if
          surface.coords[1,i]≈1.0 && abs(surface.coords[2,i])<1e-9]
    @test count(i->surface.coords[3,i]≈0.0,seam)>=1
    @test count(i->surface.coords[3,i]≈1.0,seam)>=1
    # Every wall node projects onto the rim B-spline with zero residual —
    # the unwrapped uv polygon was lifted through the NURBS evaluator, not
    # a polygonal approximation.
    curve=_rim_curve("brep_occ_nurbs_pipe.step")
    nwall=0;ncap=0
    for i in axes(surface.coords,2)
        p=surface.coords[:,i]
        r=_curve_residual(curve,p)
        if r<1e-9
            nwall+=1
            @test -1e-12<=p[3]<=1+1e-12
            @test hypot(p[1],p[2])<=1+1e-9
        elseif p[3]≈0.0 && r>0.05
            ncap+=1   # cap-interior lattice node strictly inside the rim
        end
    end
    @test nwall>20
    @test ncap>=1

    filled=import_step(joinpath(FIXTURES,"brep_occ_nurbs_pipe.step"))
    @test filled isa E.Mesh && validate(filled).ok && ntets(filled)>0
    volumes=_tet_volumes(filled)
    @test all(>(0),volumes)
    @test sum(volumes)≈2.668243 atol=1e-6
end

@testset "pillow: two B-spline patches share one sampled rim" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_nurbs_pillow.step");
                        fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test nnodes(surface)==130 && ntris(surface)==256
    @test vec(minimum(surface.coords,dims=2))≈[0.,0.,-0.125] atol=1e-12
    @test vec(maximum(surface.coords,dims=2))≈[2.,1.,0.125] atol=1e-12
    @test all(values(_edge_incidence(surface)).==(2))
    @test _signed_surface_volume(surface)≈0.215581 atol=1e-6
    @test any(i->surface.coords[3,i]≈0.125,axes(surface.coords,2))
    @test any(i->surface.coords[3,i]≈-0.125,axes(surface.coords,2))
    # Every node lies on one of the two Bézier patches |z|=x(1-x/2)·y(1-y).
    for i in axes(surface.coords,2)
        p=surface.coords[:,i]
        residual=abs(p[3])-(p[1]-p[1]^2/2)*p[2]*(1-p[2])
        @test abs(residual)<1e-9
    end

    filled=import_step(joinpath(FIXTURES,"brep_occ_nurbs_pillow.step"))
    @test filled isa E.Mesh && validate(filled).ok && ntets(filled)>0
    volumes=_tet_volumes(filled)
    @test all(>(0),volumes)
    @test sum(volumes)≈0.215581 atol=1e-6
end

@testset "trimmed-NURBS certificate rejects off-surface or broken trims" begin
    # Raising one boundary control of the top patch pulls the surface edge
    # off the unchanged rim trim — every sampled trim point must lift back
    # onto its face's patch within tolerance.
    _with_source("brep_occ_nurbs_bump_plate.step",
        "#126 = CARTESIAN_POINT('',(1.,0.,1.));"=>
        "#126 = CARTESIAN_POINT('',(1.,0.,1.4));") do path
        err=try import_step(path;fill=false);nothing catch e e end
        @test err isa ArgumentError
        @test occursin("B-spline surface",sprint(showerror,err))
    end
    # Denting one rim-curve control leaves the z=0 cap plane but pulls the
    # wall's trim off the tube patch — same on-surface certificate.
    _with_source("brep_occ_nurbs_pipe.step",
        "#12 = CARTESIAN_POINT('',(0.7071067811865476,0.7071067811865475,0.));"=>
        "#12 = CARTESIAN_POINT('',(0.5,0.5,0.));") do path
        err=try import_step(path;fill=false);nothing catch e e end
        @test err isa ArgumentError
        @test occursin("B-spline surface",sprint(showerror,err))
    end
    # Reversing one oriented edge use breaks the uv wire's traversal order.
    _with_source("brep_occ_nurbs_pillow.step",
        "#39 = ORIENTED_EDGE('',*,*,#25,.T.);"=>
        "#39 = ORIENTED_EDGE('',*,*,#25,.F.);") do path
        err=try import_step(path;fill=false);nothing catch e e end
        @test err isa ArgumentError
        @test occursin("ordered closed chain",sprint(showerror,err))
    end
end

end # module
