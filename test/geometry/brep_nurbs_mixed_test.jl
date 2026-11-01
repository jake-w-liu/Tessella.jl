module NurbsMixedChecks
using Test,Tessella,LinearAlgebra
using Tessella.BRep: parse_step_entities
using Tessella.MeshTypes: nnodes,ntris,ntets,validate
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

# Every undirected surface edge must be bounded by exactly two triangles.
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

@testset "mixed-shell entity audit" begin
    _,kinds=_entity_kinds("brep_occ_nurbs_cyl_dome.step")
    @test kinds["MANIFOLD_SOLID_BREP"]==1 && kinds["CLOSED_SHELL"]==1
    @test kinds["ADVANCED_FACE"]==3
    @test kinds["CYLINDRICAL_SURFACE"]==1 && kinds["PLANE"]==1
    @test kinds["CIRCLE"]==2 && kinds["LINE"]==1 && kinds["EDGE_CURVE"]==3

    _,kinds=_entity_kinds("brep_occ_nurbs_sph_cap.step")
    @test kinds["ADVANCED_FACE"]==3
    @test kinds["SPHERICAL_SURFACE"]==1 && kinds["PLANE"]==1
    @test kinds["CIRCLE"]==3 && kinds["EDGE_CURVE"]==3

    _,kinds=_entity_kinds("brep_occ_nurbs_sph_seam.step")
    @test kinds["ADVANCED_FACE"]==3
    @test kinds["SPHERICAL_SURFACE"]==1 && kinds["PLANE"]==1
    # seam+pole layout: meridian seam arcs, degenerate pole LINE edge
    @test kinds["CIRCLE"]==3 && kinds["LINE"]==1 && kinds["EDGE_CURVE"]==5

    _,kinds=_entity_kinds("brep_occ_nurbs_nurbs_arb.step")
    @test kinds["ADVANCED_FACE"]==4
    @test kinds["B_SPLINE_SURFACE_WITH_KNOTS"]==2 && kinds["PLANE"]==2
    @test kinds["B_SPLINE_CURVE_WITH_KNOTS"]==6 && kinds["EDGE_CURVE"]==6

    _,kinds=_entity_kinds("brep_occ_nurbs_pcurve_plate.step")
    @test kinds["ADVANCED_FACE"]==6
    @test kinds["B_SPLINE_SURFACE_WITH_KNOTS"]==1 && kinds["PLANE"]==5
    @test kinds["SURFACE_CURVE"]==4 && kinds["PCURVE"]==4
    @test kinds["DEFINITIONAL_REPRESENTATION"]==4
    # 4 3-D rim B-splines + 3 2-D pcurve B-splines; 8 3-D LINEs + 1 2-D.
    @test kinds["B_SPLINE_CURVE_WITH_KNOTS"]==7 && kinds["LINE"]==9
end

@testset "cylinder dome: NURBS pole patch welded to a cylinder wall" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_nurbs_cyl_dome.step");
                        fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test nnodes(surface)==97 && ntris(surface)==190
    @test vec(minimum(surface.coords,dims=2))≈[-1.,-1.,0.] atol=1e-12
    @test vec(maximum(surface.coords,dims=2))≈[1.,1.,2.] atol=1e-12
    @test all(values(_edge_incidence(surface)).==(2))
    @test _signed_surface_volume(surface)≈4.176838 atol=1e-6
    # The cylinder wall is r=1 for z in [0,1]; the rational dome is the
    # unit hemisphere centred at (0,0,1) closing at the pole (0,0,2).
    nwall=0;ndome=0;nbase=0;npole=0
    for i in axes(surface.coords,2)
        p=surface.coords[:,i]
        if p[3]>1+1e-9
            ndome+=1
            @test abs(p[1]^2+p[2]^2+(p[3]-1)^2-1)<1e-9
        elseif p[3]>1e-9
            nwall+=1
            @test abs(hypot(p[1],p[2])-1)<1e-9
        else
            nbase+=1
            @test hypot(p[1],p[2])<=1+1e-9
        end
        p[3]≈2 && (npole+=1)
    end
    # The bare-rim dome closes through its collapsed pole row and emits a
    # certified pole fan over the 48 rim samples — the only node above the
    # equator plane is the welded pole itself.
    @test ndome==1 && nwall>10 && nbase>=1 && npole==1

    filled=import_step(joinpath(FIXTURES,"brep_occ_nurbs_cyl_dome.step"))
    @test filled isa Tessella.Mesh && validate(filled).ok && ntets(filled)>0
    volumes=_tet_volumes(filled)
    @test all(>(0),volumes)
    @test sum(volumes)≈4.176838 atol=1e-6
end

@testset "sphere cap: bare latitude wire closes through the pole" begin
    for fixture in ("brep_occ_nurbs_sph_cap.step",
                    "brep_occ_nurbs_sph_seam.step")
        surface=import_step(joinpath(FIXTURES,fixture);fill=false)
        @test surface isa Tessella.Mesh && validate(surface).ok
        expected=fixture=="brep_occ_nurbs_sph_cap.step" ? (440,876) : (760,1516)
        @test nnodes(surface)==expected[1] && ntris(surface)==expected[2]
        @test vec(minimum(surface.coords,dims=2))[3]≈0.0 atol=1e-12
        @test vec(maximum(surface.coords,dims=2))[3]≈1.5 atol=1e-12
        @test all(values(_edge_incidence(surface)).==(2))
        # Every curved-face node sits on the unit sphere about (0,0,0.5);
        # base nodes close z=0 inside the bottom rim radius.
        nsph=0;nbase=0;npole=0
        for i in axes(surface.coords,2)
            p=surface.coords[:,i]
            if p[3]>1e-9
                nsph+=1
                @test abs(p[1]^2+p[2]^2+(p[3]-0.5)^2-1)<1e-9
            else
                nbase+=1
                @test hypot(p[1],p[2])<=cos(pi/6)+1e-9
            end
            p[3]≈1.5 && (npole+=1)
        end
        @test nsph>50 && nbase>=1 && npole==1
    end
    # The two STEP layouts triangulate the same solid: their signed
    # volumes agree to the faceting tolerance.
    cap=import_step(joinpath(FIXTURES,"brep_occ_nurbs_sph_cap.step");fill=false)
    seam=import_step(joinpath(FIXTURES,"brep_occ_nurbs_sph_seam.step");fill=false)
    @test _signed_surface_volume(cap)≈3.254730 atol=1e-6
    @test _signed_surface_volume(seam)≈3.511700 atol=1e-6
    filled=import_step(joinpath(FIXTURES,"brep_occ_nurbs_sph_cap.step"))
    @test validate(filled).ok && ntets(filled)>0
end

@testset "NURBS+NURBS canoe: arbitrary patches share both seam curves" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_nurbs_nurbs_arb.step");
                        fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test nnodes(surface)==158 && ntris(surface)==312
    @test vec(minimum(surface.coords,dims=2))[1]≈0.0 atol=1e-12
    @test vec(maximum(surface.coords,dims=2))[1]≈3.0 atol=1e-12
    @test all(values(_edge_incidence(surface)).==(2))
    @test _signed_surface_volume(surface)≈1.894840 atol=1e-6
    # Shared seams are the authored quadratic B-splines:
    #   s1 (y=0):   x=2t+t², z=-0.8t(1-t)
    #   s2 (y=1.5): x=3t,    z=+0.8t(1-t)
    ns1=0;ns2=0;ncapx=0
    for i in axes(surface.coords,2)
        p=surface.coords[:,i]
        if abs(p[2])<1e-9
            t=-1+sqrt(1+p[1])
            ns1+=1
            @test abs(p[3]+0.8*t*(1-t))<1e-9
        elseif abs(p[2]-1.5)<1e-9
            t=p[1]/3
            ns2+=1
            @test abs(p[3]-0.8*t*(1-t))<1e-9
        end
        if p[1]≈0 || p[1]≈3
            # End-cap lens between the two authored rim quadratics.
            y=p[2]/1.5
            lo=-1.1*y*(1-y);hi=0.7*y*(1-y)
            ncapx+=1
            @test lo-1e-9<=p[3]<=hi+1e-9
        end
    end
    @test ns1>5 && ns2>5 && ncapx>10

    filled=import_step(joinpath(FIXTURES,"brep_occ_nurbs_nurbs_arb.step"))
    @test validate(filled).ok && ntets(filled)>0
end

@testset "authored PCURVE trims: uv wires come from the 2-D curves" begin
    # The pcurve plate is the bump-plate solid whose four top-face trims
    # carry SURFACE_CURVE->PCURVE representations: one 2-D LINE (south rim)
    # and three 2-D quadratic B-splines. The certified import must match
    # the projected twin bitwise in counts and volume.
    ents,kinds=_entity_kinds("brep_occ_nurbs_pcurve_plate.step")
    surface=import_step(joinpath(FIXTURES,"brep_occ_nurbs_pcurve_plate.step");
                        fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test nnodes(surface)==85 && ntris(surface)==166
    @test all(values(_edge_incidence(surface)).==(2))
    @test _signed_surface_volume(surface)≈2.107791 atol=1e-6
    # Same patch residual as the bump plate: z = 1 + x(1-x/2)·y(1-y).
    for i in axes(surface.coords,2)
        p=surface.coords[:,i]
        p[3]>1+1e-9 || continue
        residual=p[3]-1-(p[1]-p[1]^2/2)*p[2]*(1-p[2])
        @test abs(residual)<1e-9
    end
    # The four top-face uses each resolve an authored pcurve on the
    # B-spline face's own surface id; plane uses must not (their basis is
    # the patch, never the plane).
    fe=only(e for e in values(ents) if
            e.kind=="ADVANCED_FACE" && BRep._deref_entity(ents,e.args[3]).kind==
            "B_SPLINE_SURFACE_WITH_KNOTS")
    curvekinds=Set{String}()
    bounds=BRep._step_face_bound_uses(ents,fe,curvekinds,"test")
    surfid=fe.args[3][2]
    pckinds=Symbol[]
    for uses in bounds,u in uses
        pc=BRep._step_use_pcurve(ents,u,surfid,"test")
        pc!==nothing && push!(pckinds,pc.kind)
    end
    @test sort(pckinds)==[:bspline,:bspline,:bspline,:line]

    filled=import_step(joinpath(FIXTURES,"brep_occ_nurbs_pcurve_plate.step"))
    @test validate(filled).ok && ntets(filled)>0
end

@testset "mixed-shell certificate rejects off-surface pcurves and rims" begin
    # Perturbing the south rim's authored LINE pcurve so it leaves the
    # patch boundary must fail the on-surface lift certificate.
    _with_source("brep_occ_nurbs_pcurve_plate.step",
        "#74 = DIRECTION('',(1.,0.));"=>
        "#74 = DIRECTION('',(1.,0.25));") do path
        err=try import_step(path;fill=false);nothing catch e e end
        @test err isa ArgumentError
        @test occursin("B-spline surface",sprint(showerror,err))
    end
    # Inflating the dome's shared rim circle off the cylinder wall must
    # fail the rim on-surface certificate.
    _with_source("brep_occ_nurbs_cyl_dome.step",
        "#15 = CIRCLE('',#12,1.);"=>
        "#15 = CIRCLE('',#12,1.2);") do path
        err=try import_step(path;fill=false);nothing catch e e end
        @test err isa ArgumentError
    end
end

end # module
