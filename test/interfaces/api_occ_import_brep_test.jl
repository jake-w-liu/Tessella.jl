module OCCImportBREPChecks
using Test, Tessella
using LinearAlgebra: det
const A=Tessella.API
const M=Tessella.Model
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function session(f)
    A.initialize()
    try
        A.option("General.Terminal",0)
        return f(A.model.occ)
    finally
        A.finalize()
    end
end

function trace_certificate(model)
    kernel=model.occ_kernel
    residual=0.0
    for ((dim,tag),fid) in kernel.bindings
        dim==2 || continue
        face=kernel.shapes[fid]
        g=face.geometry.geometry
        for wid in face.children,(eid,reversed) in kernel.shapes[wid].geometry.uses
            edge=kernel.shapes[eid]
            eg=edge.geometry.geometry
            etag=kernel.shape_tags[eid][2]
            entry=g.pcurves[etag]
            pc=reversed && entry.rev!==nothing ? entry.rev : entry.fwd
            for fraction in (0.0,.25,.5,.75,1.0)
                t=eg.t0+fraction*(eg.t1-eg.t0)
                uv=M._occ_pcurve_d0(pc,t)
                p=M.model_value(model,1,etag,[t])
                q=M.model_value(model,2,tag,collect(uv))
                residual=max(residual,sqrt(sum(abs2,p.-q)))
            end
        end
    end
    return residual
end

@testset "BREP imports actual shared analytic CAD topology" begin
    # Static fixtures were authored with the public OCCT-backed CAD writer:
    # unit box, R=.5 H=2 cylinder, R=1/.4 H=2 frustum, unit sphere,
    # and R=1 r=.25 torus. Runtime import and meshing are entirely native.
    session() do occ
        out=occ.import_shapes(joinpath(FIXTURES,"brep_native_analytic.brep"))
        @test first.(out)==fill(3,5)
        @test Tuple(length(occ.get_entities(d)) for d in 0:3)==(15,23,14,5)
        @test [occ.get_mass(d,t) for (d,t) in out]≈
            [1.,π/2,2π*(1+.4+.16)/3,4π/3,π^2/8] rtol=1e-8
        occ.synchronize()
        @test trace_certificate(A.CURRENT[])<2e-7
        occ.remove(out[2:end],true)
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(8,12,6,1)
        A.option("Mesh.Algorithm",6)
        A.option("Mesh.MeshSizeMin",2.)
        A.option("Mesh.MeshSizeMax",2.)
        A.mesh.generate(3)
        tags,xyz,_=A.mesh.get_nodes()
        coordinates=reshape(xyz,3,:)
        rows=Dict(tag=>row for (row,tag) in enumerate(tags))
        families,_,cells=A.mesh.get_elements(3)
        @test families==Int32[4]
        volumes=Float64[]
        incidence=Dict{NTuple{3,UInt64},Int}()
        for group in cells,i in 1:4:length(group)
            tet=group[i:i+3]
            p=coordinates[:,rows[tet[1]]]
            push!(volumes,det(hcat((coordinates[:,rows[tet[j]]].-p for j in 2:4)...))/6)
            for face in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
                key=Tuple(sort(UInt64[tet[k] for k in face]))
                incidence[key]=get(incidence,key,0)+1
            end
        end
        @test !isempty(volumes) && all(>(0),volumes)
        @test sum(volumes)≈1.0 rtol=1e-10
        sfamilies,_,faces=A.mesh.get_elements(2)
        @test sfamilies==Int32[2]
        surface=Set(Tuple(sort(UInt64[group[i],group[i+1],group[i+2]]))
            for group in faces for i in 1:3:length(group))
        @test surface==Set(face for (face,count) in incidence if count==1)
    end
    session() do occ
        mktempdir() do directory
            path=joinpath(directory,"renamed.dat")
            cp(joinpath(FIXTURES,"half_cyl.brep"),path)
            out=occ.import_shapes(path;format="brep")
            @test length(out)==1 && first(out[1])==3
            @test occ.get_mass(out[1]...)≈2π rtol=1e-9
            @test Tuple(length(occ.get_entities(d)) for d in 0:3)==(4,6,4,1)
            occ.synchronize()
            @test trace_certificate(A.CURRENT[])<2e-7
        end
    end
end

@testset "BREP implicit planar traces and root orientations" begin
    source=read(joinpath(FIXTURES,"brep_native_plane.brep"),String)
    for reversed in (false,true)
        session() do occ
            mktempdir() do directory
                path=joinpath(directory,"plane.brep")
                write(path,reversed ? replace(source,r"\+1 0\s*$"=>"-1 0\n") : source)
                out=occ.import_shapes(path)
                @test length(out)==1 && first(out[1])==2
                occ.synchronize()
                @test A.model.get_normal(out[1][2],[.5,.5])≈
                    (reversed ? [0.,0.,-1.] : [0.,0.,1.])
                @test occ.get_mass(out[1]...)≈1.
                @test trace_certificate(A.CURRENT[])<2e-7
            end
        end
    end
    session() do occ
        out=occ.import_shapes(joinpath(FIXTURES,"brep_native_nonlinear_pcurve.brep"))
        occ.synchronize()
        @test occ.get_mass(out[1]...)≈1.
        @test trace_certificate(A.CURRENT[])<2e-7
    end
end

function point_brep(locations,children;point="1 2 3")
    return "CASCADE Topology V1, (c) Matra-Datavision\n"*locations*
        "Curve2ds 0\nCurves 0\nPolygon3D 0\nPolygonOnTriangulations 0\n"*
        "Surfaces 0\nTriangulations 0\nTShapes 2\n"*
        "Ve\n1e-7\n$point\n0 0\n0101101\n*\n"*
        "Co\n1100000\n$children *\n+1 0\n"
end

@testset "BREP location identity preserves coincident CAD instances" begin
    identity="1\n1 0 0 0\n0 1 0 0\n0 0 1 0\n"
    for (locations,children,distinct) in (
            ("Locations 2\n"*identity*identity,"+2 1 +2 2",2),
            ("Locations 1\n"*identity,"+2 1 +2 1",1),
            ("Locations 1\n"*identity,"+2 0 +2 1",2))
        session() do occ
            mktempdir() do directory
                path=joinpath(directory,"points.brep")
                write(path,point_brep(locations,children))
                out=occ.import_shapes(path;highestDimOnly=false)
                @test length(out)==2
                @test length(Set(out))==distinct
                @test length(occ.get_entities(0))==distinct
                occ.synchronize()
                @test all(dt->A.model.get_value(dt...,Float64[])≈[1.,2.,3.],out)
            end
        end
    end
    session() do occ
        mktempdir() do directory
            path=joinpath(directory,"scaled point.brep")
            location="Locations 1\n1\n1e-200 0 0 0\n0 1e-200 0 0\n0 0 1e-200 0\n"
            write(path,point_brep(location,"+2 1";point="1e200 0 0"))
            out=occ.import_shapes(path)
            occ.synchronize()
            @test A.model.get_value(out[1]...,Float64[])≈[1.,0.,0.]
        end
    end
end

@testset "BREP compact DAG expansion is bounded before publication" begin
    session() do occ
        occ.add_point(2.,3.,4.)
        before=deepcopy(A.CURRENT[].occ_kernel)
        prefix=split(point_brep("Locations 0\n","+2 0"),"TShapes 2\n")[1]
        count=201
        source=prefix*"TShapes $count\nVe\n1e-7\n0 0 0\n0 0\n0101101\n*\n"
        for index in count-1:-1:1
            source*="Co\n1100000\n+$(index+1) 0 +$(index+1) 0 *\n"
        end
        source*="+1 0\n"
        mktempdir() do directory
            path=joinpath(directory,"compound DAG.brep")
            write(path,source)
            @test_throws ArgumentError occ.import_shapes(path;highestDimOnly=false)
            @test A.CURRENT[].occ_kernel==before
        end
    end
end

@testset "BREP rational and periodic geometry retains its parameters" begin
    session() do occ
        out=occ.import_shapes(joinpath(FIXTURES,"brep_native_spline.brep"))
        @test length(out)==1 && first(out[1])==2
        occ.synchronize()
        face=out[1][2]
        @test A.model.get_value(2,face,[.5,.5])≈[1.,1.,1.1025641025641026]
        midpoints=Vector{Float64}[]
        for (_,tag) in A.model.get_entities(1)
            lo,hi=A.model.get_parametrization_bounds(1,tag)
            push!(midpoints,A.model.get_value(1,tag,(lo.+hi)./2))
        end
        @test any(p->isapprox(p,[1.,.8235294117647058,0.]),midpoints)
        @test any(p->isapprox(p,[5.,2.,0.]),midpoints)
        @test trace_certificate(A.CURRENT[])<2e-7
        occ.translate(out,2.,3.,4.)
        occ.synchronize()
        @test A.model.get_value(2,face,[.5,.5])≈[3.,4.,5.1025641025641026]
        @test trace_certificate(A.CURRENT[])<2e-7
    end
end

@testset "BREP nonlinear periodic pcurves keep their continuous UV branch" begin
    session() do occ
        out=occ.import_shapes(joinpath(FIXTURES,"brep_native_nonlinear_cylinder.brep"))
        occ.synchronize()
        model=A.CURRENT[]
        face=model.surface_geometry[out[1][2]]
        tested=0
        for (_,tag) in A.model.get_entities(1)
            curve=model.curve_geometry[tag]
            curve.occ===:circle || continue
            pc=face.pcurves[tag].fwd
            hasproperty(pc,:brep_inverse2d) || continue
            values=[M._occ_pcurve_d0(pc,curve.t0+(curve.t1-curve.t0)*i/20)[1] for i in 0:20]
            @test all(u->-1e-9<=u<=2π+1e-9,values)
            @test issorted(values)
            @test values≈collect(range(0.,2π;length=21)) atol=2e-10
            tested+=1
        end
        @test tested==1
        @test trace_certificate(model)<2e-7
    end
end

@testset "BREP generic extrusion retains exact geometry under reflection" begin
    session() do occ
        out=occ.import_shapes(joinpath(FIXTURES,"brep_native_extrusion.brep"))
        occ.synchronize()
        @test A.model.get_value(2,out[1][2],[.25,-.75])≈[.5,.375,.75]
        @test trace_certificate(A.CURRENT[])<2e-7
        occ.mirror(out,1.,0.,0.,0.)
        occ.synchronize()
        @test A.model.get_value(2,out[1][2],[.25,-.75])≈[-.5,.375,.75]
        @test trace_certificate(A.CURRENT[])<2e-7
    end
end

@testset "Malformed BREP leaves pending native geometry unchanged" begin
    session() do occ
        occ.add_point(2.,3.,4.)
        before=deepcopy(A.CURRENT[].occ_kernel)
        source=read(joinpath(FIXTURES,"half_cyl.brep"),String)
        mktempdir() do directory
            path=joinpath(directory,"broken.brep")
            for invalid in ("not a BREP file",source[1:div(length(source),2)],
                            replace(source,"TShapes 21"=>"TShapes 2147483647"),
                            replace(source,"Curve2ds 13"=>"Curve2ds 2147483647"))
                write(path,invalid)
                @test_throws ArgumentError occ.import_shapes(path)
                @test A.CURRENT[].occ_kernel==before
            end
        end
    end
end
end
