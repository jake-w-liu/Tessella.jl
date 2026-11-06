module APIMetricSurfacePeriodicChecks
using Test, Tessella, LinearAlgebra
const A=Tessella.API

function pair(;rotated=false,recombination=nothing,smoothing=0)
    A.initialize()
    try
        for (name,value) in (("General.Terminal",0),("Mesh.Algorithm",7),
                ("Mesh.Smoothing",0),("Mesh.Renumber",0),
                ("Mesh.MeshSizeFromPoints",0),("Mesh.MeshSizeFromCurvature",0),
                ("Mesh.MeshSizeExtendFromBoundary",0))
            A.option(name,value)
        end
        A.model.add("owned metric periodic pair")
        surfaces=Int[];boundary=Vector{Int}[]
        for offset in (0.,2.)
            points=[A.model.add_point(x+offset,y,0.;meshSize=.2) for
                (x,y) in ((0.,0.),(1.,0.),(1.,1.),(0.,1.))]
            curves=[A.model.add_line(points[i],points[mod1(i+1,4)]) for i in 1:4]
            for curve in curves
                A.mesh.set_transfinite_curve(curve,recombination in (2,3) ? 7 : 6)
            end
            surface=A.model.add_plane_surface([A.model.add_curve_loop(curves)])
            A.mesh.set_smoothing(2,surface,smoothing)
            if recombination!==nothing
                A.option("Mesh.RecombinationAlgorithm",recombination)
                A.mesh.set_recombine(2,surface)
            end
            push!(surfaces,surface);push!(boundary,curves)
        end
        field=A.mesh.field.add("MathEvalAniso")
        diagonal=rotated ? (16.,16.,1.) : (16.,4.,1.)
        for (name,value) in zip(("M11","M22","M33","M12","M13","M23"),
                (diagonal...,0.,0.,0.))
            A.mesh.field.set_string(field,name,string(value))
        end
        A.mesh.field.set_as_background_mesh(field)
        affine=rotated ? [0.,-1.,0.,3.,1.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.] :
            [1.,0.,0.,2.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.]
        A.mesh.set_periodic(2,[surfaces[2]],[surfaces[1]],affine)
        A.mesh.generate(2)
        ids,xyz,_=A.mesh.get_nodes()
        coordinates=Dict(uid=>xyz[3i-2:3i] for (i,uid) in enumerate(ids))
        @test allunique(ids)
        master=A.mesh.get_elements(2,surfaces[1])
        slave=A.mesh.get_elements(2,surfaces[2])
        @test master[1]==slave[1] && !isempty(master[1])
        @test length.(master[2])==length.(slave[2])
        @test allunique(vcat(master[2]...,slave[2]...))
        link=A.mesh.get_periodic_nodes(2,surfaces[2])
        @test link.master_entity==surfaces[1]
        @test allunique(link.master_nodes) && allunique(link.slave_nodes)
        mapping=Dict(zip(link.master_nodes,link.slave_nodes))
        @test Set(link.master_nodes)==Set(vcat(master[3]...))
        @test Set(link.slave_nodes)==Set(vcat(slave[3]...))
        @test all(zip(master[1],master[3],slave[3])) do (kind,first,last)
            width=kind==2 ? 3 : 4
            mapped=Set(Tuple(sort!([mapping[uid] for uid in cell])) for
                cell in eachcol(reshape(first,width,:)))
            mapped==Set(Tuple(sort!(collect(cell))) for cell in eachcol(reshape(last,width,:)))
        end
        transform=reshape(affine,4,4)'
        # Retained Curve middles are re-evaluated on their own CAD support.
        # Their inversion accuracy is 1e-10; the ordinary copy certificate
        # permits 1e-8*lc. Keep this tiny fixture stricter than that contract.
        @test all(mapping) do (first,last)
            isapprox(coordinates[last],(transform*[coordinates[first];1.])[1:3];atol=1e-10,rtol=0)
        end
        # Point/Curve copies keep their own actual lower products. Every
        # lower Curve cell remains a literal edge of its host Surface.
        for (surface,curves) in zip(surfaces,boundary)
            kinds,_,nodes=A.mesh.get_elements(2,surface)
            incidence=Dict{Tuple{UInt64,UInt64},Int}()
            area=0.
            for (kind,flat) in zip(kinds,nodes)
                width=kind==2 ? 3 : kind==3 ? 4 : error("unexpected Surface family")
                for cell in eachcol(reshape(flat,width,:))
                    for i in 1:width
                        edge=minmax(cell[i],cell[mod1(i+1,width)])
                        incidence[edge]=get(incidence,edge,0)+1
                    end
                    p=coordinates[cell[1]]
                    for i in 2:width-1
                        area+=norm(cross(coordinates[cell[i]]-p,coordinates[cell[i+1]]-p))/2
                    end
                end
            end
            @test isapprox(area,1.;atol=2e-12) && all(in((1,2)),values(incidence))
            for curve in curves
                kinds,_,nodes=A.mesh.get_elements(1,curve)
                @test kinds==Int32[1] && all(eachcol(reshape(only(nodes),2,:))) do edge
                    get(incidence,minmax(edge...),0)==1
                end
            end
            if recombination===nothing
                @test all(keys(incidence)) do (first,last)
                    delta=coordinates[last]-coordinates[first]
                    sqrt(sum(diagonal[i]*delta[i]^2 for i in 1:3))<=1+2e-12
                end
            end
        end
        return (ids,xyz,master,slave,link)
    finally
        A.finalize()
    end
end

@testset "Metric Surface masters publish final constructors to periodic copies" begin
    pair()
    pair(;rotated=true)
    pair(;recombination=2)
    pair(;rotated=true,recombination=3)
    # Per-Surface elliptic smoothing belongs to transfinite meshes. It must
    # neither reject an ordinary metric front nor alter its owned products.
    @test pair(;smoothing=3)==pair()
end
end
