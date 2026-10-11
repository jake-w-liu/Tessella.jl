using TOML

@testset "Surface recombination uses published physical coordinates" begin
    model=Tessella.Model
    for normal in ((-1.,-1.,1.),(-.5,0.,1.),(0.,0.,1.)),u in (2.,6.,10.),v in (3.,11.,15.)
        plane=(normal=normal,axes=(1,2),k=3,anchor=(4.,8.,16.))
        frame=model._surface_metric_frame(plane)
        x,y=model._surface_metric_point(frame,u,v)
        mesh=Mesh(reshape([x,y,0.],3,1))
        canonical=Dict{NTuple{2,Float64},NTuple{3,Float64}}()
        physical=model._surface_recombine_physical_coords(mesh,frame,canonical)
        expected=(u,v,16.0-normal[1]*(u-4.0)-normal[2]*(v-8.0))
        @test all(isapprox.(vec(physical),collect(expected);rtol=2eps(),atol=0.))
        published=model._surface_metric_lift(frame,x,y)
        @test reinterpret(UInt64,vec(physical))==reinterpret(UInt64,collect(published))
        canonical[(x==0 ? 0.0 : x,y==0 ? 0.0 : y)]=expected
        restored=model._surface_recombine_physical_coords(mesh,frame,canonical)
        @test reinterpret(UInt64,vec(restored))==reinterpret(UInt64,collect(expected))
        @test restored!==mesh.coords
    end
end

@testset "Pinned native tilted-plane candidate priority" begin
    fixture=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","model_recombine_physical_priority.toml"))
    @test fixture["primary_version"]=="4.15.2"
    difference_count=0
    function selected_cells(result)
        cells=Vector{Int}[]
        for block in result.blocks,column in axes(block.nodes,2)
            push!(cells,vcat(Int(block.msh),sort(Int.(block.nodes[:,column]))))
        end
        sort!(cells;by=Tuple)
    end
    for control in fixture["controls"]
        projected=reduce(hcat,control["metric"])
        physical=reduce(hcat,control["world"])
        triangles=Int32.(reduce(hcat,control["triangles"]))
        input=Mesh(projected;tris=triangles)
        result=recombine_triangles(input;algorithm=:greedy,recombine_angle=45.,
            recombine_angle_coords=physical,recombine_priority_coords=physical)
        @test selected_cells(result)==control["physical_cells"]
        @test result.coords==projected
        @test result.coords!==input.coords
        @test validate(result).ok
        legacy=recombine_triangles(input;algorithm=:greedy,recombine_angle=45.)
        @test selected_cells(legacy)==control["projected_cells"]
        difference_count+=control["physical_cells"]!=control["projected_cells"]
    end
    @test difference_count==4
    input=Mesh(Float64[0 1 1 0;0 0 1 1;0 0 0 0];tris=Int32[1 1;2 3;3 4])
    for invalid in (zeros(2,4),fill(NaN,3,4),fill(Inf,3,4),fill(true,3,4),fill(1+im,3,4))
        @test_throws ArgumentError recombine_triangles(input;recombine_angle=45.,recombine_priority_coords=invalid)
    end
end
