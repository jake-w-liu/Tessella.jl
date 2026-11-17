using Test, Tessella

const _RQU = Tessella.API

@testset "Raw UInt64 triangle collection keeps linear allocation growth" begin
    previous = 0
    for count in (64,128,256)
        model = GeoModel()
        Tessella.Model.add_discrete_entity!(model, 2, 1)
        labels = UInt64[(UInt64(1)<<63)+UInt64(7index) for index in 1:3count]
        xyz = repeat([0.,0.,0.], 3count)
        Tessella.Model.add_discrete_nodes!(model, 2, 1, labels, xyz)
        elements = UInt64.(1:count)
        Tessella.Model.add_discrete_elements!(model, 2, 1, [2], [elements], [labels])
        collect_rows() = _RQU._entity_tri_elements(model, nothing, nothing, 2, 1)
        collect_rows(); collect_rows()
        samples = [@allocated(collect_rows()) for _ in 1:2]
        allocated = minimum(samples)
        @test allocated <= 32768 + 1024count
        previous == 0 || @test allocated <= 2.5previous + 262144
        entries, dense_columns, _ = collect_rows()
        @test length(entries) == count && isempty(dense_columns)
        @test all(entry -> entry[1] isa UInt64 && entry[2] isa Vector{UInt64}, entries)
        @info "Raw UInt64 lineage collection allocations" count samples
        previous = allocated
    end
end
