using Test, Tessella, SHA

const _RQU = Tessella.API
_rqu_bits(values) = [string(reinterpret(UInt64, Float64(x)); base=16, pad=16) for x in values]
function _rqu_box_count(value)
    value === Core.Box && return 1
    value isa GlobalRef && value.mod === Core && value.name === :Box && return 1
    if value isa Expr
        return sum(_rqu_box_count(item) for item in value.args; init=0)
    elseif value isa Core.CodeInfo
        return _rqu_box_count(value.code) + _rqu_box_count(value.slottypes) + _rqu_box_count(value.ssavaluetypes)
    elseif value isa Union{AbstractArray,Tuple,Core.SimpleVector}
        return sum(_rqu_box_count(item) for item in value; init=0)
    end
    return 0
end
function _rqu_real_box_factory()
    value = 1
    reader = () -> value
    value = 2
    return reader
end

@testset "Raw recombine actual method origin and lowered Box calibration" begin
    file = realpath(joinpath(Base.pkgdir(Tessella), "src", "interfaces", "API.jl"))
    probes = [(_RQU._entity_tri_elements, Tuple{GeoModel,Nothing,Nothing,Int,Int}),
        (_RQU._recombine, Tuple{})]
    records = Any[]
    for (callable, signature) in probes
        method = which(callable, signature)
        @test realpath(String(method.file)) == file
        lowered = Base.uncompressed_ast(method)
        @test _rqu_box_count(lowered) == 0
        @test all(pair -> _rqu_box_count(first(pair)) == 0, Base.code_typed(callable, signature; optimize=false))
        push!(records, (method=string(method), file=file, source_sha256=bytes2hex(sha256(read(file))),
            line=method.line, lowered=lowered))
    end
    positive_method = which(_rqu_real_box_factory, Tuple{})
    positive = Base.uncompressed_ast(positive_method)
    @test _rqu_box_count(positive) > 0
    reader = _rqu_real_box_factory()
    @test reader() == 2
    @test any(type -> type === Core.Box, fieldtypes(typeof(reader)))
    @info "Raw recombine actual methods and real reassigned-capture positive" records positive
end

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
