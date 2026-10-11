using Test, Tessella, SHA
include(joinpath(@__DIR__, "..", "artifacts", "native_raw_recombine_uint64_states.jl"))

const _RQU = Tessella.API
_rqu_bits(values) = [string(reinterpret(UInt64, Float64(x)); base=16, pad=16) for x in values]
function _rqu_nodes(dim=-1, tag=-1)
    ids, xyz, uv = _RQU.mesh.get_nodes(dim, tag, false, true)
    @test ids isa Vector{UInt64}
    return Dict("ids" => ids, "xyz" => xyz, "xyz_bits" => _rqu_bits(xyz),
        "uv" => uv, "uv_bits" => _rqu_bits(uv))
end
function _rqu_elements(dim=-1, tag=-1)
    types, ids, connectivity = _RQU.mesh.get_elements(dim, tag)
    @test types isa Vector{Int32}
    @test all(tags -> tags isa Vector{UInt64}, ids)
    @test all(nodes -> nodes isa Vector{UInt64}, connectivity)
    return Dict("types" => types, "ids" => ids, "connectivity" => connectivity)
end
function _rqu_snapshot()
    entities = _RQU.model.get_entities()
    groups = _RQU.model.get_physical_groups()
    return Dict("global_nodes" => _rqu_nodes(), "global_elements" => _rqu_elements(),
        "per_entity" => [Dict("dim" => dim, "tag" => tag,
            "nodes" => _rqu_nodes(dim, tag), "elements" => _rqu_elements(dim, tag),
            "physical" => _RQU.model.get_physical_groups_for_entity(dim, tag)) for (dim, tag) in entities],
        "physical_groups" => [Dict("dim" => dim, "tag" => tag,
            "entities" => _RQU.model.get_entities_for_physical_group(dim, tag)) for (dim, tag) in groups],
        "max_node" => _RQU.mesh.get_max_node_tag(), "max_element" => _RQU.mesh.get_max_element_tag())
end
function _rqu_seed(row)
    _RQU.initialize()
    _RQU.option("General.AbortOnError", 2)
    _RQU.option("Mesh.Renumber", 0)
    _RQU.option("Mesh.RecombinationAlgorithm", 0)
    _RQU.model.add_discrete_entity(2, 1)
    tags = UInt64.(row["case"]["nodes"])
    _RQU.mesh.add_nodes(2, 1, tags, [0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.],
        [0.,0.,1.,0.,1.,1.,0.,1.])
    _RQU.mesh.add_elements_by_type(1, 2, [30001,30012], tags[[1,2,3,1,3,4]])
    _RQU.model.add_physical_group(2, [1]; tag=10)
    _RQU.model.add_discrete_entity(0, 41001)
    _RQU.mesh.add_nodes(0, 41001, [20003], [12.,12.,0.])
    _RQU.mesh.add_elements_by_type(41001, 15, [row["case"]["sentinel_requested"]], [20003])
    return nothing
end
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

@testset "Raw recombine UInt64 labels and native temporary/final births" begin
    for row in RAW_RECOMBINE_NATIVE_ROWS
        try
            _rqu_seed(row)
            msg_before=(_RQU.NATIVE_LOOKUP_QUERY_ERRORS[],_RQU.NATIVE_MSG_LAST_ERROR[])
            @test _rqu_snapshot() == row["states"][1]["payload"]
            @test (_RQU.NATIVE_LOOKUP_QUERY_ERRORS[],_RQU.NATIVE_MSG_LAST_ERROR[])==msg_before
            _RQU.mesh.remove_elements(0, 41001)
            @test _rqu_snapshot() == row["states"][2]["payload"]
            @test (_RQU.NATIVE_LOOKUP_QUERY_ERRORS[],_RQU.NATIVE_MSG_LAST_ERROR[])==msg_before
            _RQU.mesh.set_recombine(2, 1, row["case"]["angle"])
            @test _RQU.mesh.recombine() === nothing
            @test _rqu_snapshot() == row["states"][3]["payload"]
            @test (_RQU.NATIVE_LOOKUP_QUERY_ERRORS[],_RQU.NATIVE_MSG_LAST_ERROR[])==msg_before

        finally
            _RQU.finalize()
        end
    end
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
