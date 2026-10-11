# Native model-data partitions belong to each populated step context.
# A kind rebuild resets all contexts; negative partition identifiers do not
# enter the set used by that step's MSH integer header.
function _native_view_partition_plan(previous::Set{Int},partition::Int;reset::Bool)
    result=reset ? Set{Int}() : copy(previous)
    partition>=0 && push!(result,partition)
    return result
end

function _native_view_partition_header(step::Int,components::Int,rows::Int,partitions::Set{Int})
    step>=0 && components>0 && rows>=0 || throw(ArgumentError(
        "model view: invalid data-section header"))
    all(partition->partition>=0,partitions) || throw(ArgumentError(
        "model view: invalid partition history"))
    header=Int64[step,components,rows]
    length(partitions)==1 && push!(header,only(partitions))
    return header
end
