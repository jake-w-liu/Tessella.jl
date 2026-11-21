# Gmsh's own .geo files, Tessella side: every file of `tutorials` and
# `examples/simple_geo` in a Gmsh source tree is opened and meshed in one,
# two and three dimensions through the session API. `gmsh_examples_gmsh.py`
# does the same through Gmsh's Python API; see README.md.
#
#     julia --project=. --startup-file=no -O2 --threads=1 benchmark/gmsh_examples.jl <gmsh source> [results.toml]
#
# Every file is run twice in the same process and the second run is reported,
# so compilation is excluded; garbage left by earlier files is collected before
# the clock starts. A file that stops with an error keeps the steps
# it completed and the first line of the error.

using Tessella
using Printf

const A = Tessella.API

const DIRECTORIES = ("tutorials", "examples/simple_geo")

# Files that need what neither driver provides: a post-processing view to
# merge (t8, t9), a surface mesh to remesh (t13) or a background mesh (t17).
const SKIPPED = ("t8.geo", "t9.geo", "t13.geo", "t17.geo")

count_cells(dim) = (blocks = A.mesh.get_elements(dim)[2]; isempty(blocks) ? 0 : sum(length, blocks))

function measure(path)
    opened = NaN
    seconds = Float64[]; nodes = Int[]; surface = Int[]; volume = Int[]
    failure = ""
    A.initialize()
    try
        A.option("General.Terminal", 0)
        GC.gc()     # the previous file's garbage is not this file's time
        opened = @elapsed A.open(path)
        for dim in 1:3
            push!(seconds, @elapsed A.mesh.generate(dim))
            push!(nodes, length(A.mesh.get_nodes()[1]))
            push!(surface, count_cells(2))
            push!(volume, count_cells(3))
        end
    catch err
        step = isnan(opened) ? "open" : "generate($(length(seconds) + 1))"
        failure = step * ": " * first(replace(sprint(showerror, err), r"\s+" => " "), 120)
    finally
        A.finalize()
    end
    return (; opened, seconds, nodes, surface, volume, failure)
end

function main(arguments)
    isempty(arguments) && error("usage: gmsh_examples.jl <gmsh source directory> [results.toml]")
    root = arguments[1]
    io = length(arguments) > 1 ? open(arguments[2], "w") : stdout
    try
        println(io, "# Tessella $(pkgversion(Tessella)), Julia $(VERSION), $(Threads.nthreads()) thread(s)")
        for directory in DIRECTORIES, name in sort!(readdir(joinpath(root, directory)))
            endswith(name, ".geo") || continue
            name in SKIPPED && continue
            path = joinpath(root, directory, name)
            measure(path)
            result = measure(path)
            println(io, "\n[[file]]")
            println(io, "name = \"$directory/$name\"")
            @printf(io, "open_seconds = %.4f\n", isnan(result.opened) ? -1.0 : result.opened)
            println(io, "generate_seconds = [", join((@sprintf("%.4f", s) for s in result.seconds), ", "), "]")
            println(io, "nodes = [", join(result.nodes, ", "), "]")
            println(io, "surface_cells = [", join(result.surface, ", "), "]")
            println(io, "volume_cells = [", join(result.volume, ", "), "]")
            # A TOML string: no quote or backslash from the message survives.
            println(io, "error = \"", replace(result.failure, '"' => '\'', '\\' => '/'), "\"")
            flush(io)
        end
    finally
        io === stdout || close(io)
    end
end

main(ARGS)
