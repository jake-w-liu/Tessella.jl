using Documenter
using Tessella

# A page per module keeps the full docstring index browsable and searchable.
function documented_modules(root)
    result = Module[root]
    for name in names(root; all=true)
        isdefined(root, name) || continue
        child = getfield(root, name)
        child isa Module && child !== root && parentmodule(child) === root || continue
        append!(result, documented_modules(child))
    end
    return unique(result)
end

modules = documented_modules(Tessella)
DocMeta.setdocmeta!(Tessella, :DocTestSetup, :(using Tessella); recursive=true)
api_dir = joinpath(@__DIR__, "src", "api")
isdir(api_dir) && rm(api_dir; recursive=true)
mkpath(api_dir)
api_pages = Pair{String,String}[]
for mod in modules
    name = string(mod)
    filename = name * ".md"
    open(joinpath(api_dir, filename), "w") do io
        println(io, "# ", name, "\n\n```@meta\nCurrentModule = ", name, "\n```")
        println(io, "\n```@autodocs\nModules = [Main.", name, "]\nPublic = true\nPrivate = true\n```")
    end
    push!(api_pages, name => "api/" * filename)
end
references = sort(filter(f -> endswith(f, ".md") && f != "index.md",
    readdir(joinpath(@__DIR__, "src", "reference"))))
reference_pages = ["Overview" => "reference/index.md"]
for file in references
    title = replace(first(readlines(joinpath(@__DIR__, "src", "reference", file))), r"^#+\s*" => "")
    push!(reference_pages, title => "reference/" * file)
end

makedocs(;
    modules=[Tessella],
    sitename="Tessella.jl",
    authors="Jake Liu",
    format=Documenter.HTML(;
        canonical="https://jake-w-liu.github.io/Tessella.jl/stable/",
        prettyurls=get(ENV, "CI", "false") == "true"),
    pages=["Getting started" => "index.md", "Feature reference" => reference_pages,
        "API" => api_pages, "Contributing" => "contributing.md"],
    checkdocs=:all,
    warnonly=false,
)

if get(ENV, "GITHUB_ACTIONS", "false") == "true"
    deploydocs(; repo="github.com/jake-w-liu/Tessella.jl.git", devbranch="main")
end
