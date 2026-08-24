"""
    CLI

Command-line façade: `tessella file.geo -2|-3` executes the bounded `.geo`
subset and writes `file.msh`. Unknown flags and OCC-only files are blockers.
"""
module CLI

using ..GeoExec: execute_geo
using ..IO: write_msh

export main

function main(args::Vector{String})
    isempty(args) && throw(ArgumentError("tessella: missing input file"))
    path=nothing; dim=0; outfile=nothing
    i=1
    while i<=length(args)
        a=args[i]
        if a in ("-2","-3")
            dim=parse(Int,a[2:2])
        elseif a=="-o"
            i+=1; i<=length(args) || throw(ArgumentError("tessella: -o needs a path"))
            outfile=args[i]
        elseif startswith(a,"-")
            throw(ArgumentError("tessella: unsupported flag $a"))
        else
            path===nothing || throw(ArgumentError("tessella: multiple input files"))
            path=a
        end
        i+=1
    end
    path===nothing && throw(ArgumentError("tessella: missing input file"))
    result=execute_geo(path; mesh_dim=dim)
    if dim>0
        result.mesh===nothing && throw(ErrorException("tessella: no mesh produced"))
        dest=outfile===nothing ? replace(path, r"\.(geo|GEO)$"=>".msh") : outfile
        write_msh(dest, result.mesh; version=4.1)
        return dest
    end
    return path
end

end # module
