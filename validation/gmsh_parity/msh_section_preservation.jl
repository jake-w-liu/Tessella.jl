#!/usr/bin/env julia
# P2: MSH ancillary/unknown-section and view-data preservation vs Gmsh 4.15.2.
#
# Gmsh writes a meshed square carrying one view of each model-data kind
# ($NodeData, $ElementData, $ElementNodeData) plus an $InterpolationScheme
# ancillary section, in both ASCII and binary form. Tessella must read every
# section, re-emit them into every supported output mode without losing the
# data or its tag bindings, and stay readable by Gmsh itself.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
const Elements = Tessella.Elements

function find_gmsh_api()
    explicit=get(ENV,"GMSH_JULIA_API","")
    !isempty(explicit) && isfile(explicit) && return explicit
    executable=Sys.which("gmsh")
    executable===nothing && error("gmsh is not on PATH")
    prefix=dirname(dirname(realpath(executable)))
    for path in (joinpath(prefix,"lib","gmsh.jl"),
                 "/opt/homebrew/opt/gmsh/lib/gmsh.jl")
        isfile(path) && return path
    end
    error("could not locate gmsh.jl")
end
include(find_gmsh_api())

const COMMENT = "\$Comments\ntessella preserved comment\n  indented\n\$EndComments\n"

# Sample points kept away from element edges so probes land unambiguously.
const PROBES = [(0.23, 0.31, 0.0), (0.61, 0.72, 0.0), (0.87, 0.11, 0.0)]

function build_fixture(dir)
    gmsh.model.add("views")
    p1=gmsh.model.geo.addPoint(0,0,0,0.35)
    p2=gmsh.model.geo.addPoint(1,0,0,0.35)
    p3=gmsh.model.geo.addPoint(1,1,0,0.35)
    p4=gmsh.model.geo.addPoint(0,1,0,0.35)
    cl=gmsh.model.geo.addCurveLoop([
        gmsh.model.geo.addLine(p1,p2), gmsh.model.geo.addLine(p2,p3),
        gmsh.model.geo.addLine(p3,p4), gmsh.model.geo.addLine(p4,p1)])
    gmsh.model.geo.addPlaneSurface([cl])
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    nx = Dict(tag => ncoords[3i-2] for (i,tag) in pairs(ntags))
    etypes, etags, enodes = gmsh.model.mesh.getElements(2)
    etags_flat = collect(etags[findfirst(==(2), etypes)])
    enodes_flat = collect(enodes[findfirst(==(2), etypes)])

    # NodeData: each node carries its x coordinate, so a probe at p returns p.x
    v_n = gmsh.view.add("nodal")
    gmsh.view.addModelData(v_n,0,"views","NodeData",
        collect(ntags),[[nx[t]] for t in ntags])
    # ElementData: each element carries its centroid x
    v_e = gmsh.view.add("elemental")
    gmsh.view.addModelData(v_e,0,"views","ElementData",etags_flat,
        [[sum(nx[enodes_flat[3i-2+k]] for k in 0:2)/3]
         for i in eachindex(etags_flat)])
    # ElementNodeData: each element-node carries its x coordinate
    v_en = gmsh.view.add("elnode")
    gmsh.view.addModelData(v_en,0,"views","ElementNodeData",etags_flat,
        [[nx[enodes_flat[3i-2+k]] for k in 0:2] for i in eachindex(etags_flat)])

    paths = Dict{Symbol,String}()
    for (binary,key) in ((false,:ascii),(true,:binary))
        gmsh.option.setNumber("Mesh.Binary",binary ? 1 : 0)
        singles = String[]
        for (v,name) in ((v_n,"nodal"),(v_e,"elemental"),(v_en,"elnode"))
            path = joinpath(dir,"single-$key-$name.msh")
            gmsh.view.write(v,path)
            push!(singles,path)
        end
        # gmsh.view.write emits mesh + one view per file; splice the data and
        # interpolation-scheme sections of the last two singles into the
        # first so one file carries all three views.
        base = read(singles[1])
        for extra in singles[2:3]
            bytes = read(extra)
            for name in ("InterpolationScheme","NodeData","ElementData",
                         "ElementNodeData")
                marker = codeunits("\$"*name)
                firstidx = findfirst(marker,bytes)
                firstidx===nothing && continue
                endmark = codeunits("\$End"*name)
                lastidx = findfirst(endmark,bytes[firstidx[end]+1:end])
                lastidx===nothing &&
                    error("unterminated \$"*name*" in $extra")
                endpos = firstidx[end]+lastidx[end]
                tail = findfirst(==(0x0a),bytes[endpos:end])
                tail===nothing &&
                    error("no newline after \$End"*name*" in $extra")
                push = bytes[firstidx[1]:endpos+tail-1]
                base = vcat(base,push)
            end
        end
        path = joinpath(dir,"views-$key.msh")
        write(path,base)
        paths[key] = path
    end
    # A text unknown section Gmsh never wrote, appended between the mesh and
    # the data sections like a pre-processing annotation.
    for key in (:ascii,:binary)
        bytes = read(paths[key])
        marker = findfirst(codeunits("\$NodeData"),bytes)
        marker===nothing && error("fixture has no \$NodeData section")
        open(paths[key],"w") do io
            write(io,bytes[1:marker[1]-1])
            write(io,COMMENT)
            write(io,bytes[marker[1]:end])
        end
    end
    return paths
end

function probe_views(path)
    gmsh.clear()
    gmsh.open(path)
    out = Dict{String,Vector{Float64}}()
    for t in gmsh.view.getTags()
        name = gmsh.view.option.getString(t,"Name")
        out[name] = [only(gmsh.view.probe(t,x,y,z)[1]) for (x,y,z) in PROBES]
    end
    return out
end

mktempdir() do dir
    gmsh.initialize(["gmsh","-v","0"])
    gmsh.GMSH_API_VERSION=="4.15.2" || error(
        "section preservation parity requires Gmsh API 4.15.2, " *
        "got $(gmsh.GMSH_API_VERSION)")
    try
        paths = build_fixture(dir)
        reference = probe_views(paths[:ascii])
        length(reference)==3 || error("reference views missing: $reference")
        for (x,y,z) in PROBES
            i = findfirst(==((x,y,z)),PROBES)
            isapprox(reference["nodal"][i],x;atol=1e-12) ||
                error("reference nodal probe wrong")
        end

        for key in (:ascii,:binary)
            mesh = Elements.read_mixed_msh(paths[key])
            # All three views parse structurally; the comment and the
            # interpolation scheme survive as ancillary sections.
            names = [s.name for s in mesh.data_sections]
            names==["NodeData","ElementData","ElementNodeData"] ||
                error("$key: data sections lost or reordered: $names")
            anc = [s.name for s in mesh.ancillary_sections]
            ("Comments" in anc && "InterpolationScheme" in anc) ||
                error("$key: ancillary sections lost: $anc")
            mesh.data_sections[3].implicit_nodes ||
                error("$key: ElementNodeData dialect not recognised")
            # Internal element indices can be renumbered across formats, so
            # compare order-independent content; the Gmsh probes below verify
            # the tag bindings land on the right mesh locations.
            signature(s) = (s.name, length(s.elements), length(s.nodes),
                            s.implicit_nodes, sort(s.row_nodes), sort(s.values))
            expected = map(signature, mesh.data_sections)
            for (version,binary,tag) in
                ((2.2,false,"v2"),(4.1,false,"v4"),(4.1,true,"v4b"))
                out = joinpath(dir,"out-$key-$tag.msh")
                Elements.write_mixed_msh(out,mesh;version=version,binary=binary)
                again = Elements.read_mixed_msh(out)
                map(signature,again.data_sections)==expected ||
                    error("$key/$tag: data sections changed on round trip")
                findfirst(codeunits("tessella preserved comment"),
                          read(out))!==nothing ||
                    error("$key/$tag: \$Comments payload lost")
                # Gmsh itself must still read the file and bind the values
                # to the right mesh locations.
                probed = probe_views(out)
                keys(probed)==keys(reference) ||
                    error("$key/$tag: Gmsh views differ: $(keys(probed))")
                for name in keys(reference)
                    isapprox(probed[name],reference[name];atol=1e-9) ||
                        error("$key/$tag: view $name moved: " *
                              "$(probed[name]) != $(reference[name])")
                end
            end
        end
    finally
        gmsh.finalize()
    end
end
println("msh_section_preservation: all probes matched across ascii/binary " *
        "sources and v2.2/v4.1 ascii/binary outputs")
