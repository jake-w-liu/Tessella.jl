# CLI UI01: `-save_parametric`/`-parametric`, `-save_topology` and the
# high-order flags (`-optimize_ho`/`-hoOptimize`/`-ho_min`/`-ho_max`/
# `-ho_nlayers`) — implemented against the verified Gmsh 4.15.2 output
# semantics; retained rejections for genuinely absent backends.
module CliSaveOptionsChecks
using Test, Tessella
using Tessella.CLI: main
using Tessella.IO: read_msh

const _SQUARE_GEO="""
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface("domain", 7) = {1};
"""

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

# Parse an ASCII `$ParametricNodes`/`$Nodes` section into per-node records.
function _parse_v2_nodes(path,section)
    lines=readlines(path)
    i=findfirst(==(section),lines)
    i===nothing && return nothing
    n=parse(Int,strip(lines[i+1]))
    return [split(strip(lines[i+1+k])) for k in 1:n]
end

# Slice one `$Start...$End` section out of an MSH file so section-scoped
# regexes cannot match records in neighbouring sections.
function _section(path,starttok,endtok)
    text=read(path,String)
    i=findfirst(starttok,text); j=findfirst(endtok,text)
    (i===nothing || j===nothing) && return ""
    return String(text[i[1]:j[end]])
end

@testset "MSH 4.1 parametric flag control" begin
    mktempdir() do dir
        geo=joinpath(dir,"sq.geo"); write(geo,_SQUARE_GEO)
        on=joinpath(dir,"on.msh"); off=joinpath(dir,"off.msh")
        # Explicit 1: parametric coordinates are written (the default output
        # already carries them).
        main([geo,"-2","-o",on,"-save_parametric"])
        main([geo,"-2","-o",off,"-setnumber","Mesh.SaveParametric","0"])
        ontext=read(on,String); offtext=read(off,String)
        @test occursin("\$Entities",ontext)
        @test occursin("\$Entities",offtext)
        onnodes=_section(on,"\$Nodes","\$EndNodes")
        offnodes=_section(off,"\$Nodes","\$EndNodes")
        # Parametric node blocks carry flag 1 and extra coordinate columns.
        @test occursin(r"^2 \d+ 1 \d+$"m,onnodes)
        @test !occursin(r"^2 \d+ 1 \d+$"m,offnodes)
        @test occursin(r"^2 \d+ 0 \d+$"m,offnodes)
        # `-parametric` is an alias of `-save_parametric`.
        aliased=joinpath(dir,"alias.msh")
        main([geo,"-2","-o",aliased,"-parametric"])
        @test occursin(r"^2 \d+ 1 \d+$"m,
            _section(aliased,"\$Nodes","\$EndNodes"))
        # Read-back preserves owners/params through the MshFile metadata.
        mf=read_msh(on)
        @test !isempty(mf.node_entities)
        @test any(p->p!==nothing && length(p)==2,mf.node_parametric)
        mfo=read_msh(off)
        @test all(p->p===nothing,mfo.node_parametric)
    end
end

@testset "MSH 2.2 ParametricNodes" begin
    mktempdir() do dir
        geo=joinpath(dir,"sq.geo"); write(geo,_SQUARE_GEO)
        out=joinpath(dir,"p.msh")
        main([geo,"-2","-o",out,"-format","msh2","-save_parametric"])
        text=read(out,String)
        @test occursin("\$ParametricNodes",text)
        @test !occursin(r"^\$Nodes$"m,text)
        records=_parse_v2_nodes(out,"\$ParametricNodes")
        @test records!==nothing && !isempty(records)
        dims=Set{Int}()
        for r in records
            length(r)>=6 || error("record $(r) too short")
            dim=parse(Int,r[5]); push!(dims,dim)
            @test parse(Int,r[6])>0
            @test length(r)-6 in (0,dim)
        end
        @test 2 in dims                      # surface-owned interior nodes
        # Read-back: owners + params land in the extended MshFile fields.
        mf=read_msh(out)
        @test length(mf.node_entities)==Tessella.MeshTypes.nnodes(mf.mesh)
        @test length(mf.node_parametric)==Tessella.MeshTypes.nnodes(mf.mesh)
        @test any(((d,t),)->d==2,mf.node_entities)
        @test any(p->p!==nothing && length(p)==2,mf.node_parametric)
    end
end

@testset "MSH 2.2 Entities topology" begin
    mktempdir() do dir
        geo=joinpath(dir,"sq.geo"); write(geo,_SQUARE_GEO)
        out=joinpath(dir,"t.msh")
        main([geo,"-2","-o",out,"-format","msh2","-save_topology"])
        text=read(out,String)
        @test occursin("\$Entities",text)
        @test occursin("\$Entities\n4 4 1 0\n",text)
        # Surface entity 1 lists its 4 boundary curves and physical tag 7
        # (upstream v2 records end with a trailing space).
        @test occursin(r"^1 4 1 2 3 4 1 7 $"m,text)
        # Points declare their coordinates-free topology record.
        @test occursin(r"^1 0 $"m,text)
        mf=read_msh(out)
        @test haskey(mf.entities,(2,1))
        @test mf.entities[(2,1)]==(Int32[1,2,3,4],Int32[7])
        @test haskey(mf.entities,(0,1))
        # Combined flags keep the upstream section order Entities→Nodes.
        both=joinpath(dir,"both.msh")
        main([geo,"-2","-o",both,"-format","msh2","-save_parametric",
              "-save_topology"])
        btext=read(both,String)
        @test occursin(r"\$Entities\n.*\$EndEntities\n\$ParametricNodes"s,btext)
        mfb=read_msh(both)
        @test haskey(mfb.entities,(2,1))
        @test !isempty(mfb.node_entities)
        # `-save_entities` is the upstream alias of `-save_topology` — it
        # stores the same `Mesh.SaveTopology` value (Gmsh 4.15.2
        # CommandLine.cpp maps both spellings to `mesh.saveTopology`).
        aliased=joinpath(dir,"entities.msh")
        main([geo,"-2","-o",aliased,"-format","msh2","-save_entities"])
        @test occursin("\$Entities",read(aliased,String))
        @test haskey(read_msh(aliased).entities,(2,1))
        # On MSH 4.1 the option is inert (entities always serialize) — the
        # alias writes the same byte stream as `-save_topology`.
        a4=joinpath(dir,"a4.msh"); t4=joinpath(dir,"t4.msh")
        main([geo,"-2","-o",a4,"-save_entities"])
        main([geo,"-2","-o",t4,"-save_topology"])
        @test read(a4)==read(t4)
    end
end

@testset "binary MSH 2.2 parametric" begin
    mktempdir() do dir
        geo=joinpath(dir,"sq.geo"); write(geo,_SQUARE_GEO)
        out=joinpath(dir,"pb.msh")
        main([geo,"-2","-o",out,"-format","msh2","-bin","-save_parametric"])
        bytes=read(out)
        @test findfirst(codeunits("\$ParametricNodes\n"),bytes)!==nothing
        @test findfirst(codeunits("\$EndParametricNodes"),bytes)!==nothing
        @test findfirst(codeunits("\$Nodes\n"),bytes)===nothing
    end
end

@testset "session path save options" begin
    mktempdir() do dir
        geo=joinpath(dir,"sq.geo"); write(geo,_SQUARE_GEO)
        # An `.msh` input routes through the API session path; `-save`
        # re-serializes the loaded model and the same file transforms apply.
        input=joinpath(dir,"in.msh")
        main([geo,"-2","-o",input,"-format","msh4","-save_parametric"])
        out=joinpath(dir,"s.msh")
        main([input,"-save","-o",out,"-format","msh2","-save_topology"])
        text=read(out,String)
        @test occursin("\$Entities",text)
        mf=read_msh(out)
        @test haskey(mf.entities,(2,1))
        # Session parametric strip on MSH 4.1 output.
        off=joinpath(dir,"soff.msh")
        main([input,"-save","-o",off,
              "-setnumber","Mesh.SaveParametric","0"])
        offtext=read(off,String)
        @test occursin("\$Entities",offtext)
        @test !occursin(r"^2 \d+ 1 \d+$"m,
            _section(off,"\$Nodes","\$EndNodes"))
        # Session MSH 2.2 parametric output.
        par=joinpath(dir,"sp.msh")
        main([input,"-save","-o",par,"-format","msh2","-save_parametric"])
        ptext=read(par,String)
        @test occursin("\$ParametricNodes",ptext)
        mfp=read_msh(par)
        @test !isempty(mfp.node_entities)
    end
end

@testset "high-order option flags" begin
    mktempdir() do dir
        geo=joinpath(dir,"sq.geo"); write(geo,_SQUARE_GEO)
        out=joinpath(dir,"ho.msh")
        # `-optimize_ho` routes through the session generation plan — an
        # order-2 triangle mesh is produced and remains valid.
        main([geo,"-2","-order","2","-optimize_ho","-o",out])
        mf=read_msh(out)
        @test Tessella.MeshTypes.ntris(mf.mesh)>0
        # `-hoOptimize` is the upstream alias.
        out2=joinpath(dir,"ho2.msh")
        main([geo,"-2","-order","2","-hoOptimize","-o",out2])
        @test isfile(out2)
        # Threshold and layer options are accepted and recorded.
        out3=joinpath(dir,"ho3.msh")
        main([geo,"-2","-order","2","-optimize_ho","-ho_min","0.2",
              "-ho_max","1.5","-ho_nlayers","4","-o",out3])
        @test isfile(out3)
        # Conflicting aliases on the same option reject.
        @test occursin("duplicate or conflicting",
            _errtext(()->main([geo,"-2","-optimize_ho","-hoOptimize"])))
        @test occursin("duplicate or conflicting",
            _errtext(()->main([geo,"-2","-save_parametric","-parametric"])))
        # Non-finite/invalid numeric arguments reject.
        @test_throws ArgumentError main([geo,"-2","-ho_min","abc"])
        @test_throws ArgumentError main([geo,"-2","-ho_nlayers","-1"])
    end
end

@testset "retained rejections" begin
    mktempdir() do dir
        geo=joinpath(dir,"sq.geo"); write(geo,_SQUARE_GEO)
        # `-0`/`-barycentric_refine` are batch actions now — combined with
        # `-2` they reject through mutual exclusion, not "not implemented".
        for flag in ("-0","-barycentric_refine")
            @test occursin("conflicting batch action",
                _errtext(()->main([geo,"-2",flag])))
        end
        # `-save_entities` is implemented (Mesh.SaveTopology alias); both
        # spellings on one option reject like the parametric pair.
        @test occursin("duplicate or conflicting",
            _errtext(()->main([geo,"-2","-save_topology","-save_entities"])))
        # Options whose only documented consumers are absent upstream
        # backends still reject — the rejection names the missing feature.
        # `-scale` is implemented (Mesh.ScalingFactor) and no longer
        # belongs on this list.
        for flag in ("-swapangle","-match")
            @test occursin("not implemented",
                _errtext(()->main([geo,"-2",flag,"1"])))
        end
    end
end

end # module
