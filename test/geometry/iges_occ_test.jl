# IGES (.igs/.iges) → native OCC kernel import: GEO `Merge`, the
# `occ.importShapes`/`ShapeFromFile` channel, MSBO solids, bounded
# faces (143/144/114/510), bare surfaces (120/122/128/198), standalone
# curves (110/100/126/102), points (116), and CSG primitives
# (150/156/158/160). Fixtures under test/fixtures cover both
# OCC-style exports (144+102+110) and the MSBO 502/504/508/510/514/186
# graph.
using Test
using Tessella
using Tessella.Model: GeoModel, _model_occ_import_shapes,
    _model_occ_synchronize!, model_entities, model_bounding_box

const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function _iges_import(name;kwargs...)
    m=GeoModel()
    roots=_model_occ_import_shapes(m,joinpath(FIXTURES,name);kwargs...)
    _model_occ_synchronize!(m)
    return m,roots
end

function _iges_merge_geo(name)
    return mktempdir() do dir
        src=joinpath(FIXTURES,name)
        dst=joinpath(dir,name)
        cp(src,dst)
        geo=joinpath(dir,"merge.geo")
        write(geo,"Merge \"$name\";\n")
        execute_geo(geo).model
    end
end

@testset "IGES OCC box faces (144+102+110)" begin
    m,roots=_iges_import("brep_occ_box.igs")
    @test roots==[(2,1),(2,2),(2,3),(2,4),(2,5),(2,6)]
    ents=model_entities(m)
    @test count(e->e[1]==0,ents)==24
    @test count(e->e[1]==1,ents)==24
    @test count(e->e[1]==2,ents)==6
end

@testset "IGES MSBO solids" begin
    m,roots=_iges_import("brep_iges_msbo_tet.igs")
    @test roots==[(3,1)]
    ents=model_entities(m)
    @test count(e->e[1]==0,ents)==4
    @test count(e->e[1]==1,ents)==6
    @test count(e->e[1]==2,ents)==4
    @test count(e->e[1]==3,ents)==1
    @test collect(model_bounding_box(m,3,1))≈[0.0,0.0,0.0,1.0,1.0,1.0]
    m,roots=_iges_import("brep_two_msbo_tets.igs")
    @test roots==[(3,1),(3,2)]
    @test count(e->e[1]==3,model_entities(m))==2
    m,roots=_iges_import("brep_planar_through_hole_msbo.igs")
    @test roots==[(3,1)]
end

@testset "IGES torus imports" begin
    m,roots=_iges_import("brep_occ_torus.igs")
    @test roots==[(2,1)]
    bb=model_bounding_box(m,2,1)
    @test bb[1]≈-2.5 atol=1e-9
    @test bb[4]≈2.5 atol=1e-9
    @test bb[3]≈-0.5 atol=1e-9
    @test bb[6]≈0.5 atol=1e-9
    m,roots=_iges_import("brep_iges_160_torus.igs")
    @test roots==[(3,1)]
    m,roots=_iges_import("brep_iges_198_torus.igs")
    @test roots==[(2,1)]
end

@testset "IGES bounded faces (143/142/190)" begin
    for (name,nfaces) in (("brep_bounded143_single_wire_box.igs",6),
                          ("brep_bounded143_through_hole.igs",10),
                          ("brep_curve142_dual_space_box.igs",6),
                          ("brep_curve142_model_space_box.igs",6),
                          ("brep_planar_through_hole_trimmed.igs",10),
                          ("brep_plane190_chain_box.igs",6),
                          ("brep_plane190_transform_exact_hole.igs",10),
                          ("brep_plane190_uv_through_hole.igs",10))
        m,roots=_iges_import(name)
        @test length(roots)==nfaces
        @test all(d==2 for (d,_) in roots)
    end
    m,roots=_iges_import("brep_plane190_uv_assembly.igs")
    @test roots==[(3,1)]   # 150 block wins highest_dim_only
end

@testset "IGES tensor/revolve/tabulated surface roots" begin
    m,roots=_iges_import("iges_occ_bspline_surface.igs")
    @test roots==[(2,1)]
    @test m.surface_types[1]==:tensor_bspline
    @test collect(model_bounding_box(m,2,1))≈[0.0,0.0,0.0,1.0,1.0,2.0] atol=1e-6
    m,roots=_iges_import("iges_occ_trimmed_bspline.igs")
    @test roots==[(2,1)]
    @test m.surface_types[1]==:tensor_bspline
    m,roots=_iges_import("iges_occ_revolve_tabulated.igs")
    @test roots==[(2,1),(2,2),(2,3)]
    # 120 full-circle generatrix about the x axis: ring torus R=2,r=.5.
    @test m.surface_types[1]==:torus
    bb=model_bounding_box(m,2,1)
    @test bb[1]≈-0.5 atol=1e-9
    @test bb[4]≈0.5 atol=1e-9
    @test bb[2]≈-2.5 atol=1e-9
    @test bb[5]≈2.5 atol=1e-9
    # 122 arc generatrix → quarter-cylinder tensor extrusion.
    bb=model_bounding_box(m,2,2)
    @test bb[2]≈3.0 atol=1e-6
    @test bb[6]≈2.0 atol=1e-6
    # 122 line generatrix → swept strip face.
    bb=model_bounding_box(m,2,3)
    @test bb[1]≈3.0 atol=1e-6
    @test bb[4]≈5.0 atol=1e-6
end

@testset "IGES standalone curves and points" begin
    m,roots=_iges_import("iges_occ_curves.igs")
    # 110 line, 100 arc, 126 spline, 102 composite → 5 dim-1 roots.
    @test roots==[(1,1),(1,2),(1,3),(1,4),(1,5)]
    m,roots=_iges_import("iges_occ_curves.igs";highest_dim_only=false)
    @test count(d->d[1]==1,roots)==5
    @test count(d->d[1]==0,roots)==1
end

@testset "IGES Merge through GEO" begin
    m=_iges_merge_geo("brep_iges_msbo_tet.igs")
    ents=model_entities(m)
    @test count(e->e[1]==3,ents)==1
    @test count(e->e[1]==2,ents)==4
    m=_iges_merge_geo("brep_occ_box.igs")
    @test count(e->e[1]==2,model_entities(m))==6
end

@testset "IGES .iges extension alias" begin
    mktempdir() do dir
        dst=joinpath(dir,"tet.iges")
        cp(joinpath(FIXTURES,"brep_iges_msbo_tet.igs"),dst)
        m=GeoModel()
        roots=_model_occ_import_shapes(m,dst)
        @test roots==[(3,1)]
    end
end

# Minimal 80-column IGES writer for one entity (S/G/D/P/T sections).
# Every record is exactly 80 columns: content in 1:72, section letter at
# 73, right-justified sequence number in 74:80.
function _mini_iges(t::Int,params::Vector)
    rec(letter,seq)=s->rpad(s,72)*letter*lpad(seq,7)
    pad8(s)=lpad(s,8)
    body=join((p isa Integer ? string(p) : string(Float64(p)) for p in params),",")*";"
    chunks=[body[i:min(i+63,end)] for i in 1:64:length(body)]
    d1=pad8(t)*pad8(1)*pad8(0)*pad8(0)*pad8(0)*pad8(0)*pad8(0)*pad8(0)*
        "00000000"
    d2=pad8(t)*pad8(0)*pad8(0)*pad8(length(chunks))*pad8(0)*pad8("")*
        pad8("")*pad8(0)
    lines=[rec("S",1)("fixture"),
           rec("G",1)("1H,,1H;"),
           rec("D",1)(d1),rec("D",2)(d2)]
    for (i,c) in enumerate(chunks)
        push!(lines,rpad(c,64)*pad8(1)*"P"*lpad(i,7))
    end
    push!(lines,rec("T",1)("S"*lpad(1,7)*"G"*lpad(1,7)*"D"*lpad(2,7)*
        "P"*lpad(length(chunks),7)))
    return join(lines,"\n")*"\n"
end

@testset "IGES unsupported-entity blockers" begin
    mktempdir() do dir
        # 104 conic arc — no supported geometry root at all.
        path=joinpath(dir,"conic.igs")
        write(path,_mini_iges(104,[104,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.]))
        m=GeoModel()
        err=try
            _model_occ_import_shapes(m,path);nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("no supported IGES geometry",sprint(showerror,err))
        # A lone 118 ruled surface is likewise outside the bound set.
        path2=joinpath(dir,"ruled.igs")
        write(path2,_mini_iges(118,[118,0,0,1,2,0,0]))
        err2=try
            _model_occ_import_shapes(m,path2);nothing
        catch e
            e
        end
        @test err2 isa ArgumentError
        @test occursin("no supported IGES geometry",sprint(showerror,err2))
        # A malformed .igs (not ASCII) rejects at the record layer.
        path3=joinpath(dir,"binary.igs")
        write(path3,"Tessella-ümlaut")
        err3=try
            _model_occ_import_shapes(m,path3);nothing
        catch e
            e
        end
        @test err3 isa ArgumentError
        @test occursin("ASCII",sprint(showerror,err3))
    end
end
