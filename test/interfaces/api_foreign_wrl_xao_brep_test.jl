module ForeignWrlXaoBrepChecks
# Coverage for the dead-agent IO increment: `IO.read_wrl` (VRML 1.0/2.0
# indexed geometry), `IO.read_xao`/`IO.read_brep` envelope census readers,
# the `.wrl` open/merge dispatch, and the honest `.xao`/`.brep` rejection.
using Test,Tessella
const A=Tessella.API
const IO=Tessella.IO
const FIX=joinpath(@__DIR__,"..","fixtures")

function session(body)
    A.initialize()
    try;body();finally;A.finalize();end
end

@testset "read_wrl parses VRML 2.0 indexed geometry" begin
    w=IO.read_wrl(joinpath(FIX,"vrml_indexed_quad.wrl"))
    @test w isa IO.WrlFile
    @test w.version==2
    @test size(w.mesh.coords,2)==4
    @test w.node_tags==UInt64[1,2,3,4]
    # The IndexedLineSet-style 2-index record and the quad's two
    # sliding-window triangles land on a curve and a surface block.
    fams=sort!(Int32[b.msh for b in w.mesh.blocks])
    @test fams==Int32[1,2]
    counts=Dict(b.msh=>size(b.nodes,2) for b in w.mesh.blocks)
    @test counts[1]==1 && counts[2]==2
end

@testset "read_wrl parses VRML 1.0 indexed geometry" begin
    w=IO.read_wrl(joinpath(FIX,"vrml_v1_tris.wrl"))
    @test w.version==1
    @test size(w.mesh.coords,2)==4
    @test all(b.msh==2 for b in w.mesh.blocks)
    @test sum(size(b.nodes,2) for b in w.mesh.blocks)==2
end

@testset "read_wrl rejects malformed and missing files" begin
    @test_throws ArgumentError IO.read_wrl(joinpath(FIX,"does_not_exist.wrl"))
    mktempdir() do dir
        bad=joinpath(dir,"bad.wrl")
        write(bad,"#VRML V2.0 utf8\nShape { geometry IndexedFaceSet {\n"*
                 "coord Coordinate { point [ 0 0 0 ] }\n"*
                 "coordIndex [ 0 1 7 -1 ] } }\n")
        @test_throws ArgumentError IO.read_wrl(bad)
        empty_=joinpath(dir,"empty.wrl")
        write(empty_,"#VRML V2.0 utf8\n")
        @test_throws ArgumentError IO.read_wrl(empty_)
    end
end

@testset ".wrl open stages entities and cells" begin
    session() do
        A.option("General.Terminal",0)
        A.open(joinpath(FIX,"vrml_indexed_quad.wrl"))
        tags,xyz,_=A.mesh.get_nodes()
        @test length(tags)==4
        fams,_,nodes=A.mesh.get_elements(2,1)
        @test fams==Int32[2]
        @test sum(length(n) for n in nodes)÷3==2    # two triangles
        fams1,_,nodes1=A.mesh.get_elements(1,1)
        @test fams1==Int32[1]
        @test sum(length(n) for n in nodes1)÷2==1   # one line
    end
end

@testset ".wrl merges onto the current model" begin
    session() do
        A.option("General.Terminal",0)
        A.open(joinpath(FIX,"vrml_indexed_quad.wrl"))
        before=length(first(A.mesh.get_nodes()))
        A.merge(joinpath(FIX,"vrml_indexed_quad.wrl"))
        tags,_,_=A.mesh.get_nodes()
        @test length(tags)>before
    end
end

@testset "read_brep reports the CASCADE TShape census" begin
    b=IO.read_brep(joinpath(FIX,"half_cyl.brep"))
    @test b isa IO.BrepFile
    @test startswith(b.header,"CASCADE Topology V")
    @test b.tshape_total>0
    # The half-cylinder fixture: one solid, one shell, four faces/wires,
    # six edges, four vertices, one compound.
    @test b.tshapes[:solid]==1 && b.tshapes[:shell]==1
    @test b.tshapes[:face]==4 && b.tshapes[:wire]==4
    @test b.tshapes[:edge]==6 && b.tshapes[:vertex]==4
    @test b.tshapes[:compound]==1
end

@testset "read_xao reports entities, groups and the embedded census" begin
    x=IO.read_xao(joinpath(FIX,"group_box.xao"))
    @test x isa IO.XaoFile
    @test x.version=="1.0" && x.format=="BREP"
    @test !isempty(x.entities)
    @test all(e->e.dim in 0:3,x.entities)
    @test length(x.groups)==2
    @test any(g->g.name=="body",x.groups)
    @test x.brep!==nothing && x.brep.tshapes[:solid]==1 &&
          x.brep.tshapes[:face]==6 && x.brep.tshapes[:edge]==12 &&
          x.brep.tshapes[:vertex]==8
end

@testset ".xao/.brep import rejects honestly at the classified boundary" begin
    session() do
        A.option("General.Terminal",0)
        for path in (joinpath(FIX,"group_box.xao"),joinpath(FIX,"half_cyl.brep"))
            err=try
                A.open(path);nothing
            catch caught;caught;end
            @test err isa ArgumentError
            text=sprint(showerror,err)
            @test occursin("OpenCASCADE BRep boundary kernel",text)
            @test occursin("read_xao",text) || occursin("read_brep",text)
            err=try
                A.merge(path);nothing
            catch caught;caught;end
            @test err isa ArgumentError
            @test occursin("OpenCASCADE BRep boundary kernel",sprint(showerror,err))
        end
    end
end

end
