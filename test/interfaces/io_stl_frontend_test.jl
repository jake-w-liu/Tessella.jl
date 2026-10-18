module IoStlFrontendChecks
# Focused checks for the `.stl` session lifecycle, the versioned `.mshN`
# write aliases, `Extrude {shapes} Using Wire {expr}` (OCC pipe), and
# BoundaryLayer field resolution in `.geo` meshing.
using Test, Tessella
const A=Tessella.API

const _TRI_STL="""
solid tile
 facet normal 0 0 1
  outer loop
   vertex 0 0 0
   vertex 1 0 0
   vertex 0 1 0
  endloop
 endfacet
 facet normal 0 0 1
  outer loop
   vertex 1 0 0
   vertex 1 1 0
   vertex 0 1 0
  endloop
 endfacet
endsolid tile
"""

function _write_binary_stl(path)
    open(path,"w") do io
        write(io,zeros(UInt8,80))
        write(io,UInt32(1))
        for v in (0f0,0f0,1f0, 0f0,0f0,0f0, 1f0,0f0,0f0, 0f0,1f0,0f0)
            write(io,Float32(v))
        end
        write(io,UInt16(0))
    end
    return path
end

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

_geo_path(path)=replace(path,'\\'=>'/')

function _execute_geo_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source);close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _geo_error(source::AbstractString;mesh_dim=0)
    try
        _execute_geo_source(source;mesh_dim)
        return nothing
    catch err
        return err
    end
end

@testset ".stl session lifecycle (open/merge/write)" begin
    A.initialize()
    try
        dir=mktempdir()
        stl=joinpath(dir,"tile.stl");write(stl,_TRI_STL)
        A.open(stl)
        @test (2,1) in A.model.get_entities()
        @test A.model.get_type(2,1)=="Discrete surface"
        types,etags,ntags=A.mesh.get_elements()
        @test Set(types)==Set([2])
        @test length(reduce(vcat,etags;init=eltype(first(etags))[]))==2
        # Re-writing the session mesh emits ASCII STL facets.
        out=joinpath(dir,"copy.stl")
        A.write(out)
        @test isfile(out)
        text=read(out,String)
        @test occursin("solid",text) && occursin("endsolid",text)
        @test length(collect(eachmatch(r"facet normal",text)))==2
        # Merging the same STL lands on the lowest free discrete surface.
        A.merge(stl)
        @test (2,2) in A.model.get_entities()
        @test A.model.get_type(2,2)=="Discrete surface"
        types2,etags2,_=A.mesh.get_elements()
        @test length(reduce(vcat,etags2;init=eltype(first(etags2))[]))==4
        # Binary STL reads through the same channel.
        binary=_write_binary_stl(joinpath(dir,"one.stl"))
        A.merge(binary)
        @test (2,3) in A.model.get_entities()
        # STL export rejects volume cells and never creates the destination.
        A.model.add("vol")
        A.model.add_discrete_entity(3,1)
        A.mesh.add_nodes(3,1,[1,2,3,4],[0.,0.,0., 1.,0.,0., 0.,1.,0., 0.,0.,1.])
        A.mesh.add_elements_by_type(1,4,[7],[1,2,3,4])
        bad=joinpath(dir,"should_not_exist.stl")
        err=_errtext(()->A.write(bad))
        @test occursin("volume",err) && occursin("stl",lowercase(err))
        @test !isfile(bad)
    finally
        A.finalize()
    end
end

@testset ".mshN write aliases" begin
    A.initialize()
    try
        dir=mktempdir()
        stl=joinpath(dir,"tile.stl");write(stl,_TRI_STL)
        A.open(stl)
        for (ext,want) in (("msh2","2.2"),("msh22","2.2"),
                           ("msh4","4.1"),("msh41","4.1"),("msh","4.1"))
            path=joinpath(dir,"out.$ext")
            A.write(path)
            @test isfile(path)
            head=readlines(path)
            @test head[1]==raw"$MeshFormat"
            @test startswith(head[2],want)
        end
        for ext in ("msh1","msh3","msh40")
            path=joinpath(dir,"out.$ext")
            err=_errtext(()->A.write(path))
            @test occursin("unsupported MSH extension",err)
            @test !isfile(path)
        end
    finally
        A.finalize()
    end
end

@testset "GEO Include/Merge/Save .stl and .mshN" begin
    dir=mktempdir()
    stl=joinpath(dir,"tile.stl");write(stl,_TRI_STL)
    out_stl=joinpath(dir,"saved.stl")
    out_msh=joinpath(dir,"saved.msh41")
    r=_execute_geo_source("""
        Merge "$(_geo_path(stl))";
        Save "$(_geo_path(out_stl))";
        Save "$(_geo_path(out_msh))";
        """)
    @test haskey(r.model.meshing.imported_entities,(2,1))
    @test r.model.meshing.imported_entities[(2,1)].kind=="Discrete surface"
    @test isfile(out_stl) && occursin("facet normal",read(out_stl,String))
    @test isfile(out_msh)
    head=readlines(out_msh)
    @test head[1]==raw"$MeshFormat" && startswith(head[2],"4.1")
    err=_geo_error("""
        Save "$(_geo_path(joinpath(dir,"bad.msh3")))";
        """)
    @test err isa ArgumentError
    @test occursin("unsupported MSH extension",sprint(showerror,err))
    @test !isfile(joinpath(dir,"bad.msh3"))
end

const _PIPE_SQUARE="""
Point(1) = {0,0,0,1};
Point(2) = {1,0,0,1};
Point(3) = {1,1,0,1};
Point(4) = {0,1,0,1};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Point(5) = {0.5,0.5,0,1};
Point(6) = {0.5,0.5,1,1};
Line(5) = {5,6};
Curve Loop(2) = {5};
"""

@testset "OCC `Extrude {shapes} Using Wire {expr}`" begin
    err=_geo_error(_PIPE_SQUARE * """
        Extrude { Surface{1}; } Using Wire {2}
        """)
    @test err isa ArgumentError
    @test occursin("pipe",sprint(showerror,err))

    r=_execute_geo_source("SetFactory(\"OpenCASCADE\");\n" * _PIPE_SQUARE * """
        v[] = Extrude { Surface{1}; } Using Wire {2};
        """)
    # The unit-square profile swept one unit up produces a cube volume plus
    # lateral faces, all published back into the GEO records.
    @test length(r.model.volumes)>=1
    @test length(r.model.surfaces)>=5
    vtags=Int.(r.lists["v"])
    @test any(t->haskey(r.model.volumes,t),vtags)

    # A `Curve{}` profile produces a piped surface rather than a solid.
    r2=_execute_geo_source("SetFactory(\"OpenCASCADE\");\n" * _PIPE_SQUARE * """
        s[] = Extrude { Curve{5}; } Using Wire {2};
        """)
    stags=Int.(r2.lists["s"])
    @test any(t->haskey(r2.model.surfaces,t),stags)
end

@testset "BoundaryLayer field resolution in .geo meshing" begin
    square="""
    Point(1) = {0,0,0,0.5};
    Point(2) = {1,0,0,0.5};
    Point(3) = {1,1,0,0.5};
    Point(4) = {0,1,0,0.5};
    Line(1) = {1,2};
    Line(2) = {2,3};
    Line(3) = {3,4};
    Line(4) = {4,1};
    Curve Loop(1) = {1,2,3,4};
    Plane Surface(1) = {1};
    """
    # A BoundaryLayer field used as `Background Field` resolves through the
    # model context (wall distance -> hwall/ratio/thickness/hfar growth) and
    # genuinely feeds the mesher: the curves grade to ~0.08 near the walls,
    # far denser than the 0.5 point lc.
    r=_execute_geo_source(square * """
        Field[1] = BoundaryLayer;
        Field[1].CurvesList = {1,2,3,4};
        Field[1].Size = 0.08;
        Field[1].Ratio = 1.4;
        Field[1].Thickness = 0.5;
        Field[1].SizeFar = 0.3;
        Background Field = 1;
        Mesh 1;
        """)
    @test r.mesh!==nothing && size(r.mesh.coords,2)>=20

    # `BoundaryLayer Field =` selects layer generation; without a layer
    # consumer the selection must reject loudly instead of meshing on.
    err=_geo_error(square * """
        Field[1] = BoundaryLayer;
        Field[1].CurvesList = {1,2,3,4};
        BoundaryLayer Field = 1;
        Mesh 2;
        """)
    @test err isa ArgumentError
    @test occursin("boundary-layer",sprint(showerror,err))

    err=_geo_error(square * """
        Field[1] = MathEval;
        Field[1].F = "0.2";
        BoundaryLayer Field = 1;
        Mesh 2;
        """)
    @test err isa ArgumentError
    @test occursin("not BoundaryLayer",sprint(showerror,err))

    err=_geo_error(square * """
        BoundaryLayer Field = 9;
        Mesh 2;
        """)
    @test err isa ArgumentError
    @test occursin("not declared",sprint(showerror,err))
end

end # module
