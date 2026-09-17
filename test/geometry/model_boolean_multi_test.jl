using Test
using Tessella
using Tessella.MeshTypes: ntets, node, tet_volume, validate
using Tessella.Model: model_bounding_box

# All expected tags/orderings below were verified against Gmsh 4.15.2
# (gmsh.model.occ boolean calls with removeObject=removeTool=true).

_xbounds(m,t)=begin
    b=model_bounding_box(m,3,t)
    (round(b[1];digits=9),round(b[4];digits=9))
end

function _exec_geo(src)
    return mktemp() do path,io
        write(io,src);close(io)
        execute_geo(path)
    end
end

@testset "N-way Boolean operands — OCC cell semantics" begin
    # difference: objects [0,1.5] and [0.75,2.25], tool [1.2,2.7] —
    # kept cells {1}=[0,0.75] (fresh) and {1,2}=[0.75,1.2] (rebinds tag 20,
    # operand 20's unique sole image); Gmsh: out=[(3,21),(3,20)]
    m=GeoModel()
    add_box!(m,0,0,0,1.5,1,1;tag=10)
    add_box!(m,0.75,0,0,1.5,1,1;tag=20)
    add_box!(m,1.2,0,0,1.5,1,1;tag=30)
    out=boolean_volumes_multi!(m,:difference,[10,20],[30];
                             remove_object=true,remove_tool=true)
    @test out==[21,20]
    @test sort!(collect(keys(m.volumes)))==[20,21]
    @test _xbounds(m,21)==(0.0,0.75)
    @test _xbounds(m,20)==(0.75,1.2)

    # fragments over three chained boxes — five cells in mask order;
    # Gmsh: out=[(3,1)..(3,5)] with the same piece spans
    m=GeoModel()
    add_box!(m,0,0,0,1.5,1,1;tag=10)
    add_box!(m,0.75,0,0,1.5,1,1;tag=20)
    add_box!(m,1.2,0,0,1.5,1,1;tag=30)
    out=boolean_volumes_multi!(m,:fragments,[10,20],[30];
                             remove_object=true,remove_tool=true)
    @test out==[1,2,3,4,5]
    @test [_xbounds(m,t) for t in out]==
          [(0.0,0.75),(0.75,1.2),(1.2,1.5),(1.5,2.25),(2.25,2.7)]

    # intersection objects {10,20} ∩ tool {30} — kept cells {1,2,3} (rebinds
    # tag 10, operand 10's unique sole image) and {2,3} (fresh);
    # Gmsh: out=[(3,10),(3,11)]
    m=GeoModel()
    add_box!(m,0,0,0,1.5,1,1;tag=10)
    add_box!(m,0.75,0,0,1.5,1,1;tag=20)
    add_box!(m,1.2,0,0,1.5,1,1;tag=30)
    out=boolean_volumes_multi!(m,:intersection,[10,20],[30];
                             remove_object=true,remove_tool=true)
    @test out==[10,11]
    @test _xbounds(m,10)==(1.2,1.5)
    @test _xbounds(m,11)==(1.5,2.25)

    # two-box fragments — three cells; Gmsh: out=[1,2,3], spans as below
    m=GeoModel()
    add_box!(m,0,0,0,1.5,1,1;tag=10)
    add_box!(m,1,0,0,1.5,1,1;tag=20)
    out=boolean_volumes_multi!(m,:fragments,[10],[20];
                             remove_object=true,remove_tool=true)
    @test out==[1,2,3]
    @test [_xbounds(m,t) for t in out]==
          [(0.0,1.0),(1.0,1.5),(1.5,2.5)]
    # the shared wall is a single surface referenced by both shells
    @test length(m.surfaces)==16
    @test length(m.curves)==28

    # disjoint tool is untouched: under fragments it rebinds its own tag;
    # Gmsh: out=[(3,31),(3,32),(3,33),(3,30)]
    m=GeoModel()
    add_box!(m,0,0,0,1,1,1;tag=10)
    add_box!(m,0.5,0,0,1,1,1;tag=20)
    add_box!(m,2,0,0,1,1,1;tag=30)
    out=boolean_volumes_multi!(m,:fragments,[10,20],[30];
                             remove_object=true,remove_tool=true)
    @test out==[31,32,33,30]
    @test _xbounds(m,30)==(2.0,3.0)

    # a disjoint tool under difference has no image — Delete removes it and
    # the three object cells take fresh tags; Gmsh: out=[(3,1),(3,2),(3,3)]
    m=GeoModel()
    add_box!(m,0,0,0,1,1,1;tag=10)
    add_box!(m,0.5,0,0,1,1,1;tag=20)
    add_box!(m,2,0,0,1,1,1;tag=30)
    out=boolean_volumes_multi!(m,:difference,[10,20],[30];
                             remove_object=true,remove_tool=true)
    @test out==[1,2,3]
    @test sort!(collect(keys(m.volumes)))==[1,2,3]
end

@testset "N-way Boolean containment — preserved operands" begin
    # contained operand is preserved under intersection — OCC IsSame;
    # Gmsh: out=[(3,2)], the inner operand keeps its tag
    m=GeoModel()
    add_box!(m,0,0,0,2,2,2;tag=1)
    add_box!(m,0.5,0.5,0.5,1,1,1;tag=2)
    out=boolean_volumes_multi!(m,:intersection,[1],[2];
                             remove_object=true,remove_tool=true)
    @test out==[2]
    @test sort!(collect(keys(m.volumes)))==[2]
    @test _xbounds(m,2)==(0.5,1.5)

    # containment union collapses to the outer operand; Gmsh: out=[(3,1)]
    m=GeoModel()
    add_box!(m,0,0,0,2,2,2;tag=1)
    add_box!(m,0.5,0.5,0.5,1,1,1;tag=2)
    out=boolean_volumes_multi!(m,:union,[1],[2];
                             remove_object=true,remove_tool=true)
    @test out==[1]
    @test sort!(collect(keys(m.volumes)))==[1]

    # containment fragments — the cavity wall is shared between the inner
    # solid and the outer shell; Gmsh (tags 10,20): out=[(3,21),(3,20)]
    m=GeoModel()
    add_box!(m,0,0,0,2,2,2;tag=10)
    add_box!(m,0.5,0.5,0.5,1,1,1;tag=20)
    out=boolean_volumes_multi!(m,:fragments,[10],[20];
                             remove_object=true,remove_tool=true)
    @test out==[21,20]
    @test sort!(collect(keys(m.volumes)))==[20,21]
    @test length(m.points)==16
    @test length(m.curves)==24
    @test length(m.surfaces)==12   # 6 outer + 6 shared inner faces
    outer_shells=[m.surface_loops[sh] for sh in m.volumes[21]]
    inner_surfs=Set(abs.(m.surface_loops[only(m.volumes[20])]))
    outer_all=Set(s for sh in outer_shells for s in abs.(sh))
    @test length(outer_shells)==2      # outer shell + cavity shell
    @test inner_surfs ⊆ outer_all      # cavity wall shares the inner faces

    # a fused piece under an untouched disjoint tool — the tool is preserved
    # in place; Gmsh (tags 10,20,30): out=[(3,31),(3,30)]
    m=GeoModel()
    add_box!(m,0,0,0,1.5,1,1;tag=10)
    add_box!(m,1,0,0,1.5,1,1;tag=20)
    add_box!(m,3,0,0,1,1,1;tag=30)
    out=boolean_volumes_multi!(m,:union,[10,20],[30];
                             remove_object=true,remove_tool=true)
    @test out==[31,30]
    @test _xbounds(m,30)==(3.0,4.0)
    @test _xbounds(m,31)==(0.0,2.5)
end

@testset "N-way Boolean multi-component cells and misc" begin
    # a slicing tool splits the object into two cells with the same
    # membership mask — each piece extracts its own mesh component
    m=GeoModel()
    add_box!(m,0,0,0,1,1,1;tag=1)
    add_box!(m,0.4,-0.5,-0.5,0.2,2,2;tag=2)
    out=boolean_volumes_multi!(m,:difference,[1],[2];
                             remove_object=true,remove_tool=true)
    @test out==[1,2]
    @test _xbounds(m,1)==(0.0,0.4)
    @test _xbounds(m,2)==(0.6,1.0)
    for t in out
        mesh=mesh_model_volume(m,t)
        @test ntets(mesh)>0
        @test validate(mesh).ok
    end
    translate_volume!(m,1,(5.0,0.0,0.0))
    @test _xbounds(m,1)==(5.0,5.4)
    @test _xbounds(m,2)==(0.6,1.0)

    # an explicit tag on a multi-piece result is an error (Gmsh: "could not
    # apply boolean operator" / single tag expected)
    m=GeoModel()
    add_box!(m,0,0,0,1,1,1;tag=1)
    add_box!(m,0.4,-0.5,-0.5,0.2,2,2;tag=2)
    @test_throws ArgumentError boolean_volumes_multi!(
        m,:difference,[1],[2];tag=9,remove_object=true,remove_tool=true)
    @test haskey(m.volumes,1) && haskey(m.volumes,2)

    # empty result — a fully-swallowed object under Delete removes both
    # operands and binds nothing; Gmsh 4.15.2: out=[], no volumes remain
    m=GeoModel()
    add_box!(m,0,0,0,1,1,1;tag=1)
    add_box!(m,-1,-1,-1,3,3,3;tag=2)
    out=boolean_volumes_multi!(m,:difference,[1],[2];
                             remove_object=true,remove_tool=true)
    @test isempty(out)
    @test isempty(m.volumes)

    # a disjoint tool under difference leaves the object untouched — the
    # operand is preserved with its tag; Gmsh: out=[(3,10)], vols=[10]
    m=GeoModel()
    add_box!(m,0,0,0,1,1,1;tag=10)
    add_box!(m,5,0,0,1,1,1;tag=20)
    out=boolean_volumes_multi!(m,:difference,[10],[20];
                             remove_object=true,remove_tool=true)
    @test out==[10]
    @test sort!(collect(keys(m.volumes)))==[10]

    # operand cap — UInt64 membership masks support at most 62 operands
    m=GeoModel()
    for i in 1:63
        add_box!(m,10.0*i,0,0,1,1,1;tag=i)
    end
    @test_throws ArgumentError boolean_volumes_multi!(
        m,:union,collect(1:63),Int[];remove_object=true,remove_tool=true)

    # unsupported operand kind fails explicitly
    m=GeoModel()
    add_box!(m,0,0,0,1,1,1;tag=1)
    add_point!(m,0,0,0;tag=50)
    @test_throws ArgumentError boolean_volumes_multi!(
        m,:union,[1],[50];remove_object=true,remove_tool=true)

    # materialized result vertices carry no explicit mesh-size constraint —
    # like OCC primitive materialization, a bogus size would feed
    # `_volume_boundary_size_field` and silently refine the volume mesh
    # (oracle: the unit difference mesh has 12 tets and volume exactly 1.0)
    m=GeoModel()
    add_box!(m,0,0,0,2,1,1;tag=1)
    add_box!(m,0,0,0,1,1,1;tag=2)
    boolean_volumes!(m,:difference,1,2;tag=3)
    @test isempty(m.point_size)
    mesh=mesh_model_volume(m,3)
    @test ntets(mesh)==12
    @test sum(tet_volume(node(mesh,mesh.tets[1,t]),
                         node(mesh,mesh.tets[2,t]),
                         node(mesh,mesh.tets[3,t]),
                         node(mesh,mesh.tets[4,t]))
              for t in 1:ntets(mesh))==1.0
end

@testset ".geo N-way Boolean statements" begin
    # v() capture + multi-operand object list + disjoint tool
    parsed=_exec_geo("""
    Box(1) = {0,0,0, 1,1,1};
    Box(2) = {0.5,0,0, 1,1,1};
    Box(3) = {2,0,0, 1,1,1};
    v() = BooleanDifference{ Volume{1,2}; Delete; }{ Volume{3}; Delete; };
    Physical Volume(50) = v[];
    """)
    @test sort!(collect(keys(parsed.model.volumes)))==[1,2,3]
    @test parsed.model.physical[(3,50)]==[1,2,3]

    # name[] capture of BooleanFragments over a two-brace form
    parsed=_exec_geo("""
    Box(1) = {0,0,0, 1.5,1,1};
    Box(2) = {1,0,0, 1.5,1,1};
    w[] = BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; };
    Physical Volume(51) = w[];
    """)
    @test sort!(collect(keys(parsed.model.volumes)))==[1,2,3]
    @test parsed.model.physical[(3,51)]==[1,2,3]

    # tagged form with a multi-operand list — single fused piece
    parsed=_exec_geo("""
    Box(1) = {0,0,0, 1,1,1};
    Box(2) = {0.5,0,0, 1,1,1};
    BooleanUnion(7) = { Volume{1,2}; Delete; }{ };
    """)
    @test sort!(collect(keys(parsed.model.volumes)))==[7]
    @test _xbounds(parsed.model,7)==(0.0,1.5)

    # an empty tool group is legal for fragments/union (union of the objects)
    parsed=_exec_geo("""
    Box(1) = {0,0,0, 1,1,1};
    Box(2) = {0.5,0,0, 1,1,1};
    v() = BooleanFragments{ Volume{1,2}; Delete; }{ };
    """)
    @test sort!(collect(keys(parsed.model.volumes)))==[1,2,3]

    # single-brace form is a syntax error in Gmsh too — rejected
    @test_throws ArgumentError _exec_geo("""
    Box(1) = {0,0,0, 1,1,1};
    Box(2) = {0.5,0,0, 1,1,1};
    v() = BooleanUnion{ Volume{1,2}; Delete; };
    """)

    # explicit tag on a multi-piece result errors inside .geo as well
    err=try
        _exec_geo("""
        Box(1) = {0,0,0, 1,1,1};
        Box(2) = {0.4,-0.5,-0.5, 0.2,2,2};
        BooleanDifference(9) = { Volume{1}; Delete; }{ Volume{2}; Delete; };
        """)
        nothing
    catch e
        e
    end
    @test err!==nothing
    @test occursin("2 result volumes", sprint(showerror,err))

    # tagged form with multi-operand groups in BOTH braces — a single fused
    # piece binds the requested tag; a later Point(newp) resolves against the
    # materialized boundary (the tracker cannot model the statement, so the
    # allocator resyncs to the model)
    parsed=_exec_geo("""
    Box(1) = {0,0,0, 1,1,1};
    Box(2) = {0.5,0,0, 1,1,1};
    Box(3) = {0.75,0,0, 1,1,1};
    BooleanUnion(9) = { Volume{1,2}; Delete; }{ Volume{3}; Delete; };
    Point(newp) = {9,9,9,1};
    """)
    @test sort!(collect(keys(parsed.model.volumes)))==[9]
    @test _xbounds(parsed.model,9)==(0.0,1.75)
    @test parsed.model.points[maximum(keys(parsed.model.points))]==
        (9.0,9.0,9.0)

    # a disjoint tool under a multi-object difference leaves three cells —
    # v() capture then Physical Volume consumes the out list
    parsed=_exec_geo("""
    Box(1) = {0,0,0, 1,1,1};
    Box(2) = {0.5,0,0, 1,1,1};
    Box(3) = {5,0,0, 1,1,1};
    v() = BooleanDifference{ Volume{1,2}; Delete; }{ Volume{3}; Delete; };
    Physical Volume(52) = v[];
    """)
    @test sort!(collect(keys(parsed.model.volumes)))==[1,2,3]
    @test parsed.model.physical[(3,52)]==[1,2,3]
end
