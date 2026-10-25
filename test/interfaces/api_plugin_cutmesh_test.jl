# `plugin("CutMesh")` — Gmsh 4.15.2 `CutMesh.cpp` →
# `GModel::buildCutGModel(ls, !split, saveTri)`: cuts the current mesh with
# the zero level set of the selected view, creating a `<model>_cut` model of
# discrete entities. Bounded subset: scalar level set read per mesh node
# (model NodeData views directly; list views by exact coordinate match),
# point/line/triangle elements cut; the upstream sign convention is kept —
# `ls < 0` ("inside") fragments keep the original entity tag and physicals,
# `ls > 0` ("outside") fragments move to fresh entity and physical tags with
# `_in`/`_out` name suffixes. Border segments on the zero set are emitted as
# ordinary lines on new curve entities (one per connected component;
# upstream assigns entity tag 0, which the model API does not admit).
# `Split != 0` partitions whole elements instead of subdividing.
# `SaveTri` is accepted but has no effect (as upstream, where the DI path
# always triangulates cut pieces). Element types other than 15/1/2, 3-D
# cuts, upstream's `makeSimplyConnected` regrouping and its polygon/sub-
# element records are outside the subset and rejected or omitted as
# documented.
module PluginCutMeshChecks
using Test,Tessella
const A=Tessella.API

# Unit square split into two triangles along the (1,3) diagonal, cut by the
# vertical line x = 0.5 (ls = x - 0.5: negative on the left = "inside").
function _square(name)
    A.model.add(name)
    A.model.add_discrete_entity(2,1)
    A.mesh.add_nodes(2,1,[1,2,3,4],
        [0.,0.,0.,  1.,0.,0.,  1.,1.,0.,  0.,1.,0.])
    A.mesh.add_elements_by_type(1,2,[10,11],[1,2,3, 1,3,4])
    view=A.view.add("ls")
    A.view.add_model_data(view,0,name,"NodeData",[1,2,3,4],
        [[-0.5],[0.5],[0.5],[-0.5]])
    return view
end

@testset "Plugin(CutMesh) option table" begin
    A.initialize()
    try
        _square("o")
        A.plugin.set_number("CutMesh","View",-1.)
        A.plugin.set_number("CutMesh","Split",0.)
        A.plugin.set_number("CutMesh","SaveTri",0.)
        A.plugin.set_number("CutMesh","Run",0.0)
        @test_throws ArgumentError A.plugin.set_number(
            "CutMesh","NoSuchOption",1.)
        @test_throws ArgumentError A.plugin.set_string(
            "CutMesh","Filename","x")
    finally
        A.finalize()
    end
end

@testset "Plugin(CutMesh) cuts a two-triangle square" begin
    A.initialize()
    try
        _square("sq")
        A.plugin.run("CutMesh")
        # New model named `<model>_cut`, selected as current.
        @test "sq_cut" in A.model.list()
        ents=A.model.get_entities()
        # Inside fragments keep entity tag 1; the outside half is moved to a
        # fresh tag; the border interface gets a new curve entity.
        @test (2,1) in ents
        @test (1,1) in ents
        outside=sort([t for (d,t) in ents if d==2 && t!=1])
        @test length(outside)==1
        # Inside half: the lone left corner + the quad-side of triangle 2.
        etypes,etags,enodes=A.mesh.get_elements(2,1)
        @test etypes==[2]
        @test length(only(etags))==3
        # Outside half: the corner + quad-side of triangle 1.
        etypes2,etags2,enodes2=A.mesh.get_elements(2,only(outside))
        @test etypes2==[2]
        @test length(only(etags2))==3
        # Border entity carries the two interface segments.
        btypes,btags,bnodes=A.mesh.get_elements(1,1)
        @test btypes==[1]
        @test length(only(btags))==2
        # All element tags are unique across entities.
        alltags=collect(Iterators.flatten(
            (etags[1],etags2[1],btags[1])))
        @test length(alltags)==length(unique(alltags))
        # The three new cut vertices lie on x = 0.5, shared by both sides.
        ntags,ncoords,_=A.mesh.get_nodes(1,1)
        @test length(ntags)==3
        for i in 1:length(ntags)
            @test ncoords[(i-1)*3+1]≈0.5 atol=1e-12
        end
        @test sort([ncoords[(i-1)*3+2] for i in 1:3])≈[0.0,0.5,1.0] atol=1e-12
    finally
        A.finalize()
    end
end

@testset "Plugin(CutMesh) physical group _in/_out split" begin
    A.initialize()
    try
        A.model.add("pg")
        A.model.add_discrete_entity(2,1)
        A.mesh.add_nodes(2,1,[1,2,3,4],
            [0.,0.,0.,  1.,0.,0.,  1.,1.,0.,  0.,1.,0.])
        A.mesh.add_elements_by_type(1,2,[10,11],[1,2,3, 1,3,4])
        A.model.add_physical_group(2,[1];tag=5,name="mat")
        view=A.view.add("ls")
        A.view.add_model_data(view,0,"pg","NodeData",[1,2,3,4],
            [[-0.5],[0.5],[0.5],[-0.5]])
        A.plugin.run("CutMesh")
        groups=A.model.get_physical_groups()
        @test (2,5) in groups
        # Inside entities keep the source physical renamed *_in; the
        # outside physical is a fresh tag renamed *_out.
        @test A.model.get_physical_name(2,5)=="mat_in"
        @test A.model.get_entities_for_physical_group(2,5)==[1]
        out=only(setdiff(groups,[(2,5)]))
        @test out[1]==2
        @test A.model.get_physical_name(out[1],out[2])=="mat_out"
    finally
        A.finalize()
    end
end

@testset "Plugin(CutMesh) Split partitions whole elements" begin
    A.initialize()
    try
        _square("sp")
        A.plugin.set_number("CutMesh","Split",1.0)
        A.plugin.run("CutMesh")
        ents=A.model.get_entities()
        @test "sp_cut" in A.model.list()
        # Both triangles are mixed-sign → all land on a fresh outside
        # entity, unsplit.
        etypes,etags,_=A.mesh.get_elements(2,2)
        @test etypes==[2]
        @test length(only(etags))==2
        # One border line: the single fully-inside edge (4,1).
        bents=[t for (d,t) in ents if d==1]
        @test length(bents)==1
        btypes,btags,bnodes=A.mesh.get_elements(1,only(bents))
        @test btypes==[1]
        flat=only(bnodes)
        @test Set(minmax(flat[1],flat[2]))==Set([UInt64(1),UInt64(4)])
    finally
        A.finalize()
    end
end

@testset "Plugin(CutMesh) point and line elements" begin
    A.initialize()
    try
        A.model.add("pl")
        A.model.add_discrete_entity(0,1)
        A.model.add_discrete_entity(0,2)
        A.model.add_discrete_entity(1,1)
        A.mesh.add_nodes(0,1,[1],[0.,0.,0.])
        A.mesh.add_nodes(0,2,[2],[1.,0.,0.])
        A.mesh.add_nodes(1,1,[3,4],[0.,1.,0., 1.,1.,0.])
        A.mesh.add_elements_by_type(1,15,[20],[1])
        A.mesh.add_elements_by_type(2,15,[21],[2])
        A.mesh.add_elements_by_type(1,1,[30],[3,4])
        view=A.view.add("ls")
        # Node 2 (x=1) outside; node 1 exactly on the level set; the line
        # crosses at x=0.5.
        A.view.add_model_data(view,0,"pl","NodeData",[1,2,3,4],
            [[0.0],[0.5],[-0.5],[0.5]])
        A.plugin.run("CutMesh")
        ents=A.model.get_entities()
        # Point 1 sits on the level set → keeps entity (0,1) and stays
        # inside; point 2 is outside → fresh point entity.
        pt_out=sort([t for (d,t) in ents if d==0 && t!=1])
        @test length(pt_out)==1
        etypes,_,_=A.mesh.get_elements(0,only(pt_out))
        @test etypes==[15]
        # The line is cut: inside half on entity (1,1), outside half on a
        # fresh curve entity.
        ltypes,ltags,_=A.mesh.get_elements(1,1)
        @test ltypes==[1]
        @test length(only(ltags))==1
        lout=[t for (d,t) in ents if d==1 && t!=1]
        @test length(lout)==1
        ltypes2,ltags2,_=A.mesh.get_elements(1,only(lout))
        @test ltypes2==[1]
        @test length(only(ltags2))==1
    finally
        A.finalize()
    end
end

@testset "Plugin(CutMesh) input validation" begin
    A.initialize()
    try
        A.model.add("err")
        A.model.add_discrete_entity(2,1)
        A.mesh.add_nodes(2,1,[1,2,3],
            [0.,0.,0.,  1.,0.,0.,  0.,1.,0.])
        A.mesh.add_elements_by_type(1,2,[10],[1,2,3])
        # No views exist: View=-1 cannot resolve.
        @test_throws ArgumentError A.plugin.run("CutMesh")
        # A NodeData view missing a mesh node errors clearly.
        view=A.view.add("ls")
        A.view.add_model_data(view,0,"err","NodeData",[1,2],
            [[-0.5],[0.5]])
        @test_throws ArgumentError A.plugin.run("CutMesh")
        # A non-scalar view is rejected.
        view2=A.view.add("ls3")
        A.view.add_model_data(view2,0,"err","NodeData",[1,2,3],
            [[-0.5,0.,0.],[0.5,0.,0.],[0.5,0.,0.]])
        A.plugin.set_number("CutMesh","View",view2)
        @test_throws ArgumentError A.plugin.run("CutMesh")
    finally
        A.finalize()
    end
end
end
