using Test
using Tessella
using Tessella.Elements: ElementBlock, MixedMesh, msh_dimension, lagrange_nodes
using Tessella.MeshTypes: nnodes, validate

function _api_refine_install(mesh)
    api=Tessella.API
    api.initialize()
    dim=maximum(msh_dimension(b.msh) for b in mesh.blocks)
    api.model.add_discrete_entity(dim,7)
    owners=Dict(Int(b.msh)=>fill(Int32(7),size(b.nodes,2)) for b in mesh.blocks)
    class=api._mixed_classification(mesh,(dim,Int32(7)),[(dim,Int32(7))],
        fill((dim,Int32(7)),nnodes(mesh)),
        Dict((dim,Int32(7))=>Int32[]),owners)
    lock(api.STATE_LOCK) do
        api._replace_mesh_cache_locked!(mesh,class)
    end
    return mesh
end

function _api_refine_single_family(msh;coords=lagrange_nodes(msh))
    cells=reshape(Int32.(1:Tessella.Elements.msh_num_nodes(msh)),:,1)
    mesh=MixedMesh(coords,[ElementBlock(msh,cells,Int32[17])];
        physical_names=Dict((msh_dimension(msh),17)=>"material"))
    return _api_refine_install(mesh)
end

function _api_refine_volume(mesh)
    return sum(Tessella.Model._extrude_cell_signed_volume(mesh.coords,b.nodes,column,b.msh)
        for b in mesh.blocks if b.msh in (4,5,6,7) for column in axes(b.nodes,2);init=0.0)
end

_api_refine_counts(mesh)=Dict(Int(b.msh)=>size(b.nodes,2) for b in mesh.blocks)

@testset "native mixed API refinement" begin
    api=Tessella.API
    for (msh,nnode,counts,volume) in (
            (1,3,Dict(1=>2),0.0),(2,6,Dict(2=>4),0.0),
            (3,9,Dict(3=>4),0.0),(4,10,Dict(4=>8),1/6),
            (5,27,Dict(5=>8),8.0),(6,18,Dict(6=>8),1.0),
            (7,14,Dict(4=>8,7=>4),4/3),(15,1,Dict(15=>1),0.0))
        try
            _api_refine_single_family(msh)
            refined=api.mesh.refine(max_nodes=nnode,max_cells=sum(values(counts)))
            @test refined isa MixedMesh
            @test nnodes(refined)==nnode
            @test _api_refine_counts(refined)==counts
            @test validate(refined).ok
            @test _api_refine_volume(refined)≈volume
            @test all(b->all(==(17),b.tags),refined.blocks)
            dim=msh_dimension(msh)
            @test api.mesh.get_element_types(dim)==sort!(Int32.(collect(keys(counts))))
            types,tags,nodes=api.mesh.get_elements(dim)
            @test sum(length,tags;init=0)==sum(values(counts))
            @test all(tag->api.mesh.get_element(tag)[3:4]==(dim,7),vcat(tags...))
            @test all(==((dim,Int32(7))),api.LAST_MESH_CLASS[].node_entities)
            before=copy(api.mesh.get().coords)
            refined.coords.=19
            @test api.mesh.get().coords==before
        finally
            api.finalize()
        end
    end
    @testset "preflight failures leave the cache and metadata untouched" begin
        try
            _api_refine_single_family(5)
            before=api.LAST_MESH[];class=api.LAST_MESH_CLASS[]
            api.mesh.create_edges()
            edges=api.LAST_MESH_EDGES[]
            for keywords in ((;max_nodes=26),(;max_cells=7),(;max_nodes=true),
                    (;max_cells=-1),(;max_nodes=1.5),(;max_cells=big(typemax(Int32))+1))
                @test_throws ArgumentError api.mesh.refine(;keywords...)
                @test api.LAST_MESH[]===before
                @test api.LAST_MESH_CLASS[]===class
                @test api.LAST_MESH_EDGES[]===edges
            end
            bad=copy(lagrange_nodes(1));bad[:,2]=bad[:,1]
            api.finalize();_api_refine_single_family(1;coords=bad)
            before=api.LAST_MESH[];class=api.LAST_MESH_CLASS[]
            @test_throws r"coordinate resolution" api.mesh.refine()
            @test api.LAST_MESH[]===before
            @test api.LAST_MESH_CLASS[]===class
        finally
            api.finalize()
        end
    end
    @testset "repeat refinement and quadratic input" begin
        for (msh,nnode,counts,volume) in (
                (5,125,Dict(5=>64),8.0),(6,75,Dict(6=>64),1.0),
                (7,55,Dict(4=>96,7=>16),4/3))
            try
                _api_refine_single_family(msh)
                api.mesh.set_order(2)
                @test api.mesh.get_element_types()==Int32[msh==5 ? 12 : msh==6 ? 13 : 14]
                first=api.mesh.refine()
                @test all(b->Tessella.Elements.msh_spec(b.msh).order==1,first.blocks)
                second=api.mesh.refine()
                @test nnodes(second)==nnode
                @test _api_refine_counts(second)==counts
                @test _api_refine_volume(second)≈volume
                @test validate(second).ok
            finally
                api.finalize()
            end
        end
    end
    @testset "generated quadrangles share graded curve nodes and owners" begin
        try
            api.initialize()
            for (x,y) in ((0,0),(1,0),(1,1),(0,1))
                api.model.add_point(x,y,0)
            end
            for (a,b) in ((1,2),(2,3),(3,4),(4,1))
                api.model.add_line(a,b)
            end
            api.model.add_curve_loop([1,2,3,4]);api.model.add_plane_surface([1])
            for curve in 1:4
                api.mesh.set_transfinite_curve(curve,3)
            end
            api.mesh.set_transfinite_surface(1);api.mesh.set_recombine(2,1)
            api.mesh.generate(2)
            refined=api.mesh.refine()
            @test nnodes(refined)==25
            @test _api_refine_counts(refined)==Dict(3=>16)
            @test length(api.mesh.get_nodes(2,1)[1])==9
            @test all(curve->length(api.mesh.get_nodes(1,curve)[1])==3,1:4)
            @test length(api.mesh.get_nodes(2,1,true)[1])==25
            @test validate(refined).ok
            @test all(==(Int32(1)),api.LAST_MESH_CLASS[].cell_entities[3])
        finally
            api.finalize()
        end
    end
    @testset "support sharing across mixed volume families" begin
        cube=[0.0 1 1 0 0 1 1 0;0 0 1 1 0 0 1 1;0 0 0 0 1 1 1 1]
        hexa_prism=MixedMesh(hcat(cube,[2.0,0,0],[2.0,0,1]),[
            ElementBlock(5,reshape(Int32.(1:8),8,1)),
            ElementBlock(6,reshape(Int32[2,9,3,6,10,7],6,1))])
        pyrcoords=lagrange_nodes(7)
        pyramid_tet=MixedMesh(hcat(pyrcoords,[0.0,-2,0]),[
            ElementBlock(7,reshape(Int32.(1:5),5,1)),
            ElementBlock(4,reshape(Int32[1,2,5,6],4,1))])
        for (mesh,nnode,counts,volume) in (
                (hexa_prism,36,Dict(5=>8,6=>8),1.5),
                (pyramid_tet,18,Dict(4=>16,7=>4),5/3))
            try
                _api_refine_install(mesh)
                refined=api.mesh.refine()
                @test nnodes(refined)==nnode
                @test _api_refine_counts(refined)==counts
                @test _api_refine_volume(refined)≈volume
                @test validate(refined).ok
                faces=Dict{Tuple,Int}()
                for b in refined.blocks,column in axes(b.nodes,2),face in Tessella.Elements._VOLUME_CELL_FACES[Int(b.msh)]
                    key=Tuple(sort!([b.nodes[k,column] for k in face]))
                    faces[key]=get(faces,key,0)+1
                end
                @test all(count->count in (1,2),values(faces))
            finally
                api.finalize()
            end
        end
    end
    @testset "order-two support creation prunes unused primary nodes" begin
        try
            coords=hcat(lagrange_nodes(3),zeros(3))
            _api_refine_single_family(3;coords=coords)
            api.mesh.set_order(2)
            @test nnodes(api.mesh.get())==9
            refined=api.mesh.refine(max_nodes=9,max_cells=4)
            @test nnodes(refined)==9
            @test _api_refine_counts(refined)==Dict(3=>4)
            @test count(i->all(iszero,refined.coords[:,i]),axes(refined.coords,2))==1
            used=Set(vec(refined.blocks[1].nodes))
            @test length(used)==9
            @test length(api.LAST_MESH_CLASS[].node_entities)==9
            @test all(==((2,Int32(7))),api.LAST_MESH_CLASS[].node_entities)
        finally
            api.finalize()
        end
    end
    @testset "refinement normalizes explicitly reversed volume winding" begin
        for msh in (4,5,6,7)
            try
                _api_refine_single_family(msh)
                api.mesh.reverse()
                before=_api_refine_volume(api.mesh.get())
                @test before<0
                refined=api.mesh.refine()
                @test _api_refine_volume(refined)≈-before
                @test all(b->all(i->Tessella.Model._extrude_cell_signed_volume(
                    refined.coords,b.nodes,i,b.msh)>0,axes(b.nodes,2)),refined.blocks)
            finally
                api.finalize()
            end
        end
    end
    @testset "discrete P2 refinement resets interpolation supports" begin
        # RefineMesh first calls SetOrder1 even for complete quadratic input;
        # discrete entities then recreate support nodes from primary corners.
        for (msh,target) in ((10,3),(11,4))
            try
                expected=lagrange_nodes(msh)
                coords=copy(expected)
                for node in 5:size(coords,2)
                    coords[3,node]+=0.02*node
                end
                input=MixedMesh(coords,[ElementBlock(msh,
                    reshape(Int32.(1:size(coords,2)),:,1))])
                _api_refine_install(input)
                before=copy(input.coords)
                refined=api.mesh.refine()
                @test refined.coords≈expected
                @test input.coords==before
                @test _api_refine_counts(refined)==Dict(target=>(target==3 ? 4 : 8))
                @test validate(refined).ok
                target==4 && @test _api_refine_volume(refined)≈1/6
            finally
                api.finalize()
            end
        end
    end
    @testset "curved CAD placement fails before cache mutation" begin
        source="""
        Point(1)={1,0,0,1};Point(2)={2,0,0,1};Point(3)={0,2,0,1};
        Point(4)={0,1,0,1};Point(5)={0,0,0,1};
        Line(1)={1,2};Circle(2)={2,5,3};Line(3)={3,4};Circle(4)={4,5,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Transfinite Curve{:}=2;Transfinite Surface{1};Recombine Surface{1};
        """
        try
            api.initialize()
            mktempdir() do directory
                path=joinpath(directory,"curved_refine.geo");write(path,source)
                api.open_geo!(path)
            end
            api.mesh.generate(2)
            before=api.LAST_MESH[];class=api.LAST_MESH_CLASS[]
            @test_throws r"quadratic CAD placement.*Curve" api.mesh.refine()
            @test api.LAST_MESH[]===before
            @test api.LAST_MESH_CLASS[]===class
            @test nnodes(before)==4
        finally
            api.finalize()
        end
    end
    @testset "unrelated curved entities do not restrict a straight mesh" begin
        source="""
        Point(1)={0,0,0,1};Point(2)={1,0,0,1};
        Point(3)={1,1,0,1};Point(4)={0,1,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Transfinite Curve{1:4}=2;Transfinite Surface{1};Recombine Surface{1};
        Point(5)={3,0,0,1};Point(6)={2,0,0,1};Point(7)={2,1,0,1};
        Circle(5)={5,6,7};
        """
        try
            api.initialize()
            mktempdir() do directory
                path=joinpath(directory,"unrelated_curve.geo");write(path,source)
                api.open_geo!(path)
            end
            api.mesh.generate(2)
            refined=api.mesh.refine()
            @test nnodes(refined)==9
            @test _api_refine_counts(refined)==Dict(3=>4)
            @test validate(refined).ok
        finally
            api.finalize()
        end
    end
    @testset "planar ruled extrusion faces admit linear refinement" begin
        source="""
        Point(1)={0,0,0,1};Point(2)={1,0,0,1};
        Point(3)={1,1,0,1};Point(4)={0,1,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Transfinite Curve{1:4}=2;Transfinite Surface{1};Recombine Surface{1};
        Extrude{0,0,1}{Surface{1};Layers{2};Recombine;}
        """
        try
            api.initialize()
            mktempdir() do directory
                path=joinpath(directory,"planar_ruled_refine.geo");write(path,source)
                api.open_geo!(path)
            end
            api.mesh.generate(3)
            refined=api.mesh.refine(max_nodes=45,max_cells=16)
            @test nnodes(refined)==45
            @test _api_refine_counts(refined)==Dict(5=>16)
            @test _api_refine_volume(refined)≈1
            @test validate(refined).ok
        finally
            api.finalize()
        end
    end
end
