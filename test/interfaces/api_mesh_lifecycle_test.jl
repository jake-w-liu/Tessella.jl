module APIMeshLifecycleTest
using Test
using Tessella
using Tessella.MeshTypes: validate, nnodes, ntets

mesh_crc(mesh)=mesh isa Tessella.MeshTypes.Mesh ? Tessella.MeshTypes.mesh_crc(mesh) : mixed_crc(mesh)

const _MESH_LIFECYCLE_API=Tessella.API

@testset "atomic cached-mesh refinement and clearing through API" begin
    _MESH_LIFECYCLE_API.finalize()
    @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine()
    @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.clear()

    try
        _MESH_LIFECYCLE_API.initialize()
        @test _MESH_LIFECYCLE_API.mesh.clear()===nothing
        @test _MESH_LIFECYCLE_API.mesh.clear([])===nothing
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine()

        @test _MESH_LIFECYCLE_API.model.add_box(
            0,0,0,1,1,1;tag=1)==1
        generated=_MESH_LIFECYCLE_API.mesh.generate(3)
        generated_crc=mesh_crc(generated)
        counts()=Dict(Int(t)=>length(e) for (t,e) in
            zip(_MESH_LIFECYCLE_API.mesh.get_elements()...))
        base=counts()
        nodes=length(_MESH_LIFECYCLE_API.mesh.get_nodes()[1])
        edges=Set{Tuple{Int,Int}}()
        for (type,_,cells) in zip(_MESH_LIFECYCLE_API.mesh.get_elements()...)
            type==4 || continue
            for k in 0:4:length(cells)-4,(a,b) in ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
                push!(edges,minmax(Int(cells[k+a]),Int(cells[k+b])))
            end
        end

        # Gmsh 4.15.2 refine(): every edge gets a midpoint and each cell splits
        # into 2, 4 or 8 (339 nodes / 1733 edges / 1125 tetrahedra become
        # 2072 nodes / 9000 tetrahedra there).
        refined=_MESH_LIFECYCLE_API.mesh.refine()
        refined_crc=mesh_crc(refined)
        @test validate(refined).ok
        after=counts()
        @test length(_MESH_LIFECYCLE_API.mesh.get_nodes()[1])==nodes+length(edges)
        @test after[4]==8base[4] && after[2]==4base[2] && after[1]==2base[1]

        stored=_MESH_LIFECYCLE_API.LAST_MESH[]
        @test stored!==nothing && stored!==refined
        @test !Base.mightalias(refined.coords,stored.coords)
        refined.coords[1,1]+=17.0
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==refined_crc

        # Resource limits fail before anything changes.
        total=nodes+length(edges)
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine(max_nodes=total)
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==refined_crc
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine(max_cells=8after[4])
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==refined_crc
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine(max_nodes=true)
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==refined_crc

        # Clearing the Volume keeps the refined boundary mesh; generate(3)
        # fills the Volume again from it, as in Gmsh (1082 nodes and 2160
        # triangles kept, then 1800 nodes / 7573 tetrahedra there).
        @test _MESH_LIFECYCLE_API.mesh.clear([(3,1)])===nothing
        leftover=counts()
        @test !haskey(leftover,4) && leftover[2]==after[2]
        @test isempty(_MESH_LIFECYCLE_API.mesh.get_nodes(3,1)[1])
        _MESH_LIFECYCLE_API.mesh.generate(3)
        again=counts()
        @test again[2]==after[2] && again[4]>after[4]÷2
        @test validate(_MESH_LIFECYCLE_API.mesh.get()).ok
        twice_crc=mesh_crc(_MESH_LIFECYCLE_API.mesh.get())
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.clear([(3,99)])
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==twice_crc
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.clear(1)
        @test mesh_crc(_MESH_LIFECYCLE_API.mesh.get())==twice_crc

        @test _MESH_LIFECYCLE_API.mesh.clear(())===nothing
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.get()
        @test_throws ArgumentError _MESH_LIFECYCLE_API.mesh.refine()
        @test _MESH_LIFECYCLE_API.model.get_entities()==vcat(
            [(0,i) for i in 1:8],[(1,i) for i in 1:12],
            [(2,i) for i in 1:6],[(3,1)])

        regenerated=_MESH_LIFECYCLE_API.mesh.generate(3)
        @test mesh_crc(regenerated)==generated_crc
        @test _MESH_LIFECYCLE_API.mesh.clear(Int32[])===nothing
        @test _MESH_LIFECYCLE_API.model.get_entities()==vcat(
            [(0,i) for i in 1:8],[(1,i) for i in 1:12],
            [(2,i) for i in 1:6],[(3,1)])
    finally
        _MESH_LIFECYCLE_API.finalize()
    end

end

@testset "periodic maps follow cached refinement" begin
    _MESH_LIFECYCLE_API.finalize()
    try
        _MESH_LIFECYCLE_API.initialize()
        for (tag,(x,y)) in enumerate(((0.0,0.0),(1.0,0.0),
                                      (1.0,1.0),(0.0,1.0)))
            @test _MESH_LIFECYCLE_API.model.add_point(
                x,y,0;tag=tag,meshSize=0.5)==tag
        end
        for (tag,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
            @test _MESH_LIFECYCLE_API.model.add_line(
                first,last;tag=tag)==tag
        end
        @test _MESH_LIFECYCLE_API.model.add_curve_loop(
            [1,2,3,4];tag=1)==1
        @test _MESH_LIFECYCLE_API.model.add_plane_surface([1];tag=1)==1
        translate_x=(1.0,0.0,0.0,1.0,
                     0.0,1.0,0.0,0.0,
                     0.0,0.0,1.0,0.0,
                     0.0,0.0,0.0,1.0)
        @test _MESH_LIFECYCLE_API.mesh.set_periodic(
            1,[2],[4],translate_x)===nothing

        generated=_MESH_LIFECYCLE_API.mesh.generate(2)
        @test validate(generated).ok
        @test length(_MESH_LIFECYCLE_API.mesh.get_periodic_nodes(
            1,2).slave_nodes)==3

        refined=_MESH_LIFECYCLE_API.mesh.refine()
        @test validate(refined).ok
        mapping=_MESH_LIFECYCLE_API.mesh.get_periodic_nodes(1,2)
        @test length(mapping.slave_nodes)==length(mapping.master_nodes)==5
        # Node tags, as in Gmsh. The refined midpoints of the two curves are
        # their own 1-D nodes' midpoints, so y differs by rounding, with the
        # same values in Gmsh 4.15.2 (slave y, master y):
        gmsh_pairs=Set([(1.0,1.0),(0.0,0.0),(0.25000000000007233,0.2500000000010296),
            (0.5,0.500000000002059),(0.75,0.7500000000008381)])
        pairs=Set{Tuple{Float64,Float64}}()
        for (slave,master) in zip(mapping.slave_nodes,mapping.master_nodes)
            s=_MESH_LIFECYCLE_API.mesh.get_node(slave)[1]
            m=_MESH_LIFECYCLE_API.mesh.get_node(master)[1]
            @test (s[1],s[3])==(m[1]+1,m[3])
            push!(pairs,(s[2],m[2]))
        end
        if Sys.ARCH===:x86_64
            @test pairs==gmsh_pairs
        else
            @test all(p->any(q->all(isapprox.(p,q;atol=1e-12)),gmsh_pairs),pairs)
        end

        # Gmsh 4.15.2 keeps the periodic link and reports no nodes after clear().
        @test _MESH_LIFECYCLE_API.mesh.clear()===nothing
        cleared=_MESH_LIFECYCLE_API.mesh.get_periodic_nodes(1,2)
        @test isempty(cleared.slave_nodes) && isempty(cleared.master_nodes)
    finally
        _MESH_LIFECYCLE_API.finalize()
    end
end

end
