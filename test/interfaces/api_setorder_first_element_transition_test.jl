module SetOrderFirstElementTransitionChecks
# HighOrder.cpp transitions each entity from the order of its first stored
# element: an entity already at the requested order keeps every cell unchanged
# (including lower-order ones), while all other entities convert their cells.
# The mixed-order guard in `_dim01_quadratic!` used to reject every mixed-order
# entity, so discrete or generated surfaces mixing Tri3 and Tri6 cells always
# threw `first-element high-order transition is not implemented`.
using Test, Tessella
const A=Tessella.API

# A discrete surface whose two triangles share edge (2,3): the Tri6 stores its
# midpoint node 6 on that edge, so a lifted Tri3 must reuse node 6 rather than
# create a duplicate support.
const _SURFACE_NODES=Float64[0,0,0, 1,0,0, 0,1,0, 1,1,0, .5,0,0, .5,.5,0, 0,.5,0]
const _SURFACE_TRI6=Int32[1,2,3,5,6,7]
const _SURFACE_TRI3=Int32[2,4,3]

function discrete_surface!(first)
    A.model.add_discrete_entity(2,7)
    A.mesh.add_nodes(2,7,collect(1:7),_SURFACE_NODES)
    if first===:linear
        A.mesh.add_elements_by_type(7,2,[101],Vector{Int64}(_SURFACE_TRI3))
        A.mesh.add_elements_by_type(7,9,[102],Vector{Int64}(_SURFACE_TRI6))
    else
        A.mesh.add_elements_by_type(7,9,[102],Vector{Int64}(_SURFACE_TRI6))
        A.mesh.add_elements_by_type(7,2,[101],Vector{Int64}(_SURFACE_TRI3))
    end
    A.mesh.generate(0)
    return nothing
end

function surface_cells()
    types,tags,connectivity=A.mesh.get_elements(2,7)
    result=Dict{Int32,Vector{UInt64}}()
    for (type,cell_tags,flat) in zip(types,tags,connectivity)
        width=Int(type)==9 ? 6 : Int(type)==2 ? 3 : error("unexpected type $type")
        for (index,tag) in enumerate(cell_tags)
            result[Int32(tag)]=flat[width*(index-1)+1:width*index]
        end
    end
    return types,result
end

function node_point(tag)
    xyz,_,_,_=A.mesh.get_node(tag)
    return collect(Float64,xyz)
end

@testset "linear-first discrete surface converts every cell" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        discrete_surface!(:linear)
        A.mesh.set_order(2)
        types,cells=surface_cells()
        # The linear cell is first: the entity converts to order 2, keeping the
        # existing Tri6's nodes and sharing its midpoint on edge (2,3).
        @test types==Int32[9]
        @test any(cell->cell==_SURFACE_TRI6,values(cells))
        lifted=only(cell for cell in values(cells) if cell[1:3]==_SURFACE_TRI3)
        @test lifted[6]==Int32(6)
        @test node_point(lifted[4])≈[1.0,0.5,0.0]
        @test node_point(lifted[5])≈[0.5,1.0,0.0]
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "quadratic-first discrete surface is already at the target" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        discrete_surface!(:quadratic)
        A.mesh.set_order(2)
        types,cells=surface_cells()
        # The Tri6 is first: the entity already runs at order 2, so Gmsh leaves
        # every stored cell untouched — the Tri3 stays linear.
        @test types==Int32[2,9]
        @test any(cell->cell==_SURFACE_TRI3,values(cells))
        @test any(cell->cell==_SURFACE_TRI6,values(cells))
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "entities transition independently" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        discrete_surface!(:quadratic)
        A.model.add_discrete_entity(2,8)
        A.mesh.add_nodes(2,8,[10,11,12],Float64[2,0,0, 3,0,0, 2,1,0])
        A.mesh.add_elements_by_type(8,2,[103],[10,11,12])
        A.mesh.generate(0)
        A.mesh.set_order(2)
        types,_=surface_cells()
        @test types==Int32[2,9]
        types8,_,conn8=A.mesh.get_elements(2,8)
        @test types8==Int32[9]
        lifted8=reduce(vcat,conn8)
        # generate(0) compacts node tags; identify the lifted corners by point.
        @test node_point(lifted8[1])≈[2.0,0.0,0.0]
        @test node_point(lifted8[2])≈[3.0,0.0,0.0]
        @test node_point(lifted8[3])≈[2.0,1.0,0.0]
        @test node_point(lifted8[4])≈[2.5,0.0,0.0]
    finally
        A.finalize()
    end
end

@testset "linear-first discrete curve converts every cell" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add_discrete_entity(1,5)
        A.mesh.add_nodes(1,5,collect(1:4),
            Float64[0,0,0, 1,0,0, 2,0,0, 1.5,0,0])
        A.mesh.add_elements_by_type(5,1,[201],[1,2])
        A.mesh.add_elements_by_type(5,8,[202],[2,3,4])
        A.mesh.generate(0)
        A.mesh.set_order(2)
        types,_,connectivity=A.mesh.get_elements(1,5)
        @test types==Int32[8]
        flat=reduce(vcat,connectivity)
        @test flat[1:2]==Int32[1,2]
        @test flat[4:6]==Int32[2,3,4]
        @test node_point(flat[3])≈[0.5,0.0,0.0]
    finally
        A.finalize()
    end
end

@testset "quadratic-first discrete curve is already at the target" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add_discrete_entity(1,5)
        A.mesh.add_nodes(1,5,collect(1:4),
            Float64[0,0,0, 1,0,0, 2,0,0, 1.5,0,0])
        A.mesh.add_elements_by_type(5,8,[202],[2,3,4])
        A.mesh.add_elements_by_type(5,1,[201],[1,2])
        A.mesh.generate(0)
        A.mesh.set_order(2)
        # The Line3 is first: the curve already runs at order 2 and the Line2
        # stays stored. API.mesh.get_elements reports only the curve's first
        # type, so inspect the mixed blocks directly.
        blocks=A.mesh.get().blocks
        @test sort!(Int32[block.msh for block in blocks])==Int32[1,8]
        @test any(block->block.msh==1 && vec(block.nodes)==Int32[1,2],blocks)
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "recombination keeps high-order Surface cells untouched" begin
    # Native recombination pairs first-order Surface elements only: the P1
    # entity's triangles still recombine while the P2 entity stays stored.
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add_discrete_entity(2,7)
        A.mesh.add_nodes(2,7,collect(1:4),Float64[0,0,0, 1,0,0, 0,1,0, 1,1,0])
        A.mesh.add_elements_by_type(7,2,[101,102],[1,2,3, 2,4,3])
        A.model.add_discrete_entity(2,8)
        A.mesh.add_nodes(2,8,collect(5:10),
            Float64[2,0,0, 3,0,0, 2,1,0, 2.5,0,0, 2.5,0.5,0, 2,0.5,0])
        A.mesh.add_elements_by_type(8,9,[110],[5,6,7,8,9,10])
        A.mesh.generate(0)
        A.mesh.set_recombine(2,7)
        A.mesh.recombine()
        types7,_,_=A.mesh.get_elements(2,7)
        @test types7==Int32[3]
        types8,tags8,conn8=A.mesh.get_elements(2,8)
        @test types8==Int32[9]
        @test sum(length,tags8)==1
        flat8=reduce(vcat,conn8)
        @test length(flat8)==6
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "recombination on an all high-order Surface mesh is a no-op" begin
    # With every stored Surface cell above order one, nothing is pairable and
    # the recombine plan leaves the cache itself untouched.
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.Renumber",0)
        A.model.add_discrete_entity(2,8)
        A.mesh.add_nodes(2,8,collect(1:7),_SURFACE_NODES)
        A.mesh.add_elements_by_type(8,9,[101],Vector{Int64}(_SURFACE_TRI6))
        A.mesh.generate(0)
        before=A.LAST_MESH[]
        A.mesh.set_recombine(2,8)
        A.mesh.recombine()
        @test A.LAST_MESH[]===before
        types,tags,_=A.mesh.get_elements(2,8)
        @test types==Int32[9]
    finally
        A.finalize()
    end
end

@testset "generation takes the general order path above the target" begin
    # A discrete entity storing a P3 cell among P1 cells is not a P2 pass:
    # like `set_order`, `generate` at order 2 routes the whole entity through
    # the general order factory, which converts the P3 cell from its corner
    # skeleton while lifting the P1 cell — no first-element throw.
    A.initialize()
    try
        A.option("General.Terminal",0)
        discrete_surface!(:linear)
        A.mesh.add_nodes(2,7,collect(8:17),
            Float64[3,0,0, 4,0,0, 3,1,0,
                    10/3,0,0, 11/3,0,0,
                    11/3,1/3,0, 10/3,2/3,0,
                    3,2/3,0, 3,1/3,0,
                    10/3,1/3,0])
        A.mesh.add_elements_by_type(7,21,[110],collect(8:17))
        A.mesh.generate(0)
        types,_,_=A.mesh.get_elements(2,7)
        @test types==Int32[2,9,21]
        A.option("Mesh.ElementOrder",2)
        A.mesh.generate(2)
        types,_,connectivity=A.mesh.get_elements(2,7)
        @test types==Int32[9]
        for conn in connectivity
            @test length(conn)%6==0
        end
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

end
