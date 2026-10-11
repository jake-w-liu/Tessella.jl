module OrientationCarrierTraceRegression
using Test,Tessella
const Model=Tessella.Model
const E=Tessella.Elements

function actual_mesh(coords,msh,cells,tags)
    bbox=ntuple(index->index<=3 ? minimum(@view coords[index,:]) :
        maximum(@view coords[index-3,:]),6)
    data=E.MixedEntityData(Dict((2,21)=>E.MixedEntity(2,21,bbox));
        node_entities=fill((2,Int32(21)),size(coords,2)),
        node_parametric=fill(nothing,size(coords,2)),
        external_node_tags=UInt64.(101:100+size(coords,2)),
        block_entities=[fill(Int32(21),size(cells,2))],external_element_tags=[UInt64.(tags)])
    return Tessella.MixedMesh(coords,[E.ElementBlock(msh,cells)];entity_data=data)
end

function orient(actual,carrier,tags,changes)
    replacement=Model._model_orientation_mixed_replacement(actual,changes)
    return Model._model_orientation_carrier(actual,replacement,carrier,
        collect(Int32,1:size(carrier.coords,2)),changes;cell_tags=tags)
end

@testset "Orientation carrier follows emitted public cell identity" begin
    xyz=[0. 1 0;0 0 1;0 0 0]
    triangles=Int32[1 1;2 2;3 3]
    actual=actual_mesh(xyz,2,triangles,[801,802])
    carrier=Tessella.Mesh(xyz;tris=triangles,segs=reshape(Int32[1,2],2,1))
    tags=[UInt64[0],UInt64[801,802],UInt64[]]
    changed=Dict(1=>Int32[1,3,2])
    oriented=orient(actual,carrier,tags,changed)
    @test oriented.tris==Int32[1 1;3 2;2 3]
    @test oriented.segs==carrier.segs
    @test oriented.coords==carrier.coords
    @test carrier.tris==triangles
    @test orient(actual,carrier,tags,Dict(1=>Int32[1,2,3]))===carrier
    @test_throws ArgumentError orient(actual,carrier,
        [UInt64[0],UInt64[0,802],UInt64[]],changed)
    @test carrier.tris==triangles
end

@testset "Repeated parent trace reverses every triangulated Quad alias" begin
    xyz=[0. 1 1 0;0 0 1 1;0 0 0 0]
    actual=actual_mesh(xyz,3,reshape(Int32[1,2,3,4],4,1),[901])
    carrier=Tessella.Mesh(xyz;tris=Int32[1 1;2 3;3 4])
    oriented=orient(actual,carrier,[UInt64[],UInt64[901,901],UInt64[]],
        Dict(1=>Int32[1,4,3,2]))
    @test oriented.tris==Int32[1 1;3 4;2 3]
    @test oriented.coords==carrier.coords
    @test carrier.tris==Int32[1 1;2 3;3 4]
end

@testset "Mixed P1 geometry keeps its actual P2 Quad parent" begin
    xyz=[0. 1 1 0 .5 1 .5 0 .5;0 0 1 1 0 .5 1 .5 .5;zeros(1,9)]
    actual=actual_mesh(xyz,10,reshape(Int32.(1:9),9,1),[1001])
    carrier=Tessella.MixedMesh(xyz,[E.ElementBlock(3,reshape(Int32[1,2,3,4],4,1))];
        physical_names=Dict((2,7)=>"geometry face"))
    oriented=orient(actual,carrier,[UInt64[1001]],Dict(1=>Int32[1,4,3,2,8,7,6,5,9]))
    @test oriented.blocks[1].nodes==reshape(Int32[1,4,3,2],4,1)
    @test oriented.coords==carrier.coords
    @test oriented.physical_names==carrier.physical_names
    @test carrier.blocks[1].nodes==reshape(Int32[1,2,3,4],4,1)
end

@testset "Unrelated unbound Surface alias survives orientation" begin
    xyz=[0. 1 0 3 4 3;0 0 1 0 0 1;zeros(1,6)]
    cells=Int32[1 4;2 5;3 6]
    actual=actual_mesh(xyz,2,cells,[1101,1102])
    carrier=Tessella.Mesh(xyz;tris=cells)
    oriented=orient(actual,carrier,[UInt64[],UInt64[1101,0],UInt64[]],Dict(1=>Int32[1,3,2]))
    @test oriented.tris==Int32[1 4;3 5;2 6]
    @test oriented.coords==carrier.coords
    @test carrier.tris==cells
end
end
