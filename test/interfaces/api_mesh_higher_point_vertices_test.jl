module HigherPointVerticesTests
using Test,Tessella
const API=Tessella.API
# Exact native generate2 rows independently recorded by Gmsh 4.15.2.
# point_primary_v4.json SHA256
# 75161be150ba8924b0a1910568adbba199309f8fcf02366420bee9d10321b141.
const ROWS=[
    (order=1,renumber=0,supplied=false,points=[101,102,103],cells=[1],point_conn=[103],lines=[101,107,107,104]),
    (order=1,renumber=0,supplied=true,points=[101,102,103],cells=[777,778],point_conn=[101,102],lines=[101,107,107,104]),
    (order=1,renumber=1,supplied=false,points=[1,2,3],cells=[1],point_conn=[3],lines=[1,7,7,4]),
    (order=1,renumber=1,supplied=true,points=[1,2,3],cells=[1,2],point_conn=[1,2],lines=[1,7,7,4]),
    (order=2,renumber=0,supplied=false,points=[103],cells=[1],point_conn=[103],lines=[101,107,112,107,104,113]),
    (order=2,renumber=0,supplied=true,points=[101,102],cells=[777,778],point_conn=[101,102],lines=[101,107,112,107,104,113]),
    (order=2,renumber=1,supplied=false,points=[1],cells=[1],point_conn=[1],lines=[5,6,7,6,2,8]),
    (order=2,renumber=1,supplied=true,points=[1,2],cells=[1,2],point_conn=[1,2],lines=[1,6,7,6,3,8]),
]

function setup(;order=1,renumber=0,only_empty=1)
    API.initialize();API.option("Mesh.Renumber",renumber)
    API.option("Mesh.ElementOrder",order);API.option("Mesh.MeshOnlyEmpty",only_empty)
    for (tag,(x,y)) in enumerate(((0.,0.),(1.,0.),(1.,1.),(0.,1.)))
        API.model.add_point(x,y,0.;tag,meshSize=1.)
    end
    for (tag,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        API.model.add_line(a,b;tag);API.mesh.set_transfinite_curve(tag,3)
    end
    API.model.add_curve_loop([1,2,3,4];tag=1)
    API.model.add_plane_surface([1];tag=1);API.mesh.set_transfinite_surface(1)
end

@testset "Actual higher generation preserves native Point vertex identities" begin
    for only_empty in (0,1),row in ROWS
        @testset "E$only_empty-P$(row.order)-R$(row.renumber)-cells$(row.supplied)" begin
            try
                setup(;order=row.order,renumber=row.renumber,only_empty)
                API.mesh.add_nodes(0,1,[101,102,103],zeros(9))
                row.supplied && API.mesh.add_elements_by_type(1,15,[777,778],[101,102])
                API.mesh.generate(2)
                tags,coords,params=API.mesh.get_nodes(0,1)
                @test tags==UInt64.(row.points) && coords==zeros(3length(row.points)) && isempty(params)
                cells,connections=API.mesh.get_elements_by_type(15,1)
                @test cells==UInt64.(row.cells) && connections==UInt64.(row.point_conn)
                types,_,line_nodes=API.mesh.get_elements(1,1)
                @test types==Int32[row.order==1 ? 1 : 8] && only(line_nodes)==UInt64.(row.lines)
                all_nodes=API.mesh.get_nodes()[1];all_cells=API.mesh.get_elements()
                @test allunique(all_nodes) && allunique(vcat(all_cells[2]...))
                @test Set(vcat(all_cells[3]...)) ⊆ Set(all_nodes)
                surface=API.mesh.get_elements(2,1)[2]
                @test sum(API.mesh.get_element_qualities(vcat(surface...),"volume"))≈1 atol=2e-11
                if row.order==1
                    product=Tessella.Model._model_mesh_product(API.CURRENT[].meshing.mesh_authority,0,1)
                    @test Tessella.MeshTypes.nnodes(product)==3
                end
            finally
                API.finalize()
            end
        end
    end
end

@testset "Actual higher products preserve distinct Point cells and coincident owners" begin
    try
        setup()
        API.mesh.add_nodes(0,1,[101,102,103],zeros(9))
        API.mesh.add_elements_by_type(1,15,[777,778],[101,101])
        API.model.add_point(0.,0.,0.;tag=5)
        API.mesh.add_nodes(0,5,[501],zeros(3))
        API.mesh.add_elements_by_type(5,15,[503],[501])
        API.mesh.generate(2)
        @test API.mesh.get_elements_by_type(15,1)==(UInt64[777,778],UInt64[101,101])
        @test API.mesh.get_elements_by_type(15,5)==(UInt64[503],UInt64[501])
        @test API.mesh.get_nodes(0,1)[1]==UInt64[101,102,103]
        @test API.mesh.get_nodes(0,5)[1]==UInt64[501]
        @test API.mesh.get_node(101)[1]==API.mesh.get_node(501)[1]==zeros(3)
        @test first(only(API.mesh.get_elements(1,1)[3]))==101
    finally
        API.finalize()
    end
end
end
