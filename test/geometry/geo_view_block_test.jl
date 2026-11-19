# Parsed `View "name" { records };` blocks inside `.geo` scripts and
# `Include`d `.pos` files (tutorials t4/t8/t9, examples/post_processing):
# names and records may use script variables and string expressions,
# `PostProcessing.NbViews` and `View[i].Max` follow the view list, and a
# color table may be set with or without views. The view data below is what
# Gmsh 4.15.2 holds for the same script (`gmsh.view.getListData`).
module GeoViewBlockTest
using Test, Tessella

@testset "geo View blocks" begin
    dir=mktempdir()
    write(joinpath(dir,"extra.pos"),"""
        View "included" {
          SL(0,0,0, 1,0,0){-4,7};
        };
        """)
    path=joinpath(dir,"views.geo")
    write(path,"""
        View.ColorTable = {Red, Green, Blue};
        a = 2;
        n0 = PostProcessing.NbViews;
        View "first" {
          ST(0,0,0, 1,0,0, 0,1,0){1,2,a+1};
          SP(a,0,0){5};
          T2(10, -10, 0){ StrCat("v", "=", Sprintf("%g", a)) };
        };
        n1 = PostProcessing.NbViews;
        For i In {1:2}
          View Sprintf("loop %g", i) {
            SQ(i,0,0, i+1,0,0, i+1,1,0, i,1,0){i,i,i,i};
          };
          View[PostProcessing.NbViews-1].ColorTable = {{i*100,0,0}};
        EndFor
        Include "extra.pos";
        n4 = PostProcessing.NbViews;
        m0 = View[0].Max;
        m3 = View[3].Min;
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test (ex.values["n0"],ex.values["n1"],ex.values["n4"])==(0.0,1.0,4.0)
    @test (ex.values["m0"],ex.values["m3"])==(5.0,-4.0)
    names=[view.name for view in ex.views]
    @test names==["first","loop 1","loop 2","included"]
    first=ex.views[1]
    @test [element.kind for element in first.elements]==[:triangle,:point]
    @test vec(first.elements[1].values)==[1.0,2.0,3.0]
    @test first.elements[2].coords[:,1]==[2.0,0.0,0.0]
    @test only(first.texts).strings==["v=2"]
    @test vec(ex.views[3].elements[1].values)==fill(2.0,4)
end
end
