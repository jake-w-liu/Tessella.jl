module APIViewTagCounterTest

using Test, Tessella

const API = Tessella.API

# Gmsh 4.15.2 numbers views with the process-wide PView::_globalTag: it starts
# at 1 in a fresh process, gmsh.finalize and gmsh.clear (DeleteAllModelsAndViews)
# reset it to 0, and gmsh.initialize leaves it alone. So a script that runs
# initialize/finalize twice gets views 1, 2, 3 the first time and 0, 1, 2 the
# second; after gmsh.clear the next view is 0 again. (Checked with the Python
# wheel; the fresh-process 1, 2, 3 is not observable inside the test suite.)
@testset "view tags follow Gmsh's global view counter" begin
    API.initialize()
    API.option("General.Terminal", 0)
    API.finalize()
    API.initialize()
    try
        API.option("General.Terminal", 0)
        @test [API.view.add("a"), API.view.add("b"), API.view.add("c")] == [0, 1, 2]
        API.clear()
        @test API.view.add("d") == 0
        @test API.view.add("e") == 1
    finally
        API.finalize()
    end
    API.initialize()
    try
        API.option("General.Terminal", 0)
        @test API.view.add("f") == 0
    finally
        API.finalize()
    end
end

end
