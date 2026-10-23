# `fltk`/`graphics` modules UI01: deterministic headless semantics — the
# capability query reports no backend, event-loop and window calls are
# validated no-ops, and selection calls return empty selections.
module FltkGraphicsChecks
using Test, Tessella
const A=Tessella.API

@testset "capability reporting is honest" begin
    # `is_available` needs no session — it is a build-capability query.
    @test A.fltk.is_available()==0
end

@testset "headless no-ops" begin
    A.initialize()
    try
        # Event-loop entry points return immediately without a display.
        @test A.fltk.run()===nothing
        @test A.fltk.run("options.opt")===nothing
        @test A.fltk.wait()===nothing
        @test A.fltk.wait(0.5)===nothing
        @test A.fltk.update()===nothing
        @test A.fltk.awake()===nothing
        @test A.fltk.awake("e")===nothing
        @test A.fltk.lock()===nothing
        @test A.fltk.unlock()===nothing
        @test A.fltk.initialize()===nothing
        @test A.fltk.finalize()===nothing
        # `graphics.draw` is a no-op like upstream with no graphics context.
        @test A.graphics.draw()===nothing
        # Selection returns empty sets — nothing is displayed.
        @test A.fltk.select_entities()==Tuple{Int,Int32}[]
        @test A.fltk.select_entities(2)==Tuple{Int,Int32}[]
        @test A.fltk.select_elements()==Int[]
        @test A.fltk.select_views()==Int[]
        # Window/status helpers validate inputs then no-op.
        @test A.fltk.split_current_window()===nothing
        @test A.fltk.split_current_window("h",0.25)===nothing
        @test A.fltk.set_current_window(0)===nothing
        @test A.fltk.set_status_message("hi")===nothing
        @test A.fltk.set_status_message("hi",true)===nothing
        @test A.fltk.show_context_window(2,5)===nothing
        @test A.fltk.open_tree_item("Modules")===nothing
        @test A.fltk.close_tree_item("Modules")===nothing
        # Bad arguments still reject.
        @test_throws ArgumentError A.fltk.wait(Inf)
        @test_throws ArgumentError A.fltk.wait("soon")
        @test_throws ArgumentError A.fltk.split_current_window("x")
        @test_throws ArgumentError A.fltk.split_current_window("v",1.5)
        @test_throws ArgumentError A.fltk.set_current_window(-1)
        @test_throws ArgumentError A.fltk.set_status_message(3)
        @test_throws ArgumentError A.fltk.show_context_window(4,1)
        @test_throws ArgumentError A.fltk.show_context_window(1,0)
        @test_throws ArgumentError A.fltk.run(3)
        @test_throws ArgumentError A.fltk.awake(3)
    finally
        A.finalize()
    end
end

@testset "fltk calls require an initialized session" begin
    A.initialize(); A.finalize()
    @test_throws Exception A.fltk.run()
    @test_throws Exception A.graphics.draw()
    @test_throws Exception A.fltk.select_elements()
end

end # module
