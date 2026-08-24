using Test
using Tessella
using Tessella.GUI: GuiState, camera_reset!, gui_command!, gui_new!

@testset "validated headless GUI state machine" begin
    selected=[1,2];log=["existing"]
    owned=GuiState(false,(1,2.0,3),selected,log,4,5)
    selected[1]=9;log[1]="changed"
    @test owned.camera==(1.0,2.0,3.0)
    @test owned.selected==[1,2] && owned.log==["existing"]
    @test owned.selected!==selected && owned.log!==log
    @test_throws ArgumentError GuiState(false,(true,0,1),Int[],String[],0,0)
    @test_throws ArgumentError GuiState(false,(Inf,0,1),Int[],String[],0,0)
    @test_throws ArgumentError GuiState(false,(0,0,1),[1,1],String[],0,0)
    @test_throws ArgumentError GuiState(false,(0,0,1),[0],String[],0,0)
    @test_throws ArgumentError GuiState(false,(0,0,1),[true],String[],0,0)
    @test_throws ArgumentError GuiState(false,(0,0,1),Int[],String[],-1,0)
    @test_throws ArgumentError GuiState(false,(0,0,1),Int[],String[],true,0)
    @test_throws ArgumentError GuiState(false,(0,0,1),Int[],String[],big(typemax(Int))+1,0)

    s=GuiState()
    @test !s.initialized && s.camera==(0.0,0.0,1.0)
    @test isempty(s.selected) && isempty(s.log)
    @test_throws ArgumentError gui_command!(s,"")
    @test_throws ArgumentError gui_command!(s,"box 0 0 0 1 1 1 1")
    @test_throws ArgumentError gui_command!(s,"new extra")
    @test_throws ArgumentError gui_command!(s,repeat("x",1_000_001))
    @test_throws ArgumentError gui_command!(s,join(fill("x",10_001),' '))

    try
        @test gui_new!(s)===s
        @test s.initialized && s.log==["new model"]
        @test_throws ArgumentError gui_command!(s,"box 0 0 0 1 1 1")
        @test_throws ArgumentError gui_command!(s,"box 0 0 0 Inf 1 1 1")
        @test_throws ArgumentError gui_command!(s,"box 0 0 0 -1 1 1 1")
        @test s.log==["new model"]

        @test gui_command!(s,"box 0 0 0 1 1 1 1")===s
        @test s.log[end]=="box 1"
        @test_throws ArgumentError gui_command!(s,"mesh 4")
        @test s.mesh_nodes==0 && s.mesh_tets==0
        @test gui_command!(s,"mesh 3")===s
        @test s.mesh_nodes>0 && s.mesh_tets>0

        previous_log=copy(s.log)
        @test_throws ArgumentError gui_command!(s,"select")
        @test_throws ArgumentError gui_command!(s,"select 1 1")
        @test_throws ArgumentError gui_command!(s,"select 0")
        @test_throws ArgumentError gui_command!(s,"select nope")
        @test s.selected==Int[] && s.log==previous_log
        gui_command!(s,"select 1 2")
        @test s.selected==[1,2]

        prior_camera=s.camera;prior_log=copy(s.log)
        @test_throws ArgumentError gui_command!(s,"camera 1 2")
        @test_throws ArgumentError gui_command!(s,"camera 1 NaN 3")
        @test_throws ArgumentError gui_command!(s,"camera 1 nope 3")
        @test s.camera==prior_camera && s.log==prior_log
        gui_command!(s,"camera 1 2 3")
        @test s.camera==(1.0,2.0,3.0)
        @test camera_reset!(s)==(0.0,0.0,1.0)
        gui_command!(s,"camera 1 2 3")
        @test_throws ArgumentError gui_command!(s,"reset extra")
        gui_command!(s,"reset")
        @test s.camera==(0.0,0.0,1.0)
        @test_throws ArgumentError gui_command!(s,"explode")

        s.selected=[1,1]
        @test_throws ArgumentError gui_command!(s,"reset")
        s.selected=Int[]
        s.camera=(NaN,0.0,1.0)
        @test_throws ArgumentError gui_command!(s,"reset")
        @test camera_reset!(s)==(0.0,0.0,1.0)

        previous_entries=length(s.log)
        gui_command!(s,"new")
        @test s.mesh_nodes==0 && s.mesh_tets==0 && isempty(s.selected)
        @test length(s.log)==previous_entries+1 && s.log[end]=="new model"
    finally
        Tessella.API.finalize()
    end

    @test isempty(Docs.undocumented_names(Tessella.GUI;private=false))
end
