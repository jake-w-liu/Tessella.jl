# Regression pin for the 22-case HFSS native-meshing set (STATUS #12): every HFSS
# User Guide case geometry (ch. 5–10) is meshed FROM SCRATCH with Tessella's own
# primitives / raw surfaces (no gmsh, no OpenCASCADE) to a valid, watertight,
# conforming tet mesh. The builders live in validation/hfss_cases/ (also a standalone
# runner); this asserts the whole set stays green.

include(joinpath(@__DIR__, "..", "validation", "hfss_cases", "hfss_case_meshes.jl"))

@testset "HFSS 22-case native meshing (no gmsh/OCC)" begin
    # box-assembly cases have EXACT analytic volumes; check representatives hard.
    exact_vol = Dict("9.1"=>2800.0, "5.3"=>1250.0, "5.7"=>1000.0)
    # multi-region conforming cases: interface faces shared by ≤2 tets, every region filled.
    conforming = Set(["5.1","6.6","9.2","5.3","5.4","6.5","7.1","7.2","8.1","8.2","8.3","8.4","10.1"])
    for id in ORDER
        @testset "case $id ($(CASE_NAMES[id]))" begin
            m = build_case(id)
            @test validate(m).ok                       # positive-volume, manifold, non-degenerate
            @test ntets(m) > 0
            @test watertight(m)                        # boundary is a closed manifold (incidence 2)
            if haskey(exact_vol, id)
                @test vol(m) ≈ exact_vol[id] rtol=1e-6 # native box CSG ⇒ exact volume
            end
            if id in conforming
                tpr = tets_per_region(m)
                @test all(v -> v > 0, values(tpr))     # every material region genuinely filled
                ftc = Dict{NTuple{3,Int32},Int}()
                for t in 1:ntets(m), k in 1:4
                    vs=(m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t])
                    f=Tuple(sort(Int32[vs[j] for j in 1:4 if j!=k])); ftc[f]=get(ftc,f,0)+1
                end
                @test maximum(values(ftc)) <= 2        # conforming: no face shared by >2 tets
            end
        end
    end
end
