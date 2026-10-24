# `plugin("BoundaryAngles")` — per-surface boundary-vertex interior angles
# into one `<Filename base>_<tag>` SP list view per dim-2 entity. Upstream
# folds the signed angle to (-2π,0]: straight boundary-edge nodes read -π
# and convex corners read -interior/2π·(2π) — a unit square gives four
# -π/2 corners and every other boundary node -π.
module PluginBoundaryAnglesChecks
using Test,Tessella
const A=Tessella.API

@testset "BoundaryAngles emits per-surface SP views" begin
    A.initialize()
    try
        A.model.add("sq")
        for (x,y) in ((0.,0.),(1.,0.),(1.,1.),(0.,1.))
            A.model.occ.add_point(x,y,0.)
        end
        for i in 1:4; A.model.occ.add_line(i,mod(i,4)+1); end
        A.model.occ.add_curve_loop([1,2,3,4])
        A.model.occ.add_plane_surface([1])
        A.model.occ.synchronize()
        A.mesh.generate(2)
        A.plugin.run("BoundaryAngles")
        tags=A.view.get_tags()
        @test length(tags)==1
        types,data=A.view.get_list_data(only(tags))
        @test types==["SP"]
        rows=only(data)
        @test length(rows)%4==0 && length(rows)>=16
        corners=Set{NTuple{2,Float64}}()
        straight=0
        for i in 1:4:length(rows)
            x,y,a=rows[i],rows[i+1],rows[i+3]
            @test rows[i+2]==0.0
            @test -2π<=a<=0
            iscorner=(x in (0.0,1.0) && y in (0.0,1.0))
            if iscorner
                push!(corners,(x,y))
                @test a≈-π/2
            else
                straight+=1
                @test a≈-π
            end
        end
        @test corners==Set([(0.,0.),(1.,0.),(1.,1.),(0.,1.)])
        @test straight>0
        # Upstream registers "Visible" default 0 despite its help text.
        @test A.view.get_number(only(tags),"Visible")==0.0
    finally
        A.finalize()
    end
end

end # module
