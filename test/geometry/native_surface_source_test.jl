module NativeSurfaceSourceChecks
using Test, Tessella
const M=Tessella.Model

# The plane producer rewrites its boundary segments into chains when a
# source vertex lies on a segment (`_gmsh_native_surface_split_constraints!`)
# or a periodic relation forces a parameter into one. Both rewrites used to
# run on every face and reproduce the segment lists they were given, at about
# twenty heap objects per boundary segment. They are now entered only when
# there is something to insert; these checks pin the decision and that the
# rewrite itself is unchanged.
@testset "segment rewrites run only when a vertex must be inserted" begin
    function square()
        source=M._ModelSurfaceSource()
        for (row,(x,y)) in enumerate(((0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)))
            push!(source.x,x);push!(source.y,y);push!(source.tags,UInt64(row))
            push!(source.owners,(0,Int32(row)));push!(source.parameters,nothing)
            push!(source.world,(x,y,0.0));source.positions[UInt64(row)]=row
        end
        for (a,b) in ((1,2),(2,3),(3,4),(4,1))
            push!(source.segments,(a,b,a));push!(source.boundaries,(a,b))
        end
        return source
    end
    scratch=Tessella.Mesh2D._GmshOrientScratch()

    plain=square()
    @test !M._gmsh_native_surface_split_needed(plain,scratch)
    @test M._gmsh_native_surface_split_constraints!(plain)===nothing
    @test plain.segments==[(1,2,1),(2,3,2),(3,4,3),(4,1,4)]
    @test plain.boundaries==[(1,2),(2,3),(3,4),(4,1)]

    # A fifth vertex in the middle of the first segment splits it; one inside
    # the square, on no segment, does not.
    inside=square()
    push!(inside.x,0.5);push!(inside.y,0.5)
    @test !M._gmsh_native_surface_split_needed(inside,scratch)
    onsegment=square()
    push!(onsegment.x,0.25);push!(onsegment.y,0.0)
    @test M._gmsh_native_surface_split_needed(onsegment,scratch)
    M._gmsh_native_surface_split_constraints!(onsegment)
    @test onsegment.segments==[(1,5,1),(5,2,1),(2,3,2),(3,4,3),(4,1,4)]
    @test onsegment.boundaries==[(1,5),(5,2),(2,3),(3,4),(4,1)]

    # Compiled, the look-ahead allocates nothing (the suite runs interpreted,
    # where every call allocates).
    if Base.JLOptions().compile_enabled==1
        M._gmsh_native_surface_split_needed(plain,scratch)
        @test (@allocated M._gmsh_native_surface_split_needed(plain,scratch))==0
    end
end
end
