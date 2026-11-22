module OCCBooleanTouchingSolidsChecks
using Test,Tessella
const A=Tessella.API

# Solids that touch only at a vertex or along an edge share that vertex or
# edge after a fragment, as OpenCASCADE's General Fuse makes them: Gmsh
# 4.15.2 gives 15 points for two unit boxes meeting at a corner (not 16) and
# 14 points, 23 curves for two meeting along an edge. Their volume meshes
# then meet at one node, and each volume's tetrahedra stay a manifold.
@testset "OCC fragment of boxes touching at a corner or an edge" begin
    for (p2,counts) in (((1,1,1),(15,24,12,2)),((1,1,0),(14,23,12,2)),
                        ((1,0.5,0),(16,26,13,2)))
        A.initialize()
        try
            A.option("General.Terminal",0)
            occ=A.model.occ
            occ.add_box(0,0,0,1,1,1)
            occ.add_box(p2...,1,1,1)
            out,images=occ.fragment([(3,1)],[(3,2)])
            occ.synchronize()
            @test out==[(3,1),(3,2)]
            @test images==[[(3,1)],[(3,2)]]
            @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==counts
            @test occ.get_mass(3,1)≈1 && occ.get_mass(3,2)≈1
            corner=[t for (_,t) in A.model.get_entities(0)
                    if A.model.get_value(0,t,Float64[])==[1.0,1.0,1.0]]
            @test length(corner)==1
            if p2==(1,1,1)
                A.option("Mesh.MeshSizeMax",0.5)
                A.mesh.generate(3)
                tags,xyz=A.mesh.get_nodes()
                at_corner=count(i->xyz[3i-2:3i]==[1.0,1.0,1.0],eachindex(tags))
                @test at_corner==1
                # Gmsh 4.15.2: 159 nodes on the surfaces at this size
                for d in 2:3
                    @test !isempty(A.mesh.get_elements(d)[2])
                end
            end
        finally
            A.finalize()
        end
    end
end
end
