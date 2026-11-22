module APINoNewP2OwnershipTest
# Order-two node ownership on QuadTriNoNewVerts sweeps, by Gmsh 4.15.2's rule:
# SetOrderN gives an edge node to the lowest-dimensional element that has the
# edge (a Curve line, then a Surface cell), the Volume owning the rest.
# Gmsh's getNodes(dim,tag) per entity agrees node for node with Tessella on
# these sources at orders 1 and 2 (cut down from the api_nonew_*_boundary
# files).
using Test, Tessella
const API=Tessella.API
include("api_nonew_gmsh_carriers.jl")

const _SQUARE="""Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={1,1,0,1};Point(4)={0,1,0,1};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
Transfinite Curve{:}=3;Transfinite Surface{1};Recombine Surface{1};
"""

@testset "QuadTri order-two nodes belong to the lowest entity having their edge" begin
    for policy in ("QuadTriNoNewVerts","QuadTriNoNewVerts RecombLaterals")
        API.initialize()
        try
            mktempdir() do directory
                path=joinpath(directory,"patch.geo")
                write(path,_SQUARE*"sweep[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;$policy;};\n")
                opened=API.open_geo!(path;mesh_dim=0)
                volume=Int(opened.lists["sweep"][2])
                API.mesh.generate(3);API.mesh.set_order(2)
                edges,_=gmsh_carriers(API,volume)
                types,_,_=API.mesh.get_elements(3,volume)
                checked=0
                for msh in types
                    for edge in eachcol(reshape(API.mesh.get_element_edge_nodes(msh,volume,false),3,:))
                        support=Tuple(sort!(Int.(collect(edge[1:2]))))
                        _,_,dim,owner=API.mesh.get_node(edge[3])
                        @test (dim,owner)==get(edges,support,(3,volume))
                        checked+=1
                    end
                end
                @test checked>0
            end
        finally
            API.finalize()
        end
    end
end

end
