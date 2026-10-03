# Native mixed topology: exact local patterns and automatic global tags.
using Tessella
using Tessella.Elements: MixedMesh, ElementBlock, msh_dimension, lagrange_nodes
if !isdefined(Main,:gmsh)
    path=get(ENV,"GMSH_JULIA_API",Sys.iswindows() ?
        joinpath(homedir(),"AppData","Local","Programs","Python","Python312","Lib","gmsh.jl") :
        "/opt/homebrew/lib/gmsh.jl")
    isfile(path) || error("set GMSH_JULIA_API to the pinned Gmsh 4.15.2 binding")
    include(path)
end
gmsh.GMSH_API_VERSION=="4.15.2" || error("Gmsh 4.15.2 required")
let topology=Tessella.MeshEntityTopology
    coords=Float64[0 1 1 0 0 1 1 0 2 2 2 2;
                   0 0 1 1 0 0 1 1 0 1 0 1;
                   0 0 0 0 1 1 1 1 0 0 1 1]
    fixtures=[(3,coords,reshape(Int32[1,2,3,4],4,1)),
        (5,coords,reshape(Int32.(1:8),8,1)),
        (6,coords,reshape(Int32[1,2,4,5,6,8],6,1)),
        (7,coords,reshape(Int32[1,2,3,4,5],5,1)),
        (5,coords,Int32[1 2;2 9;3 10;4 3;5 6;6 11;7 12;8 7])]
    for msh in 8:14
        reference=lagrange_nodes(msh)
        push!(fixtures,(msh,reference,reshape(Int32.(1:size(reference,2)),:,1)))
    end
    gmsh.initialize(String[],false,false)
    try
        gmsh.option.setNumber("General.Terminal",0)
        for (index,(msh,coordinates,cells)) in enumerate(fixtures)
            gmsh.clear();gmsh.model.add("mixed-topology-$index")
            dim=msh_dimension(msh)
            gmsh.model.addDiscreteEntity(dim,1)
            gmsh.model.mesh.addNodes(dim,1,UInt64.(1:size(coordinates,2)),vec(coordinates))
            gmsh.model.mesh.addElementsByType(1,msh,UInt64.(1:size(cells,2)),UInt64.(vec(cells)))
            mesh=MixedMesh(coordinates,[ElementBlock(msh,cells)])
            edges=topology._mesh_edge_topology(mesh)
            faces=topology._mesh_face_topology(mesh)
            gmsh.model.mesh.createEdges();gmsh.model.mesh.createFaces()
            tags,nodes=topology._mesh_all_edges(edges)
            expected,orientation=gmsh.model.mesh.getEdges(nodes)
            expected==tags || error("type $msh: global edge identifiers differ")
            topology._mesh_edges(edges,mesh,nodes)==(expected,orientation) ||
                error("type $msh: edge orientations differ")
            for kind in (3,4)
                ftags,fnodes=topology._mesh_all_faces(faces,kind)
                expected,orientation=gmsh.model.mesh.getFaces(kind,fnodes)
                expected==ftags || error("type $msh: global type-$kind face identifiers differ")
                topology._mesh_faces(faces,mesh,kind,fnodes)==(expected,orientation) ||
                    error("type $msh: face orientations differ")
            end
        end
        println("MIXED_ENTITY_TOPOLOGY_DIFFERENTIAL_OK gmsh=4.15.2 cases=$(length(fixtures))")
    finally
        gmsh.finalize()
    end
end
