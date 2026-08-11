// Reference (gmsh) input for the cylinder case. Base (0,0,0), axis +z, height 5,
// radius 2. True volume = pi*R^2*H = 62.8318... gmsh meshes the TRUE circular
// surface; Tessella meshes an inscribed N-gon prism (exact for its own model).
SetFactory("OpenCASCADE");
Cylinder(1) = {0,0,0, 0,0,5, 2};
Mesh.MeshSizeMax = 1.0;
Mesh.MeshSizeMin = 1.0;
