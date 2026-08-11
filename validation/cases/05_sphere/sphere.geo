// Reference (gmsh) input for the sphere case. Radius 1.7, true V = 4/3*pi*R^3 =
// 20.5795... gmsh meshes the TRUE sphere; Tessella meshes a twice-subdivided
// octahedron (a polyhedral inscribed approximation).
SetFactory("OpenCASCADE");
Sphere(1) = {0,0,0, 1.7};
Mesh.MeshSizeMax = 0.6;
Mesh.MeshSizeMin = 0.6;
