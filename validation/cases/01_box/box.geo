// Reference (gmsh) input for the box case. Axis-aligned [0,2]x[0,1]x[0,1], V = 2.
SetFactory("OpenCASCADE");
Box(1) = {0,0,0, 2,1,1};
Mesh.MeshSizeMax = 0.4;
Mesh.MeshSizeMin = 0.4;
