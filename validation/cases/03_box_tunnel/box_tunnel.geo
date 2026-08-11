// Reference (gmsh) input for the genus-1 box-with-rectangular-tunnel case.
// Outer [1,5]x[1,5]x[1,3] (V=32) minus a through-tunnel [2,4]x[2,4] along z (8) => 24.
SetFactory("OpenCASCADE");
Box(1) = {1,1,1, 4,4,2};
Box(2) = {2,2,0, 2,2,4};   // spans past the outer box in z: a through-hole
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };
Mesh.MeshSizeMax = 1.0;
Mesh.MeshSizeMin = 1.0;
