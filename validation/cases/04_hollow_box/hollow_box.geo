// Reference (gmsh) input for the hollow-box (Boolean difference) case.
// Outer [-1,2]x[0,3]x[1,5] (V=36) minus interior cavity [0,1]x[1,2]x[2,3] (1) => 35.
SetFactory("OpenCASCADE");
Box(1) = {-1,0,1, 3,3,4};
Box(2) = {0,1,2, 1,1,1};
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };
Mesh.MeshSizeMax = 1.0;
Mesh.MeshSizeMin = 1.0;
