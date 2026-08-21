SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 2, 1, 1};
Box(2) = {0, 0, 0, 1, 1, 1};
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };
Mesh.MeshSizeMax = 0.5;
Mesh.MeshSizeMin = 0.5;
