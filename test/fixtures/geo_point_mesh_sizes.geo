lc = 4 / 5;

Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};

MeshSize {:} = lc;
finePoints[] = {2, 4};
MeshSize {finePoints[]} = lc / 2;
Characteristic Length {Sqrt(9)} = 3 * lc / 4;

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1:4};
Plane Surface(1) = {1};

Physical Point("corners", 10) = {1:4};
Physical Surface("domain", 20) = {1};
Mesh 2;
