Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {0, 1, 0, 1};
Point(4) = {0, 0, 1, 1};
Point(5) = {0.2, 0.2, 0, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line(4) = {1, 4};
Line(5) = {2, 4};
Line(6) = {3, 4};

Curve Loop(1) = {1, 2, 3};
Curve Loop(2) = {1, 5, -4};
Curve Loop(3) = {2, 6, -5};
Curve Loop(4) = {3, 4, -6};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Surface Loop(1) = {-1, 2, 3, 4};
Volume(1) = {1};
Point {5} In Surface {1};

MeshSize {PointsOf{Volume{1};}} = 0.8;
selectedSurfaces[] = {-1, 4};
MeshSize {PointsOf{Surface{selectedSurfaces[{0}]};}} = 0.6;
MeshSize {PointsOf{Curve{2}; Point{4};}} = 0.5;
MeshSize {PointsOf{Line{-4};}} = 0.4;
Characteristic Length {PointsOf{Point{3};}} = 0.2;

Physical Point("vertices", 11) = {1:4};
Physical Volume("domain", 12) = {1};
faceTags[] = {-1, 2};
Physical Point("endpoints", 21) = CombinedBoundary{Line{1:2};};
Physical Curve("face boundary", 22) = Boundary{Surface{-1};};
Physical Surface("skin") = CombinedBoundary{Volume{:};};
Physical Curve("two face rim", 25) =
    CombinedBoundary{Surface{faceTags[]};};
Physical Curve("two face boundary", 26) =
    Boundary{Surface{faceTags[]};};
Physical Point("vertices query", 27) = PointsOf{Volume{1};};
