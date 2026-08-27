lc = 9 / 20;
p = 1.9;
c = 10.9;
loop = 30.9;
surface = 40.9;
shell = 50.9;
volume = 60.9;

Point(p + 0) = {Sin(0), Atan2(0, 1), Floor(0.9), lc};
Point(p + 1) = {Sqrt(1), 0, 0, lc};
Point(p + 2) = {1, Max(0, 1), 0, lc};
Point(p + 3) = {0, 1, 0, lc};
Point(p + 4) = {0, 0, Abs(-1), lc};
Point(p + 5) = {1, 0, 1, lc};
Point(p + 6) = {1, 1, 1, lc};
Point(p + 7) = {0, 1, 1, lc};

Line(c + 0) = {p + 0, p + 1};
Line(c + 1) = {p + 1, p + 2};
Line(c + 2) = {p + 2, p + 3};
Line(c + 3) = {p + 3, p + 0};
Line(c + 4) = {p + 4, p + 5};
Line(c + 5) = {p + 5, p + 6};
Line(c + 6) = {p + 6, p + 7};
Line(c + 7) = {p + 7, p + 4};
Line(c + 8) = {p + 0, p + 4};
Line(c + 9) = {p + 1, p + 5};
Line(c + 10) = {p + 2, p + 6};
Line(c + 11) = {p + 3, p + 7};

Curve Loop(loop + 0) = {c:c + 3};
Curve Loop(loop + 1) = {c + 4:c + 7};
Curve Loop(loop + 2) = {c + 0, c + 9, -(c + 4), -(c + 8)};
Curve Loop(loop + 3) = {c + 1, c + 10, -(c + 5), -(c + 9)};
Curve Loop(loop + 4) = {c + 2, c + 11, -(c + 6), -(c + 10)};
Curve Loop(loop + 5) = {c + 3, c + 8, -(c + 7), -(c + 11)};

Plane Surface(surface + 0) = {loop + 0};
Plane Surface(surface + 1) = {loop + 1};
Plane Surface(surface + 2) = {loop + 2};
Plane Surface(surface + 3) = {loop + 3};
Plane Surface(surface + 4) = {loop + 4};
Plane Surface(surface + 5) = {loop + 5};
Surface Loop(shell) = {surface:surface + 5};
Volume(volume) = {shell};

Point(100.9) = {1 / 2, 1 / 2, 1 / 2, lc};
Point{100.9} In Volume{volume};

Physical Point("corners", 70.9) = {p:p + 7};
Physical Point("probe", 71.9) = {100.9};
Physical Curve("edges", 72.9) = {c:c + 11};
Physical Surface("boundary", 73.9) = {surface:surface + 5};
Physical Volume("domain", 74.9) = {volume};
