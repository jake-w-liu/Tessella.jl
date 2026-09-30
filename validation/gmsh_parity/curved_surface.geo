// Curved-boundary planar surfaces: disk, annulus (curved hole), half-disk.

lc = 0.3;

// Unit disk: two semicircle arcs (built-in Circle caps at Pi).
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {-1, 0, 0, lc};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 2};
Curve Loop(1) = {1, 2};
Plane Surface(1) = {1};

// Annulus at x = 4: outer r = 1, inner r = 0.5 hole.
Point(4) = {4, 0, 0, lc};
Point(5) = {5, 0, 0, lc};
Point(6) = {3, 0, 0, lc};
Point(7) = {4.5, 0, 0, lc};
Point(8) = {3.5, 0, 0, lc};
Circle(3) = {5, 4, 6};
Circle(4) = {6, 4, 5};
Circle(5) = {7, 4, 8};
Circle(6) = {8, 4, 7};
Curve Loop(2) = {3, 4};
Curve Loop(3) = {5, 6};
Plane Surface(2) = {2, 3};

// Circular segment at x = 8: chord plus an arc bulge (center below the chord).
Point(10) = {8, 0, 0, lc};
Point(11) = {9, 0, 0, lc};
Point(12) = {8.5, -0.1, 0, lc};
Line(7) = {10, 11};
Circle(8) = {11, 12, 10};
Curve Loop(4) = {7, 8};
Plane Surface(3) = {4};
