// Periodic planar surfaces whose boundaries include a Circle-arc edge.
// Slave surface 2 is master surface 1 translated by +3x; the arc pair becomes
// a derived periodic curve relation.

lc = 0.35;

// Master flag: unit square with the top edge replaced by an arc bulge.
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
Point(5) = {0.5, 0.6, 0, lc};   // arc center (below the chord -> bulge up)
Line(1) = {1, 2};
Line(2) = {2, 3};
Circle(3) = {3, 5, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Slave flag at +3x.
Point(6) = {3, 0, 0, lc};
Point(7) = {4, 0, 0, lc};
Point(8) = {4, 1, 0, lc};
Point(9) = {3, 1, 0, lc};
Point(10) = {3.5, 0.6, 0, lc};
Line(5) = {6, 7};
Line(6) = {7, 8};
Circle(7) = {8, 10, 9};
Line(8) = {9, 6};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

Periodic Surface {2} = {1} Translate {3, 0, 0};
