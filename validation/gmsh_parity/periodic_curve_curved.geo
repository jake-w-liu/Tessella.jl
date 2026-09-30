// Curved periodic curve pairs on transfinite strip surfaces: a quarter-arc
// circle Translate pair and a spline Translate pair.

// Strip 1: the master arc dips into the strip; the slave arc is its +3y copy.
Point(1) = {0, 0, 0, 0.5};
Point(2) = {4, 0, 0, 0.5};
Point(3) = {2, -5, 0, 0.5};
Circle(1) = {1, 3, 2};
Point(4) = {0, 3, 0, 0.5};
Point(5) = {4, 3, 0, 0.5};
Point(6) = {2, -2, 0, 0.5};
Circle(2) = {4, 6, 5};
Line(3) = {2, 5};
Line(4) = {4, 1};
Curve Loop(1) = {1, 3, -2, 4};
Plane Surface(1) = {1};
Periodic Curve {2} = {1} Translate {0, 3, 0};
Transfinite Curve {1, 2} = 9;
Transfinite Curve {3, 4} = 5;
Transfinite Surface {1} = {1, 2, 5, 4};

// Strip 2: a transfinite spline pair with a density law.
Point(11) = {0, 6, 0, 0.5};
Point(12) = {4, 6, 0, 0.5};
Point(13) = {1.2, 6.4, 0, 0.5};
Point(14) = {2.8, 6.2, 0, 0.5};
Spline(11) = {11, 13, 14, 12};
Point(15) = {0, 9, 0, 0.5};
Point(16) = {4, 9, 0, 0.5};
Point(17) = {1.2, 9.4, 0, 0.5};
Point(18) = {2.8, 9.2, 0, 0.5};
Spline(12) = {15, 17, 18, 16};
Line(13) = {12, 16};
Line(14) = {15, 11};
Curve Loop(2) = {11, 13, -12, 14};
Plane Surface(2) = {2};
Periodic Curve {12} = {11} Translate {0, 3, 0};
Transfinite Curve {11, 12} = 7 Using Progression 1.3;
Transfinite Curve {13, 14} = 5;
Transfinite Surface {2} = {11, 12, 16, 15};

Mesh 2;
