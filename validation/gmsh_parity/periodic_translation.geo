// Five-node right/left curve pairing under the translation (1,0,0).
Point(1) = {0, 0, 0, 0.25};
Point(2) = {1, 0, 0, 0.25};
Point(3) = {1, 1, 0, 0.25};
Point(4) = {0, 1, 0, 0.25};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Curve {2, 4} = 5;
Periodic Curve {2} = {4} Translate {1, 0, 0};
