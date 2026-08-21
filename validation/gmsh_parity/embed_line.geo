// Unit square with an embedded interior segment. Analytic area = 1.
Point(1) = {0, 0, 0, 0.5};
Point(2) = {1, 0, 0, 0.5};
Point(3) = {1, 1, 0, 0.5};
Point(4) = {0, 1, 0, 0.5};
Point(5) = {0.25, 0.5, 0, 0.5};
Point(6) = {0.75, 0.5, 0, 0.5};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Line{5} In Surface{1};
