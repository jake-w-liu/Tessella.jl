SetFactory("Built-in");

lc = 0.6;
pointMaximum = 100;
SetMaxTag Point(pointMaximum + 0.9);
point1 = newp;
Point(point1) = {0, 0, 0, lc};
point2 = newp;
Point(point2) = {1, 0, 0, lc};
point3 = newp;
Point(point3) = {0, 1, 0, lc};
point4 = newp;
Point(point4) = {0, 0, 1, lc};

SetMaxTag Curve(2 * pointMaximum);
edge1 = newl;
Line(edge1) = {point1, point2};
edge2 = newl;
Line(edge2) = {point2, point3};
edge3 = newl;
Line(edge3) = {point3, point1};
edge4 = newl;
Line(edge4) = {point1, point4};
edge5 = newl;
Line(edge5) = {point2, point4};
edge6 = newl;
Line(edge6) = {point3, point4};

loop1 = newll;
Curve Loop(loop1) = {edge1, edge2, edge3};
loop2 = newll;
Curve Loop(loop2) = {edge1, edge5, -edge4};
loop3 = newll;
Curve Loop(loop3) = {edge2, edge6, -edge5};
loop4 = newll;
Curve Loop(loop4) = {edge3, edge4, -edge6};

SetMaxTag Surface(4 * pointMaximum);
surface1 = news;
Plane Surface(surface1) = {loop1};
surface2 = news;
Plane Surface(surface2) = {loop2};
surface3 = news;
Plane Surface(surface3) = {loop3};
surface4 = news;
Plane Surface(surface4) = {loop4};

shell = newsl;
Surface Loop(shell) = {surface1, surface2, surface3, surface4};
SetMaxTag Volume(6 * pointMaximum);
volume = newv;
Volume(volume) = {shell};
nextRegion = newreg;

physicalBase = nextRegion;
Physical Point("corners", physicalBase + 1) = {point1:point4};
Physical Curve("edges", physicalBase + 2) = {edge1:edge6};
Physical Surface("boundary", physicalBase + 3) = {surface1:surface4};
Physical Volume("domain", physicalBase + 4) = {volume};

field = newf;
Field[field] = Distance;
Field[field].PointsList = {point1, point2, point3, point4};

allocatorSnapshot[] = {
  newp, newl, newc, newll, newcl, news, newsl, newv, newreg, newf
};
