SetFactory("OpenCASCADE");

lc = 0.45;
point1 = newp;
point1Again = newp;
Point(point1) = {0, 0, 0, lc};
point2 = newp;
Point(point2) = {1, 0, 0, lc};
point3 = newp;
Point(point3) = {1, 1, 0, lc};
point4 = newp;
Point(point4) = {0, 1, 0, lc};
point5 = newp;
Point(point5) = {0, 0, 1, lc};
point6 = newp;
Point(point6) = {1, 0, 1, lc};
point7 = newp;
Point(point7) = {1, 1, 1, lc};
point8 = newp;
Point(point8) = {0, 1, 1, lc};

edge1 = newl;
edge1Again = newc;
Line(edge1) = {point1, point2};
edge2 = newc;
Line(edge2) = {point2, point3};
edge3 = newl;
Line(edge3) = {point3, point4};
edge4 = newc;
Line(edge4) = {point4, point1};
edge5 = newl;
Line(edge5) = {point5, point6};
edge6 = newc;
Line(edge6) = {point6, point7};
edge7 = newl;
Line(edge7) = {point7, point8};
edge8 = newc;
Line(edge8) = {point8, point5};
edge9 = newl;
Line(edge9) = {point1, point5};
edge10 = newc;
Line(edge10) = {point2, point6};
edge11 = newl;
Line(edge11) = {point3, point7};
edge12 = newc;
Line(edge12) = {point4, point8};

loop1 = newll;
loop1Again = newcl;
Curve Loop(loop1) = {edge1, edge2, edge3, edge4};
loop2 = newcl;
Curve Loop(loop2) = {edge5, edge6, edge7, edge8};
loop3 = newll;
Curve Loop(loop3) = {edge1, edge10, -edge5, -edge9};
loop4 = newcl;
Curve Loop(loop4) = {edge2, edge11, -edge6, -edge10};
loop5 = newll;
Curve Loop(loop5) = {edge3, edge12, -edge7, -edge11};
loop6 = newcl;
Curve Loop(loop6) = {edge4, edge9, -edge8, -edge12};

surface1 = news;
Plane Surface(surface1) = {loop1};
surface2 = news;
Plane Surface(surface2) = {loop2};
surface3 = news;
Plane Surface(surface3) = {loop3};
surface4 = news;
Plane Surface(surface4) = {loop4};
surface5 = news;
Plane Surface(surface5) = {loop5};
surface6 = news;
Plane Surface(surface6) = {loop6};

probe1 = newp;
Point(probe1) = {0, 0.5, 0.5, lc};
probe2 = newp;
Point(probe2) = {1, 0.5, 0.5, lc};
Point{probe1} In Surface{surface6};
Point{probe2} In Surface{surface4};

xSlave = surface4;
xMaster = surface6;
Periodic Surface {xSlave} = {xMaster} Translate {1, 0, 0};
ySlave = surface5;
yMaster = surface3;
Periodic Surface {ySlave} = {yMaster} Translate {0, 1, 0};

shell = newsl;
Surface Loop(shell) = {surface1, surface2, surface3, surface4, surface5, surface6};
volume = newv;
Volume(volume) = {shell};
nextRegion = newreg;

physicalBase = nextRegion + 33;
Physical Point("corners", physicalBase + 1) = {point1:point8};
Physical Curve("edges", physicalBase + 2) = {edge1:edge12};
Physical Surface("boundary", physicalBase + 3) = {surface1:surface6};
Physical Volume("domain", physicalBase + 4) = {volume};
Physical Point("face probes", physicalBase + 5) = {probe1, probe2};

field = newf;
fieldAgain = newf;
Field[field] = Distance;
Field[field].PointsList = {probe1, probe2};
nextField = newf;

allocatorSnapshot[] = {
  newp, newl, newc, newll, newcl, news, newsl, newv, newreg, newf
};
