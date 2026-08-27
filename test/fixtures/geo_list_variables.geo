lc[] = {0.45};
pointTags[] = {1:8};

origin[] = {0, 0, 0};
Point(pointTags[0]) = {origin[0], origin[1], origin[2], lc};
Point(pointTags[1]) = {1, 0, 0, lc};
Point(pointTags[2]) = {1, 1, 0, lc};
Point(pointTags[3]) = {0, 1, 0, lc};
Point(pointTags[4]) = {0, 0, 1, lc};
Point(pointTags[5]) = {1, 0, 1, lc};
Point(pointTags[6]) = {1, 1, 1, lc};
Point(pointTags[7]) = {0, 1, 1, lc};

edgeTags[] = {1:12};
Line(edgeTags[0]) = pointTags[{0, 1}];
Line(edgeTags[1]) = {pointTags[{1, 2}]};
Line(edgeTags[2]) = {pointTags[{2, 3}]};
Line(edgeTags[3]) = {pointTags[{3, 0}]};
Line(edgeTags[4]) = {pointTags[{4, 5}]};
Line(edgeTags[5]) = {pointTags[{5, 6}]};
Line(edgeTags[6]) = {pointTags[{6, 7}]};
Line(edgeTags[7]) = {pointTags[{7, 4}]};
Line(edgeTags[8]) = {pointTags[{0, 4}]};
Line(edgeTags[9]) = {pointTags[{1, 5}]};
Line(edgeTags[10]) = {pointTags[{2, 6}]};
Line(edgeTags[11]) = {pointTags[{3, 7}]};

loopTags[] = {1:6};
Curve Loop(loopTags[0]) = edgeTags[{0:3}];
Curve Loop(loopTags[1]) = {edgeTags[{4:7}]};
Curve Loop(loopTags[2]) = {edgeTags[{0, 9}], -edgeTags[{4, 8}]};
Curve Loop(loopTags[3]) = {edgeTags[{1, 10}], -edgeTags[{5, 9}]};
Curve Loop(loopTags[4]) = {edgeTags[{2, 11}], -edgeTags[{6, 10}]};
Curve Loop(loopTags[5]) = {edgeTags[{3, 8}], -edgeTags[{7, 11}]};

surfaceTags[] = {1:6};
Plane Surface(surfaceTags[0]) = {loopTags[0]};
Plane Surface(surfaceTags[1]) = {loopTags[1]};
Plane Surface(surfaceTags[2]) = {loopTags[2]};
Plane Surface(surfaceTags[3]) = {loopTags[3]};
Plane Surface(surfaceTags[4]) = {loopTags[4]};
Plane Surface(surfaceTags[5]) = {loopTags[5]};

probeTags[] = {101, 102};
pointTags[] += {probeTags[]};
Point(probeTags[0]) = {0, 0.5, 0.5, lc};
Point(probeTags[1]) = {1, 0.5, 0.5, lc};
Point{probeTags[0]} In Surface{surfaceTags[5]};
Point{probeTags[1]} In Surface{surfaceTags[3]};

xSlave[] = {999, surfaceTags[3]};
xSlave[] -= {999};
xMaster[] = {surfaceTags[0]};
xMaster[0] = surfaceTags[5];
xShift[] = {1, 0, 0};
Periodic Surface {xSlave[]} = {xMaster[]}
  Translate {xShift[0], xShift[1], xShift[2]};

ySlave[] = {surfaceTags[4]};
yMaster[] = {surfaceTags[2]};
yShift[] = {0, 1, 0};
Periodic Surface {ySlave[]} = {yMaster[]}
  Translate {yShift[0], yShift[1], yShift[2]};

mutationOracle[] = {1, 2, 3};
mutationOracle[1] *= 3;
mutationOracle = 4;
mutationOracle[] += {7};
mutationOracle[] -= {4};
mutationCopy[] = {mutationOracle[], -mutationOracle[{2:1}]};

shellTags[] = {mutationOracle[1] - 2};
Surface Loop(shellTags[0]) = surfaceTags[];
volumeTags[] = {#mutationCopy[] - 4};
Volume(volumeTags[0]) = shellTags[];

physicalTags[] = {60:64};
physicalTags[{0:4}] += {1, 1, 1, 1, 1};
cornerTags[] = {pointTags[{0:7}]};
Physical Point("corners", physicalTags[0]) = cornerTags[];
Physical Curve("edges", physicalTags[1]) = {edgeTags[]};
Physical Surface("boundary", physicalTags[2]) = {surfaceTags[]};
Physical Volume("domain", physicalTags[3]) = {volumeTags[]};
Physical Point("face probes", physicalTags[4]) = {probeTags[]};

fieldTags[] = {201};
Field[fieldTags[0]] = Distance;
Field[fieldTags[0]].PointsList = {probeTags[]};
