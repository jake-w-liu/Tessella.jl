Point(1) = {0, 0, 0, 0.5};
Point(2) = {1, 0, 0, 0.5};
Point(3) = {1, 1, 0, 0.5};
Point(4) = {0, 1, 0, 0.5};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

slaveTag = 29 / 10;
masterTag = 41 / 10;
one = Cos(0);
zero = Atan2(0, 1);
xShift = Sqrt(4) / 2;
Periodic Curve {slaveTag} = {masterTag} Affine
  {one,zero,zero,xShift,
   zero,one,zero,zero,
   zero,zero,one,zero,
   zero,zero,zero,one};
