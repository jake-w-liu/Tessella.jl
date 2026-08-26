Point(1) = {2, 1, 0, 0.5};
Point(2) = {3, 1, 0, 0.5};
Point(3) = {1, 2, 0, 0.5};
Point(4) = {1, 3, 0, 0.5};
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {3, 4};
Line(4) = {3, 1};
Curve Loop(1) = {1, 2, -3, 4};
Plane Surface(1) = {1};

slaveTag = Sqrt(9);
masterTag = Cos(0);
axisZ = Cos(0);
centerX = Cos(0);
centerY = Sqrt(1);
quarterTurn = Pi / 2;
Periodic Curve {slaveTag} = {masterTag}
  Rotate {{Atan2(0, 1), Sin(0), axisZ},
          {centerX, centerY, 0}, quarterTurn};
