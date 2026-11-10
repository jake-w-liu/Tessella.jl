
Point(1)={0,0,0,0.4};Point(2)={1,0,0,0.4};Point(3)={1,1,0,0.4};Point(4)={0,1,0,0.4};
Point(5)={0,0,1,0.4};Point(6)={1,0,1,0.4};Point(7)={1,1,1,0.4};Point(8)={0,1,1,0.4};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Line(5)={5,6};Line(6)={6,7};Line(7)={7,8};Line(8)={8,5};
Line(9)={1,5};Line(10)={2,6};Line(11)={3,7};Line(12)={4,8};
Curve Loop(1)={1,2,3,4};Curve Loop(2)={5,6,7,8};
Curve Loop(3)={1,10,-5,-9};Curve Loop(4)={2,11,-6,-10};
Curve Loop(5)={3,12,-7,-11};Curve Loop(6)={4,9,-8,-12};
Plane Surface(1)={1};Plane Surface(2)={2};Plane Surface(3)={3};
Plane Surface(4)={4};Plane Surface(5)={5};Plane Surface(6)={6};
Surface Loop(1)={1,2,3,4,5,6};Volume(1)={1};

Point(9)={0,0,2,0.45};Point(10)={1,0,2,0.45};Point(11)={1,1,2,0.45};Point(12)={0,1,2,0.45};
Line(13)={9,10};Line(14)={10,11};Line(15)={11,12};Line(16)={12,9};
Line(17)={5,9};Line(18)={6,10};Line(19)={7,11};Line(20)={8,12};
Curve Loop(7)={13,14,15,16};Curve Loop(8)={5,18,-13,-17};
Curve Loop(9)={6,19,-14,-18};Curve Loop(10)={7,20,-15,-19};Curve Loop(11)={8,17,-16,-20};
Plane Surface(7)={7};Plane Surface(8)={8};Plane Surface(9)={9};Plane Surface(10)={10};Plane Surface(11)={11};
Surface Loop(102)={2,7,8,9,10,11};Volume(102)={102};
Recombine Surface{1,3,4,5,6};Mesh.Algorithm3D=9;Recombine Volume{1};
