SetFactory("OpenCASCADE");
Geometry.OCCBoundsUseStl = 1;
Mesh.RandomSeed = 1;
Mesh.MeshSizeMax = 0.012;
Mesh.MeshSizeMin = 0.00029999999999999997;

sm_air = newv; Box(sm_air) = {-0.040500000000000001, -0.040500000000000001, -0.040500000000000001, 0.28100000000000003, 0.20100000000000001, 0.38100000000000001};
sm_air_inside = newv; Box(sm_air_inside) = {0, 0, 0, 0.22, 0.14000000000000001, 0.29999999999999999};
sm_slot = newv; Box(sm_slot) = {-0.00050000000000000001, 0.002, 0.089999999999999997, 0.00050000000000000001, 0.001, 0.12};
sm_coax_bore = newv; Cylinder(sm_coax_bore) = {0.17000000000000001, 0.14000000000000001, 0.14999999999999999, 0, 0.020500000000000001, 0, 0.0016000000000000001};
sm_case_outer = newv; Box(sm_case_outer) = {-0.00050000000000000001, -0.00050000000000000001, -0.00050000000000000001, 0.221, 0.14099999999999999, 0.30099999999999999};
sm_coax_shield_outer = newv; Cylinder(sm_coax_shield_outer) = {0.17000000000000001, 0.14000000000000001, 0.14999999999999999, 0, 0.020500000000000001, 0, 0.0020999999999999999};
sm_coax_pin = newv; Cylinder(sm_coax_pin) = {0.17000000000000001, 0.1605, 0.14999999999999999, 0, -0.15890000000000001, 0, 0.00080000000000000004};
sm_resistor = news; Rectangle(sm_resistor) = {0.16919999999999999, 0, 0.14999999999999999, 0.0016000000000000001, 0.0016000000000000001};
sm_p1_disk = news; Disk(sm_p1_disk) = {0.17000000000000001, 0, 0.14999999999999999, 0.0016000000000000001};
Rotate { {1, 0, 0}, {0.17000000000000001, 0, 0.14999999999999999}, 1.5707963267948966 } { Surface{sm_p1_disk}; }
Translate { 0, 0.1605, 0 } { Surface{sm_p1_disk}; }
sm_p1_pin_cutout = news; Disk(sm_p1_pin_cutout) = {0.17000000000000001, 0, 0.14999999999999999, 0.00080000000000000004};
Rotate { {1, 0, 0}, {0.17000000000000001, 0, 0.14999999999999999}, 1.5707963267948966 } { Surface{sm_p1_pin_cutout}; }
Translate { 0, 0.1605, 0 } { Surface{sm_p1_pin_cutout}; }

sm_void_tmp[] = BooleanUnion{ Volume{sm_air_inside}; Delete; }{ Volume{sm_slot}; Delete; };
sm_case_void[] = BooleanUnion{ Volume{sm_void_tmp[]}; Delete; }{ Volume{sm_coax_bore}; Delete; };
sm_case_with_shield[] = BooleanUnion{ Volume{sm_case_outer}; Delete; }{ Volume{sm_coax_shield_outer}; Delete; };
sm_case[] = BooleanDifference{ Volume{sm_case_with_shield[]}; Delete; }{ Volume{sm_case_void[]}; Delete; };
sm_p1_surface[] = BooleanDifference{ Surface{sm_p1_disk}; Delete; }{ Surface{sm_p1_pin_cutout}; Delete; };

sm_p1_pr = newp; Point(sm_p1_pr) = {0.1716, 0.1605, 0.14999999999999999, 0.00029999999999999997};
sm_p1_pf = newp; Point(sm_p1_pf) = {0.17080000000000001, 0.1605, 0.14999999999999999, 0.00029999999999999997};
sm_p1_pl = newl; Line(sm_p1_pl) = {sm_p1_pr, sm_p1_pf};

BooleanFragments{ Volume{sm_air, sm_coax_pin, sm_case[]}; Surface{sm_resistor, sm_p1_surface[]}; Delete; }{ Curve{sm_p1_pl}; Point{sm_p1_pr, sm_p1_pf}; Delete; }

// eps_tol = 1.2000000000000002e-06
sm_coax_pin_raw[] = {};
sm_coax_pin_raw[] += Volume In BoundingBox{0.16919880000000001, 0.0015987999999999903, 0.14919879999999999, 0.17080120000000001, 0.16050120000000001, 0.1508012};
sm_coax_pin_t[] = sm_coax_pin_raw[];
sm_resistor_raw[] = {};
sm_resistor_raw[] += Surface In BoundingBox{0.16919879999999998, -1.2000000000000002e-06, 0.14999879999999999, 0.17080119999999999, 0.0016012000000000001, 0.1500012};
sm_resistor_t[] = sm_resistor_raw[];
sm_p1_surface__hole1_raw[] = {};
sm_p1_surface__hole1_raw[] += Surface In BoundingBox{0.16919880000000001, 0.1604988, 0.14919879999999999, 0.17080120000000001, 0.16050120000000001, 0.1508012};
sm_p1_surface__hole1_t[] = sm_p1_surface__hole1_raw[];
sm_p1_surface_raw[] = {};
sm_p1_surface_raw[] += Surface In BoundingBox{0.16839880000000002, 0.1604988, 0.1483988, 0.17160120000000001, 0.16050120000000001, 0.15160119999999999};
sm_p1_surface_t[] = {};
For _ii In {0:#sm_p1_surface_raw[]-1}
  _t_ent = sm_p1_surface_raw[_ii];
  _drop = 0;
  For _jj In {0:#sm_p1_surface__hole1_t[]-1}
    If(sm_p1_surface__hole1_t[_jj] == _t_ent) _drop = 1; EndIf
  EndFor
  If(_drop == 0) sm_p1_surface_t[] += {_t_ent}; EndIf
EndFor
sm_case__hole1_raw[] = {};
sm_case__hole1_raw[] += Volume In BoundingBox{-0.00050120000000000004, -1.2000000000000002e-06, -1.2000000000000002e-06, 0.22000120000000001, 0.16050120000000001, 0.30000119999999997};
sm_case__hole1_t[] = sm_case__hole1_raw[];
sm_case_raw[] = {};
sm_case_raw[] += Volume In BoundingBox{-0.00050120000000000004, -0.00050120000000000004, -0.00050120000000000004, 0.22050120000000001, 0.16050120000000001, 0.30050119999999997};
sm_case_t[] = {};
For _ii In {0:#sm_case_raw[]-1}
  _t_ent = sm_case_raw[_ii];
  _drop = 0;
  For _jj In {0:#sm_coax_pin_t[]-1}
    If(sm_coax_pin_t[_jj] == _t_ent) _drop = 1; EndIf
  EndFor
  For _jj In {0:#sm_case__hole1_t[]-1}
    If(sm_case__hole1_t[_jj] == _t_ent) _drop = 1; EndIf
  EndFor
  If(_drop == 0) sm_case_t[] += {_t_ent}; EndIf
EndFor
sm_p1_line_raw[] = {};
sm_p1_line_raw[] += Curve In BoundingBox{0.17079594950000002, 0.16049594950000001, 0.1499959495, 0.17160405049999999, 0.1605040505, 0.15000405049999999};
sm_p1_line_t[] = sm_p1_line_raw[];
sm_air_raw[] = {};
sm_air_raw[] += Volume In BoundingBox{-0.040501200000000001, -0.040501200000000001, -0.040501200000000001, 0.24050120000000003, 0.16050120000000001, 0.3405012};
sm_air_t[] = {};
For _ii In {0:#sm_air_raw[]-1}
  _t_ent = sm_air_raw[_ii];
  _drop = 0;
  For _jj In {0:#sm_coax_pin_t[]-1}
    If(sm_coax_pin_t[_jj] == _t_ent) _drop = 1; EndIf
  EndFor
  For _jj In {0:#sm_case_t[]-1}
    If(sm_case_t[_jj] == _t_ent) _drop = 1; EndIf
  EndFor
  If(_drop == 0) sm_air_t[] += {_t_ent}; EndIf
EndFor
sm_radiation_raw[] = {};
sm_radiation_raw[] += Surface In BoundingBox{-0.040501200000000001, -0.040501200000000001, -0.040501200000000001, -0.040498800000000001, 0.16050120000000001, 0.3405012};
sm_radiation_raw[] += Surface In BoundingBox{0.24049879999999998, -0.040501200000000001, -0.040501200000000001, 0.2405012, 0.16050120000000001, 0.3405012};
sm_radiation_raw[] += Surface In BoundingBox{-0.040501200000000001, -0.040501200000000001, -0.040501200000000001, 0.2405012, -0.040498800000000001, 0.3405012};
sm_radiation_raw[] += Surface In BoundingBox{-0.040501200000000001, 0.1604988, -0.040501200000000001, 0.2405012, 0.16050120000000001, 0.3405012};
sm_radiation_raw[] += Surface In BoundingBox{-0.040501200000000001, -0.040501200000000001, -0.040501200000000001, 0.2405012, 0.16050120000000001, -0.040498800000000001};
sm_radiation_raw[] += Surface In BoundingBox{-0.040501200000000001, -0.040501200000000001, 0.34049880000000005, 0.2405012, 0.16050120000000001, 0.3405012};
sm_radiation_t[] = {};
For _ii In {0:#sm_radiation_raw[]-1}
  _t_ent = sm_radiation_raw[_ii];
  _drop = 0;
  For _jj In {0:#sm_resistor_t[]-1}
    If(sm_resistor_t[_jj] == _t_ent) _drop = 1; EndIf
  EndFor
  For _jj In {0:#sm_p1_surface_t[]-1}
    If(sm_p1_surface_t[_jj] == _t_ent) _drop = 1; EndIf
  EndFor
  If(_drop == 0) sm_radiation_t[] += {_t_ent}; EndIf
EndFor
sm_coax_pin_pec[] = CombinedBoundary{ Volume{sm_coax_pin_t[]}; };
sm_case_pec[] = CombinedBoundary{ Volume{sm_case_t[]}; };
sm_radiation_t_open[] = {};
For _ii In {0:#sm_radiation_t[]-1}
  _t_ent = sm_radiation_t[_ii];
  _drop = 0;
  For _jj In {0:#sm_coax_pin_pec[]-1}
    If(Abs(sm_coax_pin_pec[_jj]) == Abs(_t_ent)) _drop = 1; EndIf
  EndFor
  For _jj In {0:#sm_case_pec[]-1}
    If(Abs(sm_case_pec[_jj]) == Abs(_t_ent)) _drop = 1; EndIf
  EndFor
  If(_drop == 0) sm_radiation_t_open[] += {_t_ent}; EndIf
EndFor
Physical Volume("air", 1) = {sm_air_t[]};
Physical Volume("coax_pin", 2) = {sm_coax_pin_t[]};
Physical Surface("resistor", 3) = {sm_resistor_t[]};
Physical Volume("case", 4) = {sm_case_t[]};
Physical Surface("p1_surface", 5) = {sm_p1_surface_t[]};
Physical Curve("p1_line", 6) = {sm_p1_line_t[]};
Physical Surface("radiation", 7) = {sm_radiation_t_open[]};
Physical Surface("coax_pin_pecskin", 8) = {sm_coax_pin_pec[]};
Physical Surface("case_pecskin", 9) = {sm_case_pec[]};
Field[1] = Distance; Field[1].SurfacesList = {sm_coax_pin_pec[]};
Field[2] = Threshold; Field[2].InField = 1; Field[2].SizeMin = 0.00053333333333333011; Field[2].SizeMax = 0.012; Field[2].DistMin = 0.00053333333333333011; Field[2].DistMax = 0.0015999999999999903;
Field[3] = Distance; Field[3].SurfacesList = {sm_p1_surface_t[]};
Field[4] = Threshold; Field[4].InField = 3; Field[4].SizeMin = 0.00039999999999999758; Field[4].SizeMax = 0.012; Field[4].DistMin = 0.00039999999999999758; Field[4].DistMax = 0.0031999999999999806;
Field[5] = Distance; Field[5].CurvesList = {sm_p1_line_t[]};
Field[6] = Threshold; Field[6].InField = 5; Field[6].SizeMin = 0.00026669999999999933; Field[6].SizeMax = 0.012; Field[6].DistMin = 0.00026669999999999933; Field[6].DistMax = 0.00080009999999999803;
Field[7] = Box; Field[7].VIn = 0.0089999999999999993; Field[7].VOut = 0.012; Field[7].XMin = -0.018499999999999999; Field[7].XMax = 0.23849999999999999; Field[7].YMin = -0.018499999999999999; Field[7].YMax = 0.17850004999999999; Field[7].ZMin = -0.018499999999999999; Field[7].ZMax = 0.31850000000000001; Field[7].Thickness = 0.017999999999999999;
Field[8] = Min; Field[8].FieldsList = {2, 4, 6, 7};
Background Field = 8;
Mesh.MeshSizeMin = 0.00026669999999999933;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
