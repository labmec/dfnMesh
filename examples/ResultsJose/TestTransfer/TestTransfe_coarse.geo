// Gmsh project created on Thu Mar 23 15:42:51 2023
SetFactory("OpenCASCADE");
//+
h=0.0;
Point(1) = {0.0, 0.0, 0, h};
Point(2) = {10.0, 0.0, 0, h};
Point(3) = {0.0, 10.0, 0, h};
Point(4) = {10.0, 10.0, 0, h};

Point(5) = {0.0, 0.0, 15, h};
Point(6) = {10.0, 0.0, 15, h};

Point(7) = {0.0, 10.0, 15, h};
Point(8) = {10.0, 10.0, 15, h};

Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
Line(5) = {2, 6};
//+
Line(6) = {1, 5};
//+
Line(7) = {5, 6};
//+
Line(8) = {6, 8};
//+
Line(9) = {8, 7};
//+
Line(10) = {7, 5};
//+
Line(11) = {7, 3};
//+
Line(12) = {8, 4};

//+
Curve Loop(1) = {1, 5, -7, -6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, 1, 2, 3};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {6, -10, 11, 4};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {2, -12, -8, -5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {7, 8, 9, 10};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {11, -3, -12, 9};
//+
Plane Surface(6) = {6};
//+


Surface Loop(1) = {3, 1, 2, 4, 6, 5};
//+
Volume(1) = {1};

//+
Physical Volume("k33", 1) = {1};

//No Gravity
// //+
// Physical Surface("inlet", 3) = {6};
// //+
// Physical Surface("outlet", 4) = {1};
// //+
// Physical Surface("noflux", 5) = {2, 4, 5, 3};

//Gravity
Physical Surface("inlet", 3) = {5};
//+
Physical Surface("outlet", 13) = {2};
//+
Physical Surface("noflux", 14) = {1, 3, 4, 6};



Coherence Mesh;
Transfinite Curve {:} = 5;
Transfinite Curve {12,11,6,5} = 5;
Transfinite Surface{:};
Transfinite Volume{:};
Coherence Mesh;
Recombine Surface{:};
Recombine Volume{:};//+