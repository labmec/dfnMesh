// Gmsh project created on Mon Jul 18 10:53:20 2022
SetFactory("OpenCASCADE");


nely=13;
nelx=6;

//Form z=0 to z=1/3;
nfirsty=2;
//Form z=1/3 to z=2/3;
nsecondy=2;
//Form z=2/3 to z=3/3;
ftthirdy=2;


Point(1)={0,0,0};
Point(2)={1,0,0};
Point(3)={0,2.25,0};
Point(4)={1,2.25,0};

Point(5)={0,0,1};
Point(6)={1,0,1};
Point(7)={0,2.25,1};
Point(8)={1,2.25,1};

Point(9)={0,0,1/3};
Point(10)={1,0,1/3};
Point(11)={0,0,2/3};
Point(12)={1,0,2/3};

Point(13)={0,2.25,1/3};
Point(14)={1,2.25,1/3};
Point(15)={0,2.25,2/3};
Point(16)={1,2.25,2/3};
Line(1) = {2, 4};

Line(2) = {4, 3};
Line(3) = {3, 1};
Line(4) = {1, 2};
Line(5) = {4, 14};
Line(6) = {14, 16};
Line(7) = {16, 8};
Line(8) = {3, 13};
Line(9) = {13, 15};
Line(10) = {15, 7};
Line(11) = {7, 8};
Line(12) = {15, 16};
Line(13) = {13, 14};
Line(14) = {2, 10};
Line(15) = {10, 12};
Line(16) = {12, 6};
Line(17) = {1, 9};
Line(18) = {9, 11};
Line(19) = {11, 5};
Line(20) = {5, 6};
Line(21) = {11, 12};
Line(22) = {9, 10};
Line(23) = {5, 7};
Line(24) = {6, 8};
Line(25) = {16, 12};
Line(26) = {15, 11};
Line(27) = {14, 10};
Line(28) = {13, 9};


Curve Loop(1) = {23, 11, -24, -20};
Plane Surface(1) = {1};
Curve Loop(2) = {26, 21, -25, -12};
Plane Surface(2) = {2};
Curve Loop(3) = {28, 22, -27, -13};
Plane Surface(3) = {3};
Curve Loop(4) = {3, 4, 1, 2};
Plane Surface(4) = {4};
Curve Loop(5) = {7, -24, -16, -25};
Plane Surface(5) = {5};
Curve Loop(6) = {25, -15, -27, 6};
Plane Surface(6) = {6};
Curve Loop(7) = {27, -14, 1, 5};
Plane Surface(7) = {7};
Curve Loop(8) = {23, -10, 26, 19};
Plane Surface(8) = {8};
Curve Loop(9) = {26, -18, -28, 9};
Plane Surface(9) = {9};
Curve Loop(10) = {28, -17, -3, 8};
Plane Surface(10) = {10};
Curve Loop(11) = {20, -16, -21, 19};
Plane Surface(11) = {11};
Curve Loop(12) = {21, -15, -22, 18};
Plane Surface(12) = {12};
Curve Loop(13) = {22, -14, -4, 17};
Plane Surface(13) = {13};
Curve Loop(14) = {11, -7, -12, 10};
Plane Surface(14) = {14};
Curve Loop(15) = {12, -6, -13, 9};
Plane Surface(15) = {15};
Curve Loop(16) = {13, -5, 2, 8};
Plane Surface(16) = {16};
Surface Loop(1) = {1, 8, 14, 5, 11, 2};
Volume(1) = {1};
Surface Loop(2) = {2, 9, 12, 6, 15, 3};
Volume(2) = {2};
Surface Loop(3) = {3, 7, 13, 4, 10, 16};
Volume(3) = {3};

Transfinite Curve{1,3,23,24,25,26,27,28}=nely+1;
Transfinite Curve{2,13,12,11,4,22,21,20}=nelx+1;
Transfinite Curve{7,10,16,19}=nfirsty+1;
Transfinite Curve{6,9,15,18}=nsecondy+1;
Transfinite Curve{5,8,14,17}=ftthirdy+1;
Transfinite Surface{:};
Recombine Surface {:};
Transfinite Volume{:};

Physical Volume("k33") = {1, 2, 3};
Physical Surface("inlet") = {12};
Physical Surface("outlet1") = {14};
Physical Surface("outlet2") = {16};
Physical Surface("noflux") = {1,15, 4, 13, 11, 8, 9, 10, 5, 6, 7};
