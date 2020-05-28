Point(0) = {0,0,0};
Point(1) = {1,0,0};
Point(2) = {0,1,0};
Point(3) = {0,0,1};

// intruders
distortion = 0.;
Point(4) = {0.5,distortion,distortion};
Point(5) = {0.333,0.333,0};
Line(10) = {0,4};
Line(11) = {4,1};
Line(12) = {4,5};
// intruders



// Line(4) = {0,1};
Line(5) = {1,2};
Line(6) = {2,0};

Line(7) = {0,3};
Line(8) = {1,3};
Line(9) = {2,3};

Curve Loop(10) = {10,11,5,6};
Curve Loop(11) = {10,11,8,-7};
Curve Loop(12) = {5,9,-8};
Curve Loop(13) = {6,7,-9};

Plane Surface(10) = {10};
Plane Surface(11) = {11};
Surface(12) = {12};
Surface(13) = {13};

Line {12} In Surface {10};

Surface Loop(14) = {10,11,12,13};
Volume(14) = {14};


Transfinite Curve {:} = 3;
