/// @brief: Unstructured coarse mesh for Frac3D Benchmark #4 
/// @authors: Pedro Lima

// Mesh size
h = 450;

// Bottom points
Point( 1) = {-500, 100,-100, h};
Point( 2) = { 350, 100,-100, h};
Point( 3) = { 350, 400,-100, h};
Point( 4) = { 350,1500,-100, h};
Point( 5) = {-500,1500,-100, h};
Point( 6) = {-500, 400,-100, h};

// Mid lower points
Point( 7) = {-500, 100, 100, h};
Point( 8) = { 350, 100, 100, h};
Point( 9) = { 350, 400, 100, h};
Point(10) = {-500, 400, 100, h};

// Mid upper points
Point(11) = {-200,1500, 300, h};
Point(12) = {-500,1500, 300, h};
Point(13) = {-500,1200, 300, h};

// Top points
Point(14) = {-500, 100, 500, h};
Point(15) = { 350, 100, 500, h};
Point(16) = { 350,1500, 500, h};
Point(17) = {-200,1500, 500, h};
Point(18) = {-500,1500, 500, h};
Point(19) = {-500,1200, 500, h};

// XY Bottom lines
Line( 1) = { 1, 2};
Line( 2) = { 2, 3};
Line( 3) = { 3, 4};
Line( 4) = { 4, 5};
Line( 5) = { 5, 6};
Line( 6) = { 6, 1};

// XY Mid lower lines
Line( 7) = { 7,10};
Line( 8) = { 9, 8};

// XY Mid upper lines
Line( 9) = {11,12};
Line(10) = {12,13};

// XY Top lines
Line(11) = {14,19};
Line(12) = {19,18};
Line(13) = {18,17};
Line(14) = {17,16};
Line(15) = {16,15};
Line(16) = {15,14};

// Z lower lines
Line(17) = { 1, 7};
Line(18) = { 2, 8};
Line(19) = { 3, 9};
Line(20) = { 4,16};
Line(21) = { 5,12};
Line(22) = { 6,10};

// Z upper lines
Line(23) = { 7,14};
Line(24) = { 8,15};
Line(25) = {11,17};
Line(26) = {12,18};
Line(27) = {13,19};


// NoFlux
// Bottom NoFlux
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};
// Front NoFlux
Curve Loop(2) = {1, 18, 24, 16, -23, -17};
Plane Surface(2) = {2};
// Back NoFlux
Curve Loop(3) = {4, 21, -9, 25, 14, -20};
Plane Surface(3) = {3};
// Right NoFlux
Curve Loop(4) = {3, 20, 15, -24, -8, -19};
Plane Surface(4) = {4};
// Left NoFlux
Curve Loop(5) = {5, 22, -7, 23, 11, -27, -10, -21};
Plane Surface(5) = {5};
// Top NoFlux
Curve Loop(6) = {11, 12, 13, 14, 15, 16};
Plane Surface(6) = {6};

// outlet1
Curve Loop(7) = {6, 17, 7, -22};
Plane Surface(7) = {7};

// outlet2
Curve Loop(8) = {2, 19, 8, -18};
Plane Surface(8) = {8};

// inlet
Curve Loop(9) = {10, 27, 12, -26};
Plane Surface(9) = {9};
Curve Loop(10) = {9, 26, 13, -25};
Plane Surface(10) = {10};




// Volume
Surface Loop(1) = Surface{:};
Volume(1) = {1};



// PHYSICAL GROUPS
inlet[] = {9,10};
outlet1 = {7};
outlet2 = {8};
outlet[] = {outlet1,outlet2};
// noflux[] = {1,2,3,4,5,6};
noflux[] = Surface{:};
noflux[] -= outlet[];
noflux[] -= inlet[];

Physical Surface("NoFlux",-1) = noflux[];
Physical Surface("inlet",-2) = inlet[];
Physical Surface("outlet",-3) = outlet[];

Physical Volume("rock",4) = Volume{:};