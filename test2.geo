dx = 4;

Point(0) = {-1.0,-1.0,-1.0,dx};
Point(1) = {1.0,-1.0,-1.0,dx};
Point(2) = {1.0,1.0,-1.0,dx};
Point(3) = {-1.0,1.0,-1.0,dx};

Point(4) = {-1.0,-1.0,1.0,dx};
Point(5) = {1.0,-1.0,1.0,dx};
Point(6) = {1.0,1.0,1.0,dx};
Point(7) = {-1.0,1.0,1.0,dx};

Point(8) = {0.0,0.0,0.0, dx};

Line(0) = {0,1};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,0};

Line(4) = {0,4};
Line(5) = {1,5};
Line(6) = {2,6};
Line(7) = {3,7};

Line(8) = {4,5};
Line(9) = {5,6};
Line(10) = {6,7};
Line(11) = {7,4};

// test lines
Line(12) = {0,8};
Line(13) = {8,2};
Line(14) = {0,2};

// loops
Line Loop(16) = {0,1,2,3};
Plane Surface(17) = {16};
Line Loop(18) = {0,5,-8,-4};
Plane Surface(19) = {18};
Line Loop(20) = {1,6,-9,-5};
Plane Surface(21) = {20};
Line Loop(22) = {2,7,-10,-6};
Plane Surface(23) = {22};
Line Loop(24) = {3,4,-11,-7};
Plane Surface(25) = {24};
Line Loop(26) = {8,9,10,11};
Plane Surface(27) = {26};

Line Loop(28) = {12,13,-14};
Plane Surface(29) = {28};


Surface Loop(30) = {17,19,21,23,25,27};
Volume(31) = {30};

Line{14} In Surface{17};
Surface{29} In Volume{31};