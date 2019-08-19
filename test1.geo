dx = 0.25;

Point(100) = {0.0,-1.0,-0.5,dx};
Point(101) = {0.0, 1.0,-0.5,dx};
Point(102) = {0.0, 1.0, 0.5,dx};
Point(103) = {0.0,-1.0, 0.5,dx};

Point(104) = {4.0,-1.0,-0.5,dx};
Point(105) = {4.0, 1.0,-0.5,dx};
Point(106) = {4.0, 1.0, 0.5,dx};
Point(107) = {4.0,-1.0, 0.5,dx};


Line(1) = {101, 100};
Line(2) = {100, 103};
Line(3) = {103, 102};
Line(4) = {102, 101};
Line(5) = {105, 104};
Line(6) = {104, 107};
Line(7) = {107, 106};
Line(8) = {106, 105};
Line(9) = {103, 107};
Line(10) = {100, 104};
Line(11) = {102, 106};
Line(12) = {101, 105};


// Internal surface

Point(200) = {1.0,-0.5, 0.0,dx};
Point(201) = {1.0, 0.5, 0.0,dx};
Point(202) = {3.0,-0.5, 0.0,dx};
Point(203) = {3.0, 0.5, 0.0,dx};
Line(13) = {201, 200};
Line(14) = {200, 202};
Line(15) = {202, 203};
Line(16) = {203, 201};

// Line loops

Line Loop(17) = {3, 4, 1, 2};
Plane Surface(18) = {17};
Line Loop(19) = {8, 5, 6, 7};
Plane Surface(20) = {19};
Line Loop(21) = {11, 8, -12, -4};
Plane Surface(22) = {21};
Line Loop(23) = {11, -7, -9, 3};
Plane Surface(24) = {23};
Line Loop(25) = {9, -6, -10, 2};
Plane Surface(26) = {25};
Line Loop(27) = {10, -5, -12, 1};
Plane Surface(28) = {27};
Line Loop(29) = {13, 14, 15, 16};
Plane Surface(30) = {29};

// Volume
Surface Loop(31) = {24, 22, 20, 28, 26, 18};
Volume(32) = {31};


Transfinite Line {13, 15} = 8 Using Progression 1;
Transfinite Line {16, 14} = 21 Using Progression 1;
Transfinite Surface {30} = {201, 200, 202, 203};

Surface{30} In Volume {32};