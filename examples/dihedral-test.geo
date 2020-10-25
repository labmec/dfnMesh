x = -1.;
y5 = 0.;
y6 = 1.;
z = 0.;

Point(1) = {0,0,0};
Point(2) = {1,0,0};
Point(3) = {1,1,0};
Point(4) = {0,1,0};
Point(5) = {x,y5,z};
Point(6) = {x,y6,z};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,1};
Line(6) = {4,5};
// Line(7) = {6,5};

Curve Loop(1) = {1,2,3,4};
Curve Loop(2) = {5,-4,6};
Surface(1) = {1};
Surface(2) = {2};

Transfinite Curve{:} = 2;
Transfinite Surface{:};
Recombine Surface{:};

Physical Surface(1) = {1,2};