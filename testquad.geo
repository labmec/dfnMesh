Point(0) = {-1,-1,0};
Point(1) = {1,-1,0};
Point(2) = {1,1,0};
Point(3) = {-1,1,0};


Line(4) = {0,1};
Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,0};

Line Loop(1) = {4,5,6,7};
Plane Surface(8) = {1};

Physical Curve(1) = {4,5};
Physical Curve(1) = {6,7}; 