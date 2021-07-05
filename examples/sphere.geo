/// @brief: Semi-spherical domain
/// @authors: Pedro Lima

Point(1) = {0,0,0};

Point(2) = {-1,0,0};
Point(3) = {0,-1,0};
Point(4) = {1,0,0};
Point(5) = {0,1,0};

Point(6) = {0,0,1};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {2, 1, 6};
Circle(6) = {3, 1, 6};
Circle(7) = {4, 1, 6};
Circle(8) = {5, 1, 6};

Curve Loop(1) = {1,2,3,4};
Surface(1) = {1};
Curve Loop(2) = {1, 6, -5};
Surface(2) = {2};
Curve Loop(3) = {2, 7, -6};
Surface(3) = {3};
Curve Loop(4) = {3, 8, -7};
Surface(4) = {4};
Curve Loop(5) = {4, 5, -8};
Surface(5) = {5};

Surface Loop(1) = {1,2,3,4,5};
Volume(1) = {1};

Physical Volume("rock",1) = Volume{:};