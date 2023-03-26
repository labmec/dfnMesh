z1 = 10;
z2 = 32.1;
z3 = 67.9;
z4 = 90;

Point( 1) = {  0,0,0};
Point( 2) = { 20,0,0};
Point( 3) = { 50,0,0};
Point( 4) = { 80,0,0};
Point( 5) = {100,0,0};

Point( 6) = {  0,0,z1};
Point( 7) = { 20,0,z1};
Point( 8) = { 50,0,z1};
Point( 9) = { 80,0,z1};
Point(10) = {100,0,z1};

Point(11) = {  0,0,z2};
Point(12) = { 20,0,z2};
Point(13) = { 50,0,z2};
Point(14) = { 80,0,z2};
Point(15) = {100,0,z2};

Point(16) = {  0,0,z3};
Point(17) = { 20,0,z3};
Point(18) = { 50,0,z3};
Point(19) = { 80,0,z3};
Point(20) = {100,0,z3};

Point(21) = {  0,0,z4};
Point(22) = { 20,0,z4};
Point(23) = { 50,0,z4};
Point(24) = { 80,0,z4};
Point(25) = {100,0,z4};

Point(26) = {  0,0,100};
Point(27) = { 20,0,100};
Point(28) = { 50,0,100};
Point(29) = { 80,0,100};
Point(30) = {100,0,100};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};

Line(5) = {6,7};
Line(6) = {7,8};
Line(7) = {8,9};
Line(8) = {9,10};

Line(9) =  {11,12};
Line(10) = {12,13};
Line(11) = {13,14};
Line(12) = {14,15};

Line(13) = {16,17};
Line(14) = {17,18};
Line(15) = {18,19};
Line(16) = {19,20};

Line(17) = {21,22};
Line(18) = {22,23};
Line(19) = {23,24};
Line(20) = {24,25};

Line(21) = {26,27};
Line(22) = {27,28};
Line(23) = {28,29};
Line(24) = {29,30};

Line(25) = {1,6};
Line(26) = {2,7};
Line(27) = {3,8};
Line(28) = {4,9};
Line(29) = {5,10};

Line(30) = {6,11};
Line(31) = {7,12};
Line(32) = {8,13};
Line(33) = {9,14};
Line(34) = {10,15};

Line(35) = {11,16};
Line(36) = {12,17};
Line(37) = {13,18};
Line(38) = {14,19};
Line(39) = {15,20};

Line(40) = {16,21};
Line(41) = {17,22};
Line(42) = {18,23};
Line(43) = {19,24};
Line(44) = {20,25};

Line(45) = {21,26};
Line(46) = {22,27};
Line(47) = {23,28};
Line(48) = {24,29};
Line(49) = {25,30};

Curve Loop(1) = {1,26,-5,-25};
Curve Loop(2) = {2,27,-6,-26};
Curve Loop(3) = {3,28,-7,-27};
Curve Loop(4) = {4,29,-8,-28};

Curve Loop(5) = {5,31,-9,-30};
Curve Loop(6) = {6,32,-10,-31};
Curve Loop(7) = {7,33,-11,-32};
Curve Loop(8) = {8,34,-12,-33};

Curve Loop(9) = {9,36,-13,-35};
Curve Loop(10) = {10,37,-14,-36};
Curve Loop(11) = {11,38,-15,-37};
Curve Loop(12) = {12,39,-16,-38};

Curve Loop(13) = {13,41,-17,-40};
Curve Loop(14) = {14,42,-18,-41};
Curve Loop(15) = {15,43,-19,-42};
Curve Loop(16) = {16,44,-20,-43};

Curve Loop(17) = {17,46,-21,-45};
Curve Loop(18) = {18,47,-22,-46};
Curve Loop(19) = {19,48,-23,-47};
Curve Loop(20) = {20,49,-24,-48};

Surface(1)  = {1};
Surface(2)  = {2};
Surface(3)  = {3};
Surface(4)  = {4};

Surface(5)  = {5};
Surface(6)  = {6};
Surface(7)  = {7};
Surface(8)  = {8};

Surface(9)  = {9};
Surface(10) = {10};
Surface(11) = {11};
Surface(12) = {12};

Surface(13) = {13};
Surface(14) = {14};
Surface(15) = {15};
Surface(16) = {16};

Surface(17) = {17};
Surface(18) = {18};
Surface(19) = {19};
Surface(20) = {20};

// Physical Surface(1) = Surface{:};
Transfinite Curve{:} = 2;
Transfinite Surface{:};
Recombine Surface{:};
// Mesh 2;
Extrude {0, 100, 0} {
    Surface{:}; 
    Layers{7}; 
    Recombine;
}
// Physical Volume(1) = Volume{:};

Physical Surface("outlet",5) = {1, 2, 3, 4};
Physical Surface("inlet",4) = {422};
Physical Surface("noflux",6) = {17, 18, 19, 20, 13, 14, 15, 16, 12, 11, 10, 9, 5, 6, 7, 8, 480, 392, 304, 216, 128, 418, 440, 462, 484, 124, 102, 80, 58, 137, 225, 313, 401, 489, 115, 203, 291, 379, 467, 93, 181, 269, 357, 445, 71, 159, 247, 335, 423, 70, 158, 246, 334};

Physical Volume("k33", 3) = {1, 2, 3, 4};
Physical Volume("k31", 10) = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
