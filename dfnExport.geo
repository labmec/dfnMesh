//  Geo file generated by Discrete Fracture Network methods 
// Fracture #1 

// POINTS DEFINITION 

h = 0;

Point(1) = {-1,-1,-1, h};
Point(2) = {1,-1,-1, h};
Point(3) = {1,1,-1, h};
Point(4) = {-1,1,-1, h};
Point(5) = {-1,-1,1, h};
Point(6) = {1,-1,1, h};
Point(7) = {1,1,1, h};
Point(8) = {-1,1,1, h};
Point(9) = {-1,-0.424365,0.424365, h};
Point(10) = {-1,-0.424365,-1, h};
Point(11) = {0.209483,-0.209483,-1, h};
Point(12) = {0.209483,-0.209483,0.209483, h};
Point(13) = {0.209483,-0.209483,1, h};
Point(14) = {-1,-0.424365,1, h};
Point(15) = {1,-0.0690355,0.0690355, h};
Point(16) = {1,-0.0690355,-1, h};
Point(17) = {1,-0.0690355,1, h};


// LINES DEFINITION 

Line(8) = {5,2};
Line(9) = {1,5};
Line(12) = {1,2};
Line(13) = {5,6};
Line(15) = {6,2};
Line(16) = {4,8};
Line(20) = {4,3};
Line(22) = {8,3};
Line(23) = {7,3};
Line(25) = {8,7};
Line(44) = {4,9};
Line(45) = {9,5};
Line(46) = {4,10};
Line(47) = {10,1};
Line(48) = {2,11};
Line(49) = {11,4};
Line(50) = {6,12};
Line(51) = {12,4};
Line(52) = {6,13};
Line(53) = {13,8};
Line(54) = {5,14};
Line(55) = {14,8};
Line(56) = {6,15};
Line(57) = {15,3};
Line(58) = {2,16};
Line(59) = {16,3};
Line(60) = {7,17};
Line(61) = {17,6};
Line(90) = {10,9};
Line(91) = {11,9};
Line(92) = {11,10};
Line(93) = {12,9};
Line(94) = {11,12};
Line(95) = {14,9};
Line(96) = {13,12};
Line(97) = {13,14};
Line(98) = {15,16};
Line(99) = {12,15};
Line(100) = {16,11};
Line(101) = {15,13};
Line(102) = {15,17};
Line(103) = {17,13};


// FACES DEFINITION 

Curve Loop(28) = {12,-8,-9};
Surface(28) = {28};
Curve Loop(30) = {-8,13,15};
Surface(30) = {30};
Curve Loop(40) = {16,22,-20};
Surface(40) = {40};
Curve Loop(42) = {25,23,-22};
Surface(42) = {42};
Curve Loop(62) = {-90,-46,44};
Surface(62) = {62};
Curve Loop(63) = {90,45,-9,-47};
Surface(63) = {63};
Curve Loop(64) = {-91,49,44};
Surface(64) = {64};
Curve Loop(65) = {91,45,8,48};
Surface(65) = {65};
Curve Loop(66) = {-92,49,46};
Surface(66) = {66};
Curve Loop(67) = {92,47,12,48};
Surface(67) = {67};
Curve Loop(68) = {-93,51,44};
Surface(68) = {68};
Curve Loop(69) = {93,45,13,50};
Surface(69) = {69};
Curve Loop(70) = {-94,49,-51};
Surface(70) = {70};
Curve Loop(71) = {94,-50,15,48};
Surface(71) = {71};
Curve Loop(72) = {-95,-54,-45};
Surface(72) = {72};
Curve Loop(73) = {95,-44,16,-55};
Surface(73) = {73};
Curve Loop(74) = {-96,-52,50};
Surface(74) = {74};
Curve Loop(75) = {96,51,16,-53};
Surface(75) = {75};
Curve Loop(76) = {-97,53,-55};
Surface(76) = {76};
Curve Loop(77) = {97,-54,13,52};
Surface(77) = {77};
Curve Loop(78) = {-98,57,-59};
Surface(78) = {78};
Curve Loop(79) = {98,-58,-15,56};
Surface(79) = {79};
Curve Loop(80) = {-99,-50,56};
Surface(80) = {80};
Curve Loop(81) = {99,57,-20,-51};
Surface(81) = {81};
Curve Loop(82) = {-100,-58,48};
Surface(82) = {82};
Curve Loop(83) = {100,49,20,-59};
Surface(83) = {83};
Curve Loop(84) = {-101,-56,52};
Surface(84) = {84};
Curve Loop(85) = {101,53,22,-57};
Surface(85) = {85};
Curve Loop(86) = {-102,-56,-61};
Surface(86) = {86};
Curve Loop(87) = {102,-60,23,-57};
Surface(87) = {87};
Curve Loop(88) = {-103,61,52};
Surface(88) = {88};
Curve Loop(89) = {103,53,25,60};
Surface(89) = {89};
Curve Loop(104) = {90,-91,92};
Surface(104) = {104};
Curve Loop(105) = {91,-93,-94};
Surface(105) = {105};
Curve Loop(106) = {93,-95,-97,96};
Surface(106) = {106};
Curve Loop(107) = {94,99,98,100};
Surface(107) = {107};
Curve Loop(108) = {96,99,101};
Surface(108) = {108};
Curve Loop(109) = {-101,102,103};
Surface(109) = {109};

Physical Surface("Fractures",2) = {104,105,106,107,108,109};



// VOLUMES DEFINITION 

Surface Loop(8) = {62,64,66,104};
Volume(8) = {8};
Surface Loop(9) = {28,63,65,67,104};
Volume(9) = {9};
Surface Loop(10) = {64,68,70,105};
Volume(10) = {10};
Surface Loop(11) = {30,65,69,71,105};
Volume(11) = {11};
Surface Loop(12) = {68,73,75,76,106};
Volume(12) = {12};
Surface Loop(13) = {69,72,74,77,106};
Volume(13) = {13};
Surface Loop(14) = {70,78,81,83,107};
Volume(14) = {14};
Surface Loop(15) = {71,79,80,82,107};
Volume(15) = {15};
Surface Loop(16) = {74,80,84,108};
Volume(16) = {16};
Surface Loop(17) = {40,75,81,85,108};
Volume(17) = {17};
Surface Loop(18) = {84,86,88,109};
Volume(18) = {18};
Surface Loop(19) = {42,85,87,89,109};
Volume(19) = {19};


// COARSE ELEMENTS GROUPING

Physical Volume("c1",1) = {8,9};
Physical Volume("c2",2) = {10,11};
Physical Volume("c3",3) = {12,13};
Physical Volume("c4",4) = {14,15};
Physical Volume("c5",5) = {16,17};
Physical Volume("c6",6) = {18,19};



// BOUNDARY CONDITIONS

Physical Surface("bc1",1) = {28,30,40,42};
Physical Surface("bc3",3) = {62,63,66,67,72,73,76,77,78,79,82,83,86,87,88,89};


// OPTIONS

Coherence Mesh;
// Transfinite Surface{:};
// Transfinite Volume{:};
// Recombine Surface{:};
// Recombine Volume{:};