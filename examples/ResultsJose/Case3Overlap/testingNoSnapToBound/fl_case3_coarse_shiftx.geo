ht=250;
Point(9)= {0., 2.25, 0., ht};
Point(10)= {0.05, 2.25, 0., ht};
Point(11)= {0.95, 2.25, 0., ht};
Point(12)= {1., 2.25, 0., ht};
Point(13)= {0., 0., 0., ht};
Point(14)= {0.05, 0., 0., ht};
Point(15)= {0.95, 0., 0., ht};
Point(16)= {1., 0., 0., ht};
Point(17)= {0., 2.25, 0.333333, ht};
Point(18)= {0.05, 2.25, 0.333333, ht};
Point(19)= {0.95, 2.25, 0.333333, ht};
Point(20)= {1., 2.25, 0.333333, ht};
Point(21)= {0., 0., 0.333333, ht};
Point(22)= {0.05, 0., 0.333333, ht};
Point(23)= {0.95, 0., 0.333333, ht};
Point(24)= {1., 0., 0.333333, ht};
Point(25)= {0., 2.25, 0.666667, ht};
Point(26)= {0.05, 2.25, 0.666667, ht};
Point(27)= {0.95, 2.25, 0.666667, ht};
Point(28)= {1., 2.25, 0.666667, ht};
Point(29)= {0., 0., 0.666667, ht};
Point(30)= {0.05, 0., 0.666667, ht};
Point(31)= {0.95, 0., 0.666667, ht};
Point(32)= {1., 0., 0.666667, ht};
Point(33)= {0., 2.25, 1., ht};
Point(34)= {0.05, 2.25, 1., ht};
Point(35)= {0.95, 2.25, 1., ht};
Point(36)= {1., 2.25, 1., ht};
Point(37)= {0., 0., 1., ht};
Point(38)= {0.05, 0., 1., ht};
Point(39)= {0.95, 0., 1., ht};
Point(40)= {1., 0., 1., ht};
Line(1)= {9, 10};
Line(2)= {10, 14};
Line(3)= {14, 13};
Line(4)= {13, 9};
Line(5)= {10, 11};
Line(6)= {11, 15};
Line(7)= {15, 14};
Line(8)= {11, 12};
Line(9)= {12, 16};
Line(10)= {16, 15};
Line(11)= {17, 18};
Line(12)= {18, 22};
Line(13)= {22, 21};
Line(14)= {21, 17};
Line(15)= {18, 19};
Line(16)= {19, 23};
Line(17)= {23, 22};
Line(18)= {19, 20};
Line(19)= {20, 24};
Line(20)= {24, 23};
Line(21)= {25, 26};
Line(22)= {26, 30};
Line(23)= {30, 29};
Line(24)= {29, 25};
Line(25)= {26, 27};
Line(26)= {27, 31};
Line(27)= {31, 30};
Line(28)= {27, 28};
Line(29)= {28, 32};
Line(30)= {32, 31};
Line(31)= {33, 34};
Line(32)= {34, 38};
Line(33)= {38, 37};
Line(34)= {37, 33};
Line(35)= {34, 35};
Line(36)= {35, 39};
Line(37)= {39, 38};
Line(38)= {35, 36};
Line(39)= {36, 40};
Line(40)= {40, 39};
Line Loop(1)= {1, 2,3,4};
Surface(1)= {1};
Line Loop(2)= {5, 6,7,-2};
Surface(2)= {2};
Line Loop(3)= {8, 9,10,-6};
Surface(3)= {3};
Line Loop(4)= {11, 12,13,14};
Surface(4)= {4};
Line Loop(5)= {15, 16,17,-12};
Surface(5)= {5};
Line Loop(6)= {18, 19,20,-16};
Surface(6)= {6};
Line Loop(7)= {21, 22,23,24};
Surface(7)= {7};
Line Loop(8)= {25, 26,27,-22};
Surface(8)= {8};
Line Loop(9)= {28, 29,30,-26};
Surface(9)= {9};
Line Loop(10)= {31, 32,33,34};
Surface(10)= {10};
Line Loop(11)= {35, 36,37,-32};
Surface(11)= {11};
Line Loop(12)= {38, 39,40,-36};
Surface(12)= {12};
Line(41)= {10, 18};
Line(42)= {17, 9};
Line(43)= {14, 22};
Line(44)= {13, 21};
Line(45)= {11, 19};
Line(46)= {15, 23};
Line(47)= {12, 20};
Line(48)= {16, 24};
Line(49)= {18, 26};
Line(50)= {25, 17};
Line(51)= {22, 30};
Line(52)= {21, 29};
Line(53)= {19, 27};
Line(54)= {23, 31};
Line(55)= {20, 28};
Line(56)= {24, 32};
Line(57)= {26, 34};
Line(58)= {33, 25};
Line(59)= {30, 38};
Line(60)= {29, 37};
Line(61)= {27, 35};
Line(62)= {31, 39};
Line(63)= {28, 36};
Line(64)= {32, 40};
Line Loop(13)= {1, 41,-11,42};
Surface(13)= {13};
Line Loop(14)= {2, 43,-12,-41};
Surface(14)= {14};
Line Loop(15)= {3, 44,-13,-43};
Surface(15)= {15};
Line Loop(16)= {4, -42,-14,-44};
Surface(16)= {16};
Line Loop(17)= {5, 45,-15,-41};
Surface(17)= {17};
Line Loop(18)= {6, 46,-16,-45};
Surface(18)= {18};
Line Loop(19)= {7, 43,-17,-46};
Surface(19)= {19};
Line Loop(20)= {8, 47,-18,-45};
Surface(20)= {20};
Line Loop(21)= {9, 48,-19,-47};
Surface(21)= {21};
Line Loop(22)= {10, 46,-20,-48};
Surface(22)= {22};
Line Loop(23)= {11, 49,-21,50};
Surface(23)= {23};
Line Loop(24)= {12, 51,-22,-49};
Surface(24)= {24};
Line Loop(25)= {13, 52,-23,-51};
Surface(25)= {25};
Line Loop(26)= {14, -50,-24,-52};
Surface(26)= {26};
Line Loop(27)= {15, 53,-25,-49};
Surface(27)= {27};
Line Loop(28)= {16, 54,-26,-53};
Surface(28)= {28};
Line Loop(29)= {17, 51,-27,-54};
Surface(29)= {29};
Line Loop(30)= {18, 55,-28,-53};
Surface(30)= {30};
Line Loop(31)= {19, 56,-29,-55};
Surface(31)= {31};
Line Loop(32)= {20, 54,-30,-56};
Surface(32)= {32};
Line Loop(33)= {21, 57,-31,58};
Surface(33)= {33};
Line Loop(34)= {22, 59,-32,-57};
Surface(34)= {34};
Line Loop(35)= {23, 60,-33,-59};
Surface(35)= {35};
Line Loop(36)= {24, -58,-34,-60};
Surface(36)= {36};
Line Loop(37)= {25, 61,-35,-57};
Surface(37)= {37};
Line Loop(38)= {26, 62,-36,-61};
Surface(38)= {38};
Line Loop(39)= {27, 59,-37,-62};
Surface(39)= {39};
Line Loop(40)= {28, 63,-38,-61};
Surface(40)= {40};
Line Loop(41)= {29, 64,-39,-63};
Surface(41)= {41};
Line Loop(42)= {30, 62,-40,-64};
Surface(42)= {42};
Surface Loop(1001)= {1, 4,13,14,15,16};
Volume(1001)= {1001};
Surface Loop(1002)= {2, 5,17,18,19,14};
Volume(1002)= {1002};
Surface Loop(1003)= {3, 6,20,21,22,18};
Volume(1003)= {1003};
Surface Loop(1004)= {4, 7,23,24,25,26};
Volume(1004)= {1004};
Surface Loop(1005)= {5, 8,27,28,29,24};
Volume(1005)= {1005};
Surface Loop(1006)= {6, 9,30,31,32,28};
Volume(1006)= {1006};
Surface Loop(1007)= {7, 10,33,34,35,36};
Volume(1007)= {1007};
Surface Loop(1008)= {8, 11,37,38,39,34};
Volume(1008)= {1008};
Surface Loop(1009)= {9, 12,40,41,42,38};
Volume(1009)= {1009};
Physical Volume("k33")=Volume{:};
Physical Surface("noflux")= {1, 15, 16, 2, 19, 3, 22, 21, 23, 26, 27, 30, 31, 10, 35, 36, 11, 39, 12, 42, 41};
Physical Surface("inlet") = {32, 29, 25};
Physical Surface("outlet") = {13, 17, 20, 33, 37, 40};

Transfinite Curve{:}=2;
Transfinite Curve{5, 7, 15, 17, 25, 27, 35, 37}=9;
Transfinite Curve{2, 4, 6, 2, 9, 6, 12, 14, 16, 12, 19, 16, 22, 24, 26, 22, 29, 26, 32, 34, 36, 32, 39, 36}=30;
Transfinite Curve{41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64}=4;
Transfinite Surface{:};
Transfinite Volume{:};
Recombine Surface{:};
Recombine Volume{:};//+
