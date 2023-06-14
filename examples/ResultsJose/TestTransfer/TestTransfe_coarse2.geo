// Definición de puntos
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0.15, 0, 0, 1.0};
Point(3) = {0.15, 0.25, 0, 1.0};
Point(4) = {0, 0.25, 0, 1.0};
Point(5) = {0.05, 0.1, 0, 1.0};
Point(6) = {0.1, 0.1, 0, 1.0};
Point(7) = {0.075, 0.15, 0, 1.0};

// Definición de líneas
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 5};

// Definición de superficies
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Line Loop(2) = {5, 6, 7};
Plane Surface(2) = {2};
Line Loop(3) = {1, 5, 7, 4};
Plane Surface(3) = {3};
Line Loop(4) = {2, 6, 7, 3};
Plane Surface(4) = {4};

// Visualización
Physical Surface("top") = {1};
Physical Surface("bottom") = {2};
Physical Surface("front") = {3};
Physical Surface("back") = {4};