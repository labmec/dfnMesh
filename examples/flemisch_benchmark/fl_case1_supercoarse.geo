/// @brief: Supercoarse mesh for Flemisch Benchmark 3D #1 
/// @authors: Pedro Lima
/// @authors: Jose Villegas


Point( 1) = {  0,0,0};

Extrude {100, 0, 0} {
    Point{1}; 
    Layers{1}; 
    Recombine;
}
Extrude {0, 100, 0} {
    Curve{1}; 
    Layers{1}; 
    Recombine;
}


// Base volume
Extrude {0, 0, 10} {
    Surface{:}; 
    Layers{1}; 
    Recombine;
}
// Mid volume
Extrude {0, 0, 80} {
    Surface{27}; 
    Layers{1}; 
    Recombine;
}
// Top volume
Extrude {0, 0, 10} {
    Surface{49}; 
    Layers{1}; 
    Recombine;
}



// Physical Surface(1) = Surface{:};
// Transfinite Curve{:} = 2;
// Transfinite Surface{:};
// Recombine Surface{:};
// Mesh 2;
// Physical Volume(1) = Volume{:};
inletbc[] = {70};
outletbc[] = {14};
innersurf[] = {27,49};

Physical Surface("outlet",5) = {14};
Physical Surface("inlet",4) = {70};

nofluxbc[] = Surface{:};
nofluxbc[] -= {outletbc[]};
nofluxbc[] -= {inletbc[]};
nofluxbc[] -= {innersurf[]};

Physical Surface("noflux",6) = {nofluxbc[]};

Physical Volume("k31",10) = {2, 3};
Physical Volume("k33",3) = {1};
