/// @brief: A transfinite coarse mesh for Frac3D Benchmark #4 
/// @authors: Pedro Lima


Point( 1) = {-500,100,-100};

// Extrude point through x to create lines
Extrude {300, 0, 0} {
    Point{1}; 
    Layers{1}; 
}
Extrude {250, 0, 0} {
    Point{2}; 
    Layers{1}; 
}
Extrude {300, 0, 0} {
    Point{3}; 
    Layers{1}; 
}

// Extrude lines through y to create surfaces
Extrude{0,300,0}{
    Curve{:};
    Layers{1};
    // Recombine; 
}
Extrude{0,800,000}{
    Curve{4,8,12};
    Layers{4};
    // Recombine; 
}
Extrude{0,300,000}{
    Curve{16,20,24};
    Layers{1};
    // Recombine; 
}
// Extrude surfaces through z to create volumes
Extrude{0,0,200}{
    Surface{:};
    Layers{1};
    // Recombine; 
}
Extrude{0,0,200}{
    Surface{61,83,105,127,149,171,193,215,237};
    Layers{1};
    // Recombine; 
}
Extrude{0,0,200}{
    Surface{259,281,303,325,347,369,391,413,435};
    Layers{1};
    // Recombine; 
}

// Physical Surface(1) = Surface{:};

// Transfinite Surface{:};
// Recombine Surface{:};
// Mesh 2;
// Physical Volume(1) = Volume{:};
// inletbc[] = {70};
// outletbc[] = {14};
// innersurf[] = {27,49};

// Physical Surface("outlet",5) = {14};
// Physical Surface("inlet",4) = {70};

// nofluxbc[] = Surface{:};
// nofluxbc[] -= {outletbc[]};
// nofluxbc[] -= {inletbc[]};
// nofluxbc[] -= {innersurf[]};

// Physical Surface("noflux",6) = {nofluxbc[]};

// Physical Volume("k31",10) = {2, 3};
// Physical Volume("k33",3) = {1};
