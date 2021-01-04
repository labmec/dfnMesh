nx = 8;
ny = 8;
nz = 8;

Point( 1) = {  0,0,0};
Extrude {1, 0, 0} {
    Point{1}; 
    Layers{nx}; 
    // Recombine;
}

Extrude {0, 1, 0} {
    Line{1}; 
    Layers{ny}; 
    Recombine;
}

Extrude {0, 0, 1} {
    Surface{:}; 
    Layers{nz}; 
    Recombine;
}

Physical Volume(1) = Volume{:};