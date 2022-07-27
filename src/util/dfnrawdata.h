#ifndef DFNRAWDATA_H
#define DFNRAWDATA_H

#include "DFNNamespace.h"

/// Stores data related to each fracture
struct DFNRawData {
    
    /// this data structure defines the fractures which will cut the mesh
    /// each matrix is dimension (3xn) where n is the number of vertices
    TPZFMatrix<REAL> fpolygonmatrices;
    
    /// Fracture matid
    int fmatid = DFNMaterial::Efracture;
    
    /// Directive to define the fracture limit. For instance, should the code extend the fracture to the nearest element of further refine it to match the boundary of the fracture?
    FracLimit flimit_directives = FracLimit::Etruncated;
    
    /// Refines the the borders of the fracture fnrefborder times. CANNOT be used together with fSizeOfElementsTouchFracBorder
    int fnrefborder = -1;
    
    /// Size of the edges that touch a point on the fracture border. If NEGATIVE value, it is NOT USED
    REAL fSizeOfEdgesTouchFracBorder = -1.;
};


#endif
