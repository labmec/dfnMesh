#ifndef DFNRAWDATA_H
#define DFNRAWDATA_H

#include "DFNNamespace.h"

struct DFNRawData {
    /// this data structure defines the fractures which will cut the mesh
    // each matrix is dimension (3xn) where n is the number of vertices
    TPZFMatrix<REAL> fpolygonmatrices;
    int fmatid = DFNMaterial::Efracture;
    FracLimit flimit_directives = FracLimit::Etruncated;
};

#endif
