//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif


class TPZGeoMesh;

template<class T>
class TPZVec;

/// create a mesh with a slanted well
TPZGeoMesh *SlantedWell(int reservoirmatid, int wellmatid);

/// remove the father elements (the mesh has no hanging sides)
void RemoveFatherElements(TPZGeoMesh *gmesh);

/// remove the well elements and replace them by the wellsurfacematid
void RemoveWell(TPZGeoMesh *gmesh, int wellmatid, int wellsurfacematid);

/// switch the elements which touch the wellsurface to blend elements
void SwitchToBlend(TPZGeoMesh *gmesh, int wellsurfacematid);

/// add ellipse elements for the surface sides aligned with the plane
void AddEllipseElements(TPZGeoMesh *gmesh, int wellsurfacematid, const TPZVec<REAL> &normal);
