//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgeoel.h"
#include "pzgeoelbc.h"
#include "tpzgeoblend.h"
#include "tpzgeoelrefpattern.h"

#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzbiharmonic.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeopyramid.h"
#include "TPZGeoLinear.h"
#include "pzgeotetrahedra.h"


#include "pzmaterial.h"
#include "pzelasmat.h"
#include "pzlog.h"

#include "pzgengrid.h"

#include <time.h>
#include <stdio.h>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include <set>
#include <map>
#include <vector>
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "TPZExtendGridDimension.h"
// #ifdef LOG4CXX
// #include <log4cxx/logger.h>
// #include <log4cxx/basicconfigurator.h>
// #include <log4cxx/propertyconfigurator.h>
// #endif

TPZGeoMesh *SlantedWell(int reservoirmatid, int wellmatid)
{
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.);
    x1[2] = 0.;
    TPZManVector<int,2> nelx(2,4);
    nelx[0] = 8;
    TPZGenGrid gengrid(nelx,x0,x1);
    gengrid.SetElementType(ETriangle);

    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gengrid.Read(gmesh,reservoirmatid);
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    }
    TPZExtendGridDimension extend(gmesh,0.5);
    extend.SetElType(1);
    TPZGeoMesh *gmesh3d = extend.ExtendedMesh(4);
    {
        std::ofstream out("gmesh3d.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh3d, out, true);
    }
    {
        int numel = nelx[1];
        for(int el=0; el<numel; el++)
        {
            TPZManVector<int64_t, 3> nodes(2);
            nodes[0] = 2+(nelx[0]+1)*el+el;
            nodes[1] = 2+(nelx[0]+1)*(el+1)+1+el;
            int matid = 2;
            int64_t index;
            gmesh3d->CreateGeoElement(EOned, nodes, matid, index);

        }
        {
            std::ofstream out("gmesh3dbis.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh3d, out, true);
        }
    }
    gmesh3d->BuildConnectivity();
    gRefDBase.InitializeRefPatterns();
    std::set<int> matids;
    matids.insert(2);
    TPZRefPatternTools::RefineDirectional(gmesh3d, matids);
    {
        std::ofstream out("gmesh3dtris.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh3d, out, true);
    }
    {
        int64_t nel = gmesh3d->NElements();
        // all divided elements have material id 3
        // this will be removed when working with reservoirs
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh3d->Element(el);
            if(gel->Father()) gel->SetMaterialId(3);
        }
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh3d->Element(el);
            // all elements that touch the one dimensional well will have wellmatid
            if(gel->Dimension() == 1)
            {
                for(int side=0; side < gel->NSides(); side++)
                {
                    TPZGeoElSide gelside(gel,side);
                    TPZGeoElSide neighbour(gelside.Neighbour());
                    while (neighbour != gelside) {
                        // this avoid changing the material ids of neighbouring 1d elements
                        if(neighbour.Element()->Father())
                        {
                            neighbour.Element()->SetMaterialId(wellmatid);
                        }
                        neighbour = neighbour.Neighbour();
                    }
                }
            }
        }

    }
    {
        std::ofstream out("gmesh3dtris.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh3d, out, true);
    }
    return gmesh3d;
}

/// remove the father elements (the mesh has no hanging sides)
void RemoveFatherElements(TPZGeoMesh *gmesh)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        
        if (!gel) {
            continue;
        }
        int nsubel = gel->NSubElements();
        if (nsubel == 0) {
            continue;
        }
        for (int sub=0; sub<nsubel; sub++) {
            TPZGeoEl *subel = gel->SubElement(sub);
            if (subel)
            {
                subel->SetFather(int64_t(-1));
            }
        }
        gel->RemoveConnectivities();
        delete gel;
    }
}


/// remove the well elements and replace them by the wellsurfacematid
void RemoveWell(TPZGeoMesh *gmesh, int wellmatid, int wellsurfacematid)
{
    int64_t nel = gmesh->NElements();
    int dim = gmesh->Dimension();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) continue;
        if (gel->MaterialId() == wellmatid) {
            int nsides = gel->NSides();
            for (int is=0; is<nsides; is++) {
                if(gel->SideDimension(is) != dim-1) continue;
                TPZGeoElSide neighbour = gel->Neighbour(is);
                if (neighbour.Element()->MaterialId() != wellmatid) {
                    TPZGeoElBC gbc(neighbour,wellsurfacematid);
                }
            }
            gel->RemoveConnectivities();
            delete gel;
        }
    }
}

using namespace pzgeom;
/// switch the elements which touch the wellsurface to blend elements
void SwitchToBlend(TPZGeoMesh *gmesh, int wellsurfacematid)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) continue;
        if (gel->MaterialId() == wellsurfacematid) {
            MElementType tp = gel->Type();
            TPZManVector<int64_t,8> nodeindices(gel->NNodes());
            gel->GetNodeIndices(nodeindices);
            int matid = gel->MaterialId();
            gel->RemoveConnectivities();
            delete gel;
            
            switch (tp) {
                case EOned:
                    new TPZGeoElRefPattern<TPZGeoBlend<TPZGeoLinear> >(nodeindices,matid,*gmesh);
                    break;
                case EQuadrilateral:
                    new TPZGeoElRefPattern<TPZGeoBlend<TPZGeoQuad> >(nodeindices,matid,*gmesh);
                    break;
                case ETriangle:
                    new TPZGeoElRefPattern<TPZGeoBlend<TPZGeoTriangle> >(nodeindices,matid,*gmesh);
                    break;
                case EPrisma:
                    new TPZGeoElRefPattern<TPZGeoBlend<TPZGeoPrism> >(nodeindices,matid,*gmesh);
                    break;
                case ECube:
                    new TPZGeoElRefPattern<TPZGeoBlend<TPZGeoCube> >(nodeindices,matid,*gmesh);
                    break;
                case EPiramide:
                    new TPZGeoElRefPattern<TPZGeoBlend<TPZGeoPyramid> >(nodeindices,matid,*gmesh);
                    break;
                case ETetraedro:
                    new TPZGeoElRefPattern<TPZGeoBlend<TPZGeoTetrahedra> >(nodeindices,matid,*gmesh);
                    break;
                default:
                    break;
            }
        }
    }
    gmesh->BuildConnectivity();
}

/// add ellipse elements for the surface sides aligned with the plane
void AddEllipseElements(TPZGeoMesh *gmesh, int wellsurfacematid, const TPZVec<REAL> &normal)
{
    
}
