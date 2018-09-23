//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgeoel.h"
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

#include <TPZRefPattern.h>

#include "TPZMaterial.h"
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

#include <opencv2/opencv.hpp>

#include "TRSLinearInterpolator.h"
#include "TPZMatLaplacian.h"

using namespace std;
using namespace cv;

// Creating the computational mesh
TPZCompMesh *DarcyCompMesh(TPZGeoMesh *gmesh);


TPZCompMesh *MalhaCompU(TPZGeoMesh * gmesh,int pOrder);
TPZCompMesh *MalhaCompP(TPZGeoMesh * gmesh, int pOrder);

void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);

void PosProcess(TPZAnalysis &an, std::string plotfile);

void RefinamentoUniforme(TPZGeoMesh  *gMesh, int nh);

void TransferFromMeshes(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);

void TransferFromMultiPhysics(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);

void BuildHybridMesh(TPZCompMesh *cmesh, std::set<int> &MaterialIDs, int LagrangeMat, int InterfaceMat);

int main()
{
    const int matIdB = 0;
    const int matIdN = 1;

    const int dirichlet = 0;
    const int neumann = 1;
    
    const int bcDL = -1;
    const int bcB = -2;
    const int bcDR = -3;
    const int bcDT = -4;

    
    Mat image = imread("normal.png",IMREAD_GRAYSCALE);
    int k=0;
    int px=image.size[0];
    int py=image.size[1];
    int p =px*py;
    vector<int> vec(p,0);
    
        for (int i = 0; i<px; ++i) {
            for (int j = py  ; j>0; --j) {
                vec[p-k]=(int)image.at<uchar>(Point(j, i))/255;
                k++;
                
            }
        }
    
    
    // Creating the Geo mesh
        TPZManVector<REAL,3> x0(3,0.),x1(3,px);
        x1[2] = 0.;
        TPZManVector<int,2> nelx(2,py);
        nelx[0] = px;
        TPZGenGrid gengrid(nelx,x0,x1);
        gengrid.SetElementType(EQuadrilateral);
        TPZGeoMesh *gmesh = new TPZGeoMesh;
        gmesh->SetDimension(2);
        gengrid.Read(gmesh);
    //gengrid.Read(gmesh,2);

        //MatsID
        int nels = gmesh->NElements();
    
        for (int i=0; i<nels; i++) {
            TPZGeoEl *gel =gmesh->Element(i);
            gel->SetMaterialId(vec[i]);
            
        }
    
    
    gengrid.SetBC(gmesh, 4, bcDL);
    gengrid.SetBC(gmesh, 5, bcB);
    gengrid.SetBC(gmesh, 6, bcDR);
    gengrid.SetBC(gmesh, 7, bcDT);
    
        {
            std::ofstream out("LaberintoTestt.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
        }
    gmesh->BuildConnectivity();
    
    //Creando a malla computacional
  
    TPZCompMesh *cmesh = DarcyCompMesh(gmesh);
   
   // cmesh->Print();
    //cmesh->AutoBuild();
    std::ofstream out("ComTest.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out, true);
   
//    TPZAnalysis *an = new TPZAnalysis(cmesh);
//    an->Run();
//
    return 0;
}
TPZCompMesh *DarcyCompMesh(TPZGeoMesh *gmesh){
    int nel=gmesh->NElements();
    
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
    TPZFMatrix<STATE> val2s(2, 1, 0.);
    val2s(0, 0) = 10.0; // vx -> 0
    val2s(1, 0) = 0.0; // vy -> 0
    
    //Materiais 1
    TPZMaterial *material = new TPZMatLaplacian;
    TPZMatLaplacian *material2 = new TPZMatLaplacian;

    
    
//    material->SetPermeability(500);
//    material->SetId(0);
    material2->SetPermeability(50);
    material2->SetId(1);
//
//    TPZMaterial *bc1 = material->CreateBC(material, -1, 1, val2s, val1);
//
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
//    cmesh->InsertMaterialObject(material);
//    cmesh->InsertMaterialObject(material2);
//    cmesh->InsertMaterialObject(bc1);
//    cmesh->InsertMaterialObject(material4);
//    cmesh->InsertMaterialObject(material5);
//    cmesh->InsertMaterialObject(material6);

    cmesh->SetAllCreateFunctionsHDiv();
   //  cmesh->InsertMaterialObject(material);
    cmesh->SetName("LaberintoTest");
  //  cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(1);
//    cmesh->SetAllCreateFunctionsContinuous();
//   int nels= cmesh->NElements();
//    cmesh->SetReference(gmesh);
    cmesh->AutoBuild();

    
   
    return cmesh;
}

TPZCompMesh *MalhaCompU(TPZGeoMesh * gmesh,int pOrder){
    
    //comp mesh
    TPZCompMesh *cmesh=new TPZCompMesh(gmesh);
    TPZCompEl::SetgOrder(1);
    cmesh->SetDimModel(2);
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    
    

}
TPZCompMesh *MalhaCompP(TPZGeoMesh * gmesh, int pOrder){
    
}


