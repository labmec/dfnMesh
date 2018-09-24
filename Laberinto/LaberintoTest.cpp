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
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"

using namespace std;
using namespace cv;

// Creating the computational mesh
TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int p_order);


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
   
    const int bcDL = -1;
    const int bcB = -2;
    const int bcDR = -3;
    const int bcDT = -4;

    
    Mat image = imread("normal.png",IMREAD_GRAYSCALE);
//    Mat image = imread("small.png",IMREAD_GRAYSCALE);
//    Mat image = imread("single_quad.png",IMREAD_GRAYSCALE);
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
    gmesh->BuildConnectivity();
    
    {
        
#ifdef PZDEBUG
        std::ofstream file("maze.txt");
        gmesh->Print(file);
        
        std::ofstream out("maze.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
#endif
    }

    //Creando a malla computacional
    int p_order = 2;
    int number_threads = 8;
    bool must_opt_band_width_Q = true;
    TPZCompMesh *cmesh = CMeshH1(gmesh,p_order);
    TPZAnalysis *an = new TPZAnalysis(cmesh,must_opt_band_width_Q);
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix struct_mat(cmesh);
    struct_mat.SetNumThreads(number_threads);
    an->SetStructuralMatrix(struct_mat);
#else
    TPZSkylineStructMatrix struct_mat(cmesh);
    struct_mat.SetNumThreads(number_threads);
    an->SetStructuralMatrix(struct_mat);
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    
    
    // Solving the LS
    an->Assemble();
    an->Solve();
    
    // post-processing step
    {
        const int dim = an->Mesh()->Dimension();
        int div = 2;
        std::string plotfile = "cg_approximation.vtk";
        TPZStack<std::string> scalar_names, vec_names;
        
        scalar_names.push_back("Solution");
        vec_names.push_back("MinusKGradU");
        an->DefineGraphMesh(dim,scalar_names,vec_names,plotfile);
        an->PostProcess(div,dim);
        std::cout << "Standard post-processing finished." << std::endl;
    }
    
    return 0;
}
TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int p_order){
    
    int impervious_mat = 0;
    int permeable_mat = 1;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0e-16;
    REAL perm_1 = 1.0e-13;
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(impervious_mat,dim);
    TPZMatPoisson3d *mat_1 = new TPZMatPoisson3d(permeable_mat,dim);
    mat_0->SetParameters(perm_0, conv, convdir);
    mat_1->SetParameters(perm_1, conv, convdir);
    
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    cmesh->InsertMaterialObject(mat_1);
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    // Insert boundary conditions
    int right_bc_id = -2;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_N, val1, val2);
    cmesh->InsertMaterialObject(right_bc);
    
    int left_bc_id = -4;
    val2(0,0) = 0.0;
    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_N, val1, val2);
    cmesh->InsertMaterialObject(left_bc);
    
    int bottom_bc_id = -1;
    val2(0,0) = 1.0;
    TPZMaterial * bottom_bc = mat_0->CreateBC(mat_0, bottom_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc);
    
    int top_bc_id = -3;
    val2(0,0) = 0.0;
    TPZMaterial * top_bc = mat_0->CreateBC(mat_0, top_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(top_bc);

    cmesh->SetName("LaberintoTest");
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(p_order);
    cmesh->AutoBuild();

#ifdef PZDEBUG
    std::ofstream file("cmesh_h.txt");
    cmesh->Print(file);
#endif
    
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


