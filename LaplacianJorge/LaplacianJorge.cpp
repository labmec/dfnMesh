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
#include "TPZHybridizeHDiv.h"
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
#include "mixedpoisson.h"
#include "pzbndcond.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TRSRibFrac.h"
#include "TRSRibs.h"
#include "TRSFace.h"
#include "TPZPointCloud.h"
#include "bicgstab.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzsolve.h"
#include "TPZPersistenceManager.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"


#ifdef _AUTODIFF
#include "tfad.h"
#include "fad.h"
#include "pzextractval.h"
#endif

//Creating geometric mesh
TPZGeoMesh * GenerateGmesh(int nx, int ny, double l, double h);

void Ladoderecho (const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
//double RunErrors(int n, int order);
void SolExact(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);


//Creating pressure mesh
TPZCompMesh * GeneratePressureCmesh(TPZGeoMesh *Gmesh, int order);

//Creating flux mesh
TPZCompMesh * GenerateFluxCmesh(TPZGeoMesh *Gmesh, int order);

//Creating mixed mesh
TPZMultiphysicsCompMesh * GenerateMixedCmesh(TPZVec<TPZCompMesh *> fvecmesh, int order);

//Test for the HDiv case
void HDivTest(int nx, int ny, int order1, int order2);

using namespace std;

int main(){

   HDivTest(8, 8, 1, 2);
}

/**
 * @brief Generates a Hdiv test
 * @param nx: number of partions on x
 * @param ny: number of partions on y
 * @param order for the first mesh
 * @param order for the second mesh
 * @return Two computational meshes with different orders
 */
void HDivTest(int nx, int ny, int order1, int order2){
    
//Generating first mesh
    TPZGeoMesh *gmesh = GenerateGmesh(nx, ny, 1, 1);
    TPZCompMesh  *pmesh = GeneratePressureCmesh(gmesh, order1);
    TPZCompMesh *qmesh = GenerateFluxCmesh(gmesh, order1);
    TPZVec<TPZCompMesh *> fmesh(2);
    fmesh[0] = qmesh;
    fmesh[1] = pmesh;
    TPZMultiphysicsCompMesh *MixedMesh = GenerateMixedCmesh(fmesh, 1);
//    MixedMesh->SetDefaultOrder(order1);
    
//Generating second mesh
    TPZCompMesh *qmesh_2 = GenerateFluxCmesh(gmesh, order2);
    TPZCompMesh  *pmesh_2 = GeneratePressureCmesh(gmesh, order2);
    TPZVec<TPZCompMesh *> fmesh_2(2);
    fmesh[0] = qmesh_2;
    fmesh[1] = pmesh_2;
    TPZMultiphysicsCompMesh *MixedMesh_2 = GenerateMixedCmesh(fmesh, 1);
//    MixedMesh_2->SetDefaultOrder(order2);
    
//Solving the system:
    MixedMesh->InitializeBlock();
    bool must_opt_band_width_Q = true;
    int number_threads = 4;
    
//Analysis
    TPZAnalysis *an = new TPZAnalysis(MixedMesh,must_opt_band_width_Q);
    TPZSkylineStructMatrix sparse_matrix(MixedMesh);
    TPZStepSolver<STATE> step;
    sparse_matrix.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    an->SetStructuralMatrix(sparse_matrix);
    an->SetSolver(step);
    an->Assemble();
    an->Solve();
    
//PostProcess
    TPZStack<std::string> scalar, vectors;
    
    TPZManVector<std::string,10> scalnames(2), vecnames(1);
    vecnames[0]  = "Flux";
    scalnames[0] = "Pressure";
    scalnames[1] = "Permeability";
    
    std::ofstream filePrint("MixedHdiv1.txt");
    MixedMesh->Print(filePrint);
    std::string name = "MixedHdiv1.vtk";
  
    an->DefineGraphMesh(2, scalnames, vecnames, name);
    an->PostProcess(0,2);
}

/**
 * @brief Generates the geometric mesh
 * @param nx: number of partions on x
 * @param ny: number of partions on y
 * @param l: lenght
 * @param h: height
 * @return Geometric mesh
 */
TPZGeoMesh * GenerateGmesh(int nx, int ny, double l, double h){
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
  
    TPZVec<int> nels(3,0);
    nels[0]=nx;
    nels[1]=ny;
   
    TPZVec<REAL> x0(3,0.0);
    TPZVec<REAL> x1(3,l);
    x1[1]=h;
    x1[2]=0;
    
//Setting boundary conditions (negative numbers to recognize them)
    TPZGenGrid gen(nels,x0,x1);
    gen.SetElementType(EQuadrilateral);
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, -1);
    gen.SetBC(gmesh, 5, -2);
    gen.SetBC(gmesh, 6, -3);
    gen.SetBC(gmesh, 7, -4);
    return gmesh;
}

/**
 * @brief Generates the pressure computational mesh
 * @param Geometric mesh
 * @param Order
 * @return Pressure computational mesh
 */
TPZCompMesh * GeneratePressureCmesh(TPZGeoMesh *Gmesh, int order){
    
    TPZCompMesh *Cmesh= new TPZCompMesh (Gmesh);
 
    Cmesh->SetDimModel(Gmesh->Dimension());
    Cmesh->SetDefaultOrder(order);
    Cmesh->SetAllCreateFunctionsDiscontinuous();
    Cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
//Add material to the mesh
    int dimen = Gmesh->Dimension();
    int MaterialId = 1;
    STATE Permeability=1;
    
    TPZMatPoisson3d *mat =new TPZMatPoisson3d(MaterialId, dimen);
    
//No convection
    REAL conv=0;
    TPZVec<REAL> convdir(dimen, 0);
    mat->SetParameters(Permeability, conv, convdir);

//Insert material to mesh
    Cmesh->InsertMaterialObject(mat);
    
//Autobuild
    Cmesh->AutoBuild();
    
    int ncon = Cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = Cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    return Cmesh;
}

/**
 * @brief Generates the flux computational mesh
 * @param Geometric mesh
 * @param Order
 * @return Flux computational mesh
 */
TPZCompMesh * GenerateFluxCmesh(TPZGeoMesh *mesh, int order){
    
    int dimen = mesh->Dimension();
    TPZCompMesh *Cmesh = new TPZCompMesh(mesh);
    Cmesh->SetDimModel(dimen);
    Cmesh->SetDefaultOrder(order);
    
//Definition of the approximation space
    int perm=1;
    REAL conv=0;
    REAL perme=1;
    TPZVec<REAL> convdir(dimen , 0.0);
    TPZMatPoisson3d *mat = new TPZMatPoisson3d(perm , dimen);
    mat->SetParameters(perme, conv, convdir);
    
//Inserting volumetric materials objects
    Cmesh->InsertMaterialObject(mat);
    
//Create H(div) functions
    Cmesh->SetAllCreateFunctionsHDiv();
    
//Insert boundary conditions
        int D=0;
        int BC1=-1;
        TPZFMatrix<STATE> val1(1,1,0.0);
        TPZFMatrix<STATE> val2(1,1,0.0);
        TPZMaterial *bc1 = mat->CreateBC(mat, BC1, D, val1, val2);
        Cmesh->InsertMaterialObject(bc1);

        int BC2=-2;
        val2(0,0)=0;
        TPZMaterial *bc2 = mat->CreateBC(mat, BC2, D, val1, val2);
        Cmesh->InsertMaterialObject(bc2);

        int BC3=-3;
        val2(0,0)=0;
        TPZMaterial *bc3 = mat->CreateBC(mat, BC3, D, val1, val2);
        Cmesh->InsertMaterialObject(bc3);

        int BC4=-4;
        val2(0,0)=0;
        TPZMaterial *bc4 = mat->CreateBC(mat, BC4, D, val1, val2);
        Cmesh->InsertMaterialObject(bc4);
    
    Cmesh->AutoBuild();
    
    return Cmesh;
}

/**
 * @brief Generates the mixed computational mesh
 * @param Vector thats contains flux and pressure computational mesh
 * @param Order
 * @return Mixed computational mesh
 */
TPZMultiphysicsCompMesh * GenerateMixedCmesh(TPZVec<TPZCompMesh *> fvecmesh, int order){
    TPZGeoMesh *gmesh = fvecmesh[1]->Reference();
    TPZMultiphysicsCompMesh *MixedMesh = new TPZMultiphysicsCompMesh(gmesh);
    
//Definition of the approximation space
    int dimen= gmesh->Dimension();
    int matnum=1;
    REAL perm=1;
    REAL conv=0;
    TPZVec<REAL> convdir(dimen , 0.0);
    
//Inserting material
    TPZMixedPoisson *mat = new TPZMixedPoisson(matnum, dimen);
   
    mat->SetPermeability(perm);
    mat->SetParameters(perm, conv, convdir);
    
    TPZAutoPointer<TPZFunction<STATE> > sourceterm = new TPZDummyFunction<STATE>(Ladoderecho, 5);

    mat->SetForcingFunction(sourceterm);

//Inserting volumetric materials objects
    MixedMesh->InsertMaterialObject(mat);
    
//Boundary conditions
    int D=0;
    int BC1=-1;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZFMatrix<STATE> val2(1,1,0.0);
    TPZMaterial *bc1 = mat->CreateBC(mat, BC1, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc1);
    
    int BC2=-2;
    val2(0,0)=0;
    TPZMaterial *bc2 = mat->CreateBC(mat, BC2, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc2);
    
    int BC3=-3;
    val2(0,0)=0;
    TPZMaterial *bc3 = mat->CreateBC(mat, BC3, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc3);
    
    int BC4=-4;
    val2(0,0)=0;
    TPZMaterial *bc4 = mat->CreateBC(mat, BC4, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc4);
    
    MixedMesh->SetAllCreateFunctionsMultiphysicElem();
    MixedMesh->SetDimModel(dimen);
    
//Autobild
     TPZManVector<int,5> active_approx_spaces(2); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    MixedMesh->BuildMultiphysicsSpace(active_approx_spaces,fvecmesh);
    
    TPZBuildMultiphysicsMesh::AddElements(fvecmesh, MixedMesh);
    TPZBuildMultiphysicsMesh::AddConnects(fvecmesh,MixedMesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fvecmesh, MixedMesh);
    
    std::cout<<"n connects: "<<MixedMesh->NEquations()<<std::endl;
    std::cout<<"n connects: "<<fvecmesh[0]->NEquations()<<std::endl;
    std::cout<<"n connects: "<<fvecmesh[1]->NEquations()<<std::endl;
    return MixedMesh;
};

/**
 * @brief Generates the force function
 * @param Points values
 * @return Force function value
 */
void Ladoderecho (const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
  
    STATE x = pt[0];
    STATE y = pt[1];
    
//Force function definition
    double fx= -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    
    disp[0]=fx;
 //   disp[0]=0.0;
}
