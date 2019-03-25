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
#include "mixedpoisson.h"
#include "pzbndcond.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZHybridizeHDiv.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphasicFlowAnalysis.h"
#include "TPZDarcyFlow.h"
#include "TPZMixedDarcyFlow.h"
#include "TPZPrimalPoisson.h"
#include "TPZGmshReader.h"

using namespace std;
using namespace cv;

TPZGeoMesh *CreateGeoMesh(int order, int nx, int ny);
TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, int order);
TPZCompMesh *CreateCompMeshNewMat(TPZGeoMesh *gmesh, int order);
TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
TPZCompMesh *CMeshMultphysicsNewMat(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);

TPZCompMesh *CreatePressureCmesh(TPZGeoMesh *gmesh, int order, bool newmat);
TPZCompMesh *CreateFluxCmesh(TPZGeoMesh *gmesh, int order, bool newmat);

void H1TestMat();
void MixedTest();
void MixedTestNewMat();

void sourceterm( const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
int main(){
    
  
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    std::string dirname = PZSOURCEDIR;
    std::string grid;
    grid = dirname + "/Projects/DrillPipe/mesh/case_1.msh";
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    geomesh = Geometry.GeometricGmshMesh(grid);
    const std::string name("Drill pipe section");
    geomesh->SetName(name);
    std::ofstream file("case2.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, file);
    
    
    
  //  MixedTest();
    MixedTestNewMat();
    return 0;
    
  
}

void H1TestMat(){
    TPZGeoMesh *gmesh = CreateGeoMesh(1, 10,10);
    std::ofstream out("testprint.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    //   auto cmesh = CreateCompMesh(gmesh,1);
    auto cmesh = CreateCompMeshNewMat(gmesh, 2);
    std::ofstream outcmesh("testprintcmesh.vtk");
    
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, outcmesh);
    
    TPZAnalysis *an = new TPZAnalysis(cmesh,true);
    
#ifdef USING_MKL2
    TPZSymetricSpStructMatrix struct_mat(cmesh);
    struct_mat.SetNumThreads(4);
    an->SetStructuralMatrix(struct_mat);
#else
    TPZSkylineStructMatrix struct_mat(cmesh);
    struct_mat.SetNumThreads(4);
    an->SetStructuralMatrix(struct_mat);
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an->SetSolver(step);
    
    
    // Solving the LS
    an->Assemble();
    an->Solve();
    
    // post-processing step
    {
        
        const int dim = an->Mesh()->Dimension();
        int div = 0;
        std::string plotfile = "h1_approximation.vtk";
        TPZStack<std::string> scalar_names, vec_names;
        scalar_names.push_back("Solution");
        vec_names.push_back("Derivative");
        an->DefineGraphMesh(dim,scalar_names,vec_names,plotfile);
        an->PostProcess(div,dim);
        std::cout << "Standard post-processing finished." << std::endl;
        
    }
}
void MixedTest(){
    
    TPZGeoMesh *gmesh = CreateGeoMesh(10,10,1);
    int flux_order = 1;
    int p_order = 1;
    
    {
#ifdef PZDEBUG
        std::ofstream file("maze.txt");
        gmesh->Print(file);
        
        std::ofstream out("maze.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
#endif
    }
    
    TPZCompMesh *cmesh_flux =CreateFluxCmesh(gmesh,flux_order, false);
    TPZCompMesh *cmesh_presure =CreatePressureCmesh(gmesh,p_order, false);
    
    
    TPZVec<TPZCompMesh *> fmeshvec(2);
    fmeshvec[0]=cmesh_flux;
    fmeshvec[1]=cmesh_presure;
    gmesh->ResetReference();
    
    TPZCompMesh *MixedMesh = CMeshMultphysics(gmesh,fmeshvec);
    
    std::ofstream file("MixedCMesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(MixedMesh, file);
    
    std::ofstream out("MixedCMesh.txt");
    //    MixedMesh->Print(out);
    
    
    
    std::cout << "number of equations = " << MixedMesh->NEquations() << std::endl;
    
    
    
    //Solving the system:
    bool optimizeBandwidth = true;
    MixedMesh->InitializeBlock();
    
//    TPZCompMesh * cmesh_m_Hybrid;
//    TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
//    TPZHybridizeHDiv hybridizer;
//    tie(cmesh_m_Hybrid, meshvector_Hybrid) = hybridizer.Hybridize(MixedMesh, fmeshvec, true, -1.);
//    cmesh_m_Hybrid->InitializeBlock();
    
    bool must_opt_band_width_Q = true;
    int number_threads = 4;
    TPZAnalysis *an = new TPZAnalysis(MixedMesh,must_opt_band_width_Q);
    TPZElementMatrix mat;
    TPZElementMatrix ef;
    
    MixedMesh->Element(0)->CalcStiff(mat, ef);
    mat.Print(std::cout);
    
    //
    TPZSymetricSpStructMatrix sparse_matrix(MixedMesh);
    TPZStepSolver<STATE> step;
    sparse_matrix.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    an->SetStructuralMatrix(sparse_matrix);
    an->SetSolver(step);
    an->Assemble();
    an->Solve();
    
    
    
    //POS
    TPZManVector<std::string,10> scalnames(2), vecnames(1);
    vecnames[0]  = "Flux";
    
    scalnames[0] = "Pressure";
    scalnames[1] = "Permeability";
    
    
    const int dim = an->Mesh()->Dimension();
    int div = 0;
      std::string plotfile = "hdiv_approximation.vtk";
    an->DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an->PostProcess(div,dim);
    std::cout << "Standard post-processing finished." << std::endl;
    
    return 0;
    
}
TPZGeoMesh *CreateGeoMesh(int order, int nx, int ny){
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    // Creating the Geo mesh
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.0);
    x1[0] = 1;
    x1[1] = 1;
    TPZManVector<int,2> nel(2,2);
    nel[0] = nx;
    nel[1] = ny;
    TPZGenGrid gengrid(nel,x0,x1);
    gengrid.SetElementType(EQuadrilateral);
   
    gmesh->SetDimension(2);
    gengrid.Read(gmesh);
    //gengrid.Read(gmesh,2);
    
  //  gengrid.SetBC(TPZGeoMesh *gr, int side, int bc)
    gengrid.SetBC(gmesh, 4, -1);
    gengrid.SetBC(gmesh, 5, -2);
    gengrid.SetBC(gmesh, 6, -3);
    gengrid.SetBC(gmesh, 7, -4);
    
    gmesh->BuildConnectivity();
    return gmesh;
   
}
TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, int order){
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0e-3;
    REAL perm_1 = 1000.0;
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(impervious_mat,dim);
    TPZMatPoisson3d *mat_1 = new TPZMatPoisson3d(permeable_mat,dim);
    
    
//    TPZDarcyFlow *mat_0 = new TPZDarcyFlow(impervious_mat,dim);
//    TPZDarcyFlow *mat_1 = new TPZDarcyFlow(permeable_mat,dim);
//
    
    mat_0->SetParameters(perm_0, conv, convdir);
    mat_1->SetParameters(perm_1, conv, convdir);
    
    
    
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
  //  cmesh->InsertMaterialObject(mat_1);
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    int right_bc_id = -2;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(right_bc);
    
    int left_bc_id = -4;
    val2(0,0) = 100.0;
    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(left_bc);
    
    int bottom_bc_1id = -1;
    val2(0,0) = 0.0;
    TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc_1);
    
    int top_bc_1id = -3;
    val2(0,0) = 0.0;
    TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(top_bc_1);
    
    
    cmesh->SetName("LaberintoTest");
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(order);
    cmesh->AutoBuild();
    cmesh->Print();
    
    std::ofstream file("cmesh_h.txt");
    cmesh->Print(file);
    
    
    return cmesh;
    
    
}
TPZCompMesh *CreateCompMeshNewMat(TPZGeoMesh *gmesh, int order){
    
    int mat_id = 1;
    
    int dim = gmesh->Dimension();
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    TPZDarcyFlow *mat_0 = new TPZDarcyFlow (mat_id);
    
    TPZAutoPointer<TPZFunction<STATE> > sourcetermm = new TPZDummyFunction<STATE>(sourceterm, 5);
    
    mat_0->SetForcingFunction(sourcetermm);
    //    TPZDarcyFlow *mat_0 = new TPZDarcyFlow(impervious_mat,dim);
    //    TPZDarcyFlow *mat_1 = new TPZDarcyFlow(permeable_mat,dim);
    //
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    int right_bc_id = -1;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(right_bc);
    
     right_bc_id = -2;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc1 = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(right_bc1);
    
     right_bc_id = -3;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc2 = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(right_bc2);
     right_bc_id = -4;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc3 = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(right_bc3);
    
    
    cmesh->SetName("LaberintoTest");
    cmesh->SetDefaultOrder(order);
    cmesh->AutoBuild();
    cmesh->SetAllCreateFunctionsContinuous();
    
    
    std::cout<<"N elements: "<<cmesh->NElements()<<std::endl;
    std::cout<<"N elements: "<<gmesh->NElements()<<std::endl;
    
    //  cmesh->InsertMaterialObject(mat_1);
    return cmesh;
}
void sourceterm( const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    STATE x = pt[0];
    STATE y = pt[1];
    STATE z = pt[2];
    double fx= -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    
    disp[0]=fx;
}

TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec){
    
    //Creating computational mesh for multiphysic elements
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int impervious_mat = 1;
  
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0;
 
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    std::cout<<mphysics->NMaterials();
    
    TPZMixedPoisson *mat_0 = new TPZMixedPoisson(impervious_mat,dim);
    mat_0->SetPermeability(perm_0);
    
   // mat->SetPermeability(1.e-3);
    TPZFNMatrix<9,REAL> K(3,3,0.),KInv(3,3,0.);
    K(0,0) = 1.;
    K(1,1) = 1.;
    K(2,2) = 1.;
    KInv(0,0) = 1.;
    KInv(1,1) = 1;
    KInv(2,2) = 1;
    mat_0->SetPermeabilityTensor(K, KInv);
   
    

    
    mat_0->SetParameters(perm_0, conv, convdir);
   
    mphysics->InsertMaterialObject(mat_0);
    
    //Inserir condicoes de contorno
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    
    int bottom_bc_1id = -1;
    val2(0,0) = 0;
    TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
    mphysics->InsertMaterialObject(bottom_bc_1);
    
    int right_bc_id = -2;
    val2(0,0) = 100.0;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    mphysics->InsertMaterialObject(right_bc);
    

    int top_bc_1id = -3;
    val2(0,0) = 0.0;
    TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
    mphysics->InsertMaterialObject(top_bc_1);
    
    
    int left_bc_id = -4;
    val2(0,0) = 10.0;
    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_D, val1, val2);
    mphysics->InsertMaterialObject(left_bc);
    

    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->SetDimModel(gmesh->Dimension());
    mphysics->AutoBuild();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
#ifdef PZDEBUG
    std::ofstream file("cmesh_mphysics.txt");
    mphysics->Print(file);
#endif
    
    return mphysics;
}
TPZCompMesh *CreatePressureCmesh(TPZGeoMesh *gmesh, int pOrder, bool newmat){
    int MatId = 1;
    int dim = gmesh->Dimension();
    REAL perm_0 = 1.0;
   
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    if (newmat==false){
        TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(MatId,dim);
        mat_0->SetParameters(perm_0, conv, convdir);
        
        
        //  inserting volumetric materials objects
        cmesh->InsertMaterialObject(mat_0);
        
    }
    else{
         TPZDarcyFlow *mat_0 = new TPZDarcyFlow(MatId);
        //  inserting volumetric materials objects
        cmesh->InsertMaterialObject(mat_0);
        
    }
    
    
  
    
    cmesh->SetName("Pressure");
    cmesh->AutoBuild();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef PZDEBUG
    std::ofstream out("cmeshPress.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}
TPZCompMesh *CreateFluxCmesh(TPZGeoMesh *gmesh, int pOrder, bool newmat){
    int Mat_Id = 1;
 
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0;
 
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    
    //Definition of the approximation space:
    
    
    if (newmat==false){
        TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(Mat_Id,dim);
        mat_0->SetParameters(perm_0, conv, convdir);
        
        
        //  inserting volumetric materials objects
        cmesh->InsertMaterialObject(mat_0);
        
        
        cmesh->SetAllCreateFunctionsHDiv(); //Creating H(div) functions
        
        
        int type_D = 0;
        int type_N = 1;
        TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
        
        // Insert boundary conditions
        //Neumann boundary conditions (flux = 0)
        int right_bc_id = -2;
        TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_N, val1, val2);
        cmesh->InsertMaterialObject(right_bc);
        
        int left_bc_id = -4;
        TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_D, val1, val2);
        cmesh->InsertMaterialObject(left_bc);
        
        int bottom_bc_1id = -1;
        TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
        cmesh->InsertMaterialObject(bottom_bc_1);
        
        int top_bc_1id = -3;
        TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_D, val1, val2);
        cmesh->InsertMaterialObject(top_bc_1);
        
        
    }
    else{
        TPZDarcyFlow *mat_0 = new TPZDarcyFlow(Mat_Id);
        //  inserting volumetric materials objects
        cmesh->InsertMaterialObject(mat_0);
        
        cmesh->SetAllCreateFunctionsHDiv(); //Creating H(div) functions
        
        
        int type_D = 0;
        int type_N = 1;
        TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
        
        // Insert boundary conditions
        //Neumann boundary conditions (flux = 0)
        int right_bc_id = -2;
        TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_N, val1, val2);
        cmesh->InsertMaterialObject(right_bc);
        
        int left_bc_id = -4;
        TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_D, val1, val2);
        cmesh->InsertMaterialObject(left_bc);
        
        int bottom_bc_1id = -1;
        TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
        cmesh->InsertMaterialObject(bottom_bc_1);
        
        int top_bc_1id = -3;
        TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_D, val1, val2);
        cmesh->InsertMaterialObject(top_bc_1);
        
    }
    
    
    

    cmesh->SetName("LaberintoTest");
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream file("cmesh_flux.txt");
    cmesh->Print(file);
#endif
    
    return cmesh;
}
void MixedTestNewMat(){
    
    TPZGeoMesh *gmesh = CreateGeoMesh(2,10,10);
    int flux_order = 1;
    int p_order = 1;
    
    {
#ifdef PZDEBUG
        std::ofstream file("maze.txt");
        gmesh->Print(file);
        
        std::ofstream out("maze.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
#endif
    }
    
    TPZCompMesh *cmesh_flux =CreateFluxCmesh(gmesh,flux_order, true);
    TPZCompMesh *cmesh_presure =CreatePressureCmesh(gmesh,p_order, true);
    
    
    TPZVec<TPZCompMesh *> fmeshvec(2);
    fmeshvec[0]=cmesh_flux;
    fmeshvec[1]=cmesh_presure;
    gmesh->ResetReference();
    
    TPZCompMesh *MixedMesh = CMeshMultphysicsNewMat(gmesh,fmeshvec);
    
    std::ofstream file("MixedCMesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(MixedMesh, file);
    
    std::ofstream out("MixedCMesh.txt");
    //    MixedMesh->Print(out);
    
    
    
    std::cout << "number of equations = " << MixedMesh->NEquations() << std::endl;
    
    
    
    //Solving the system:
    bool optimizeBandwidth = true;
    MixedMesh->InitializeBlock();
    
//    TPZCompMesh * cmesh_m_Hybrid;
//    TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
//    TPZHybridizeHDiv hybridizer;
//    tie(cmesh_m_Hybrid, meshvector_Hybrid) = hybridizer.Hybridize(MixedMesh, fmeshvec, true, -1.);
//    cmesh_m_Hybrid->InitializeBlock();
//
    bool must_opt_band_width_Q = true;
    int number_threads = 4;
    TPZAnalysis *an = new TPZAnalysis(MixedMesh,must_opt_band_width_Q);
    
    //
    TPZSymetricSpStructMatrix sparse_matrix(MixedMesh);
    TPZStepSolver<STATE> step;
    sparse_matrix.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    an->SetStructuralMatrix(sparse_matrix);
    an->SetSolver(step);
    an->Assemble();

    
    TPZElementMatrix mat;
    TPZElementMatrix ef;
    
    MixedMesh->Element(0)->CalcStiff(mat, ef);
    mat.Print(std::cout);
    
    an->Solve();
    
    
    
    //POS
    TPZManVector<std::string,10> scalnames(2), vecnames(1);
    vecnames[0]  = "Flux";
    
    scalnames[0] = "Pressure";
    scalnames[1] = "Permeability";
    
    
    const int dim = an->Mesh()->Dimension();
    int div = 0;
    std::string plotfile = "hdiv_approximation.vtk";
    an->DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an->PostProcess(div,dim);
    std::cout << "Standard post-processing finished." << std::endl;
    
    return 0;
    
    
}

TPZCompMesh *CMeshMultphysicsNewMat(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec){
    
    //Creating computational mesh for multiphysic elements
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int impervious_mat = 1;
    
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0;
    
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    std::cout<<mphysics->NMaterials();
    
    TPZMixedDarcyFlow *mat_0 = new TPZMixedDarcyFlow(impervious_mat,dim);
    mat_0->SetPermeability(perm_0);
    
    
    TPZAutoPointer<TPZFunction<STATE> > sourcetermm = new TPZDummyFunction<STATE>(sourceterm, 5);
    
    mat_0->SetForcingFunction(sourcetermm);
    
    
 //   mat_0->SetParameters(perm_0, conv, convdir);
    
    mphysics->InsertMaterialObject(mat_0);
    
    //Inserir condicoes de contorno
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    
    int bottom_bc_1id = -1;
    val2(0,0) = 0;
    TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_D, val1, val2);
    mphysics->InsertMaterialObject(bottom_bc_1);
    
    int right_bc_id = -2;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    mphysics->InsertMaterialObject(right_bc);
    
    
    int top_bc_1id = -3;
    val2(0,0) = 0.0;
    TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_D, val1, val2);
    mphysics->InsertMaterialObject(top_bc_1);
    
    
    int left_bc_id = -4;
    val2(0,0) = 0.0;
    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_D, val1, val2);
    mphysics->InsertMaterialObject(left_bc);
    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->SetDimModel(gmesh->Dimension());
    mphysics->AutoBuild();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
#ifdef PZDEBUG
    std::ofstream file("cmesh_mphysics.txt");
    //mphysics->Print(file);
#endif
    
    return mphysics;
}
