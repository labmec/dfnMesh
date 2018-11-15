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
#include "Interpol3RelPerm.h"
#include "TRSLinearInterpolator.h"
#include "pzl2projection.h"
//solucao exata

void LoadRelativePermeabilities(TPZCompMesh *cmesh);
double calcKro(TPZCompEl *cel, int locPoint);
int main(){
    
    TRSLinearInterpolator Krw;
    Krw.ReadData2("krw.txt");
    Krw.GetData().Print(std::cout);
    std::cout<<Krw.Val(0.55)<<std::endl;
    
    TRSLinearInterpolator Krow;
    Krow.ReadData2("krow.txt");
    Krow.GetData().Print(std::cout);
    std::cout<<Krow.Val(0.425)<<std::endl;

    TRSLinearInterpolator Krg;
    Krg.ReadData2("krg.txt");
    Krg.GetData().Print(std::cout);
    std::cout<<Krow.Val(0.375)<<std::endl;

    TRSLinearInterpolator Krog;
    Krog.ReadData2("krog.txt");
    Krog.GetData().Print(std::cout);
    std::cout<<Krog.Val(0.375)<<std::endl;
    
    
    Interpol3RelPerm Test;
    
    Test.SetData(Krw.GetData(), Interpol3RelPerm::Kralpha::EKrw);
    Test.SetData(Krow.GetData(), Interpol3RelPerm::Kralpha::EKrow);
    Test.SetData(Krg.GetData(), Interpol3RelPerm::Kralpha::EKrg);
    Test.SetData(Krog.GetData(), Interpol3RelPerm::Kralpha::EKrog);
    Test.Val(0.2, 0.3);
   
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    
    TPZVec<TPZGeoNode> Node(4);
    TPZVec<REAL> nod(3,0);
    nod[0]=0;
    nod[1]=0;
    nod[2]=0;
    Node[0].SetCoord(nod);
    
    nod[0]=1;
    nod[1]=0;
    nod[2]=0;
    Node[1].SetCoord(nod);
    
    nod[0]=0.5;
    nod[1]=1.0;
    nod[2]=0;
    Node[2].SetCoord(nod);
    
    gmesh->NodeVec().Resize(3);
    gmesh->NodeVec()[0]=Node[0];
    gmesh->NodeVec()[1]=Node[1];
    gmesh->NodeVec()[2]=Node[2];
 
    TPZVec<int64_t> cords(4);
    cords[0]=0;
    cords[1]=1;
    cords[2]=2;
   
    int64_t Nels = gmesh->NElements();
    gmesh->CreateGeoElement(ETriangle, cords, 1, Nels);

    
    gmesh->Element(0)->Initialize();
     gmesh->BuildConnectivity();
    TPZVec<TPZGeoEl *> elements;
    int nref = 1;
    for(int i=0; i<nref; i++){
        int nel = gmesh->NElements();
        for(int i=0; i<nel; i++){
            TPZGeoEl * gel = gmesh->Element(i);
            if (!gel or gel->HasSubElement()) {
                continue;
            }
            gel->Divide(elements);
        }
    }
    std::ofstream out("gmesh.txt");
    gmesh->Print(out);
    
    
//    gmesh->Element(0)->Divide(elements);
   // gmesh->BuildConnectivity();
    gmesh->Print();
    for(int i = 0; i<gmesh->NElements(); i++){
        gmesh->Element(i)->SetMaterialId(1);
    }
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetAllCreateFunctionsContinuous();
    TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(1,2);
    TPZVec<STATE> fake_sol;
//    TPZL2Projection *mat_0 = new TPZL2Projection(1,2,1,fake_sol);
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    cmesh->SetDefaultOrder(1);
    cmesh->AutoBuild();
    cmesh->Print(std::cout);
    std::ofstream out23("Triangle.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out23, true);
    
    
     TPZFMatrix<STATE> sol;
    sol.Resize(12, 1);
    sol(0,0) =0.1;
    sol(1,0) =0.2;
    sol(2,0) =0.3;
    sol(3,0) =0.4;
    sol(4,0) =0.5;
    sol(5,0) =0.6;
    sol(6,0) =0.1;
    sol(7,0) =0.2;
    sol(8,0) =0.3;
    sol(9,0) =0.4;
    sol(10,0) =0.5;
    sol(11,0) =0.6;
    
    
    cmesh->AutoBuild();
  
   // cmesh->SetDimModel(2);
    std::ofstream file("MALLA.txt");
    cmesh->Print(file);
//    cmesh->Element(1)->Connect(0).fSequenceNumber;
//    TPZBlock<STATE> BLOCK_SOL = cmesh->Block();
//    BLOCK_SOL.Print();
   // cmesh->LoadSolution(sol);
    LoadRelativePermeabilities(cmesh);
    std::cout<<cmesh->NEquations()<<std::endl;
    std::cout<<cmesh->NConnects()<<std::endl;
    
    
    std::ofstream out2("Triangle.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out2, true);
    
    TPZAnalysis *an = new TPZAnalysis(cmesh,true);
  
//    for(int i=0; i< cmesh->NElements() - 1; i++){
//    cmesh->SetElementSolution(i, sol2);
//    }

   
    {
        const int dim = an->Mesh()->Dimension();
        int div = 0;
        std::string plotfile = "Kro_test.vtk";
        TPZStack<std::string> scalar_names, vec_names;
        
        scalar_names.push_back("Solution");
        vec_names.push_back("MinusKGradU");
        an->DefineGraphMesh(dim,scalar_names,vec_names,plotfile);
        an->PostProcess(div,dim);
        std::cout << "Standard post-processing finished." << std::endl;
    }
    
    

}
void LoadRelativePermeabilities(TPZCompMesh *cmesh){
     int nel = cmesh->NElements();
     int nconects = cmesh->NConnects();
     int gdl = cmesh->Solution().Rows();
    
  
     TPZFMatrix<STATE> sol;
     sol.Resize(gdl, 1);
     TPZFMatrix<STATE> sol2(nconects-nel,1,-1);
    
        for (int i=0; i<nel; i++) {
            TPZCompEl *cel = cmesh->Element(i);
            for(int j=0; j<3; j++){
            
                int sec_number = cel->Reference()->NodeIndex(j);
                std::cout<<"Element: "<<i<<" order: "<<cel->GetgOrder()<<std::endl;
                sol(sec_number,0)=calcKro(cel, j);
       }
    }
    std::cout<<cmesh->Solution().Rows()<<std::endl;
    std::cout<<cmesh->Solution().Cols()<<std::endl;
    std::cout<<cmesh->NElements()<<std::endl;
     cmesh->LoadSolution(sol);

}
double calcKro(TPZCompEl *cel, int locPoint){
    TPZFMatrix<REAL> cooridnates;
    cel->Reference()->NodesCoordinates(cooridnates);
    //Calc Kro
    
    TRSLinearInterpolator Krw;
    Krw.ReadData2("krw.txt");
    Krw.GetData().Print(std::cout);
    std::cout<<Krw.Val(0.55)<<std::endl;
    
    TRSLinearInterpolator Krow;
    Krow.ReadData2("krow.txt");
    Krow.GetData().Print(std::cout);
    std::cout<<Krow.Val(0.425)<<std::endl;
    
    TRSLinearInterpolator Krg;
    Krg.ReadData2("krg.txt");
    Krg.GetData().Print(std::cout);
    std::cout<<Krow.Val(0.375)<<std::endl;
    
    TRSLinearInterpolator Krog;
    Krog.ReadData2("krog.txt");
    Krog.GetData().Print(std::cout);
    std::cout<<Krog.Val(0.375)<<std::endl;
    
    
    Interpol3RelPerm Test;
    
    Test.SetData(Krw.GetData(), Interpol3RelPerm::Kralpha::EKrw);
    Test.SetData(Krow.GetData(), Interpol3RelPerm::Kralpha::EKrow);
    Test.SetData(Krg.GetData(), Interpol3RelPerm::Kralpha::EKrg);
    Test.SetData(Krog.GetData(), Interpol3RelPerm::Kralpha::EKrog);
    Test.Val(0.2, 0.3);
    
    double x = cooridnates(0,locPoint);
    double y = cooridnates(1,locPoint);
    double Sw = x - 0.5*y;
    double kro = Test.Val(Sw, x);
    
    std::cout<<"Element: "<<cel->Reference()->Index()<<std::endl;
        std::cout<<"Node: "<<locPoint<<std::endl;
        std::cout<<"x: "<<x<<std::endl;
        std::cout<<"y: "<<y<<std::endl;
        std::cout<<"Sw: "<<Sw<<std::endl;
        std::cout<<"Sg: "<<y<<std::endl;
    std::cout<<"kro: "<<kro<<std::endl;
   
    
    return kro;
    
    // calculat Kro
}
