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

double calcKro(TPZCompEl *cel, Interpol3RelPerm Test,int locPoint);
int main(){
    
   
   
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
    int nref = 5;
    for(int i=0; i<nref; i++){
        int nel = gmesh->NElements();
        for(int j=0; j<nel; j++){
            TPZGeoEl * gel = gmesh->Element(j);
            if (!gel or gel->HasSubElement()) {
                continue;
            }
            gel->Divide(elements);
        }
    }
//    std::ofstream out("gmesh.txt");
//    gmesh->Print(out);
    
    
//    gmesh->Element(0)->Divide(elements);
   // gmesh->BuildConnectivity();
   
    for(int i = 0; i<gmesh->NElements(); i++){
        gmesh->Element(i)->SetMaterialId(1);
    }
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(1,2);
    TPZVec<STATE> fake_sol;
//    TPZL2Projection *mat_0 = new TPZL2Projection(1,2,1,fake_sol);
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    cmesh->SetDefaultOrder(1);
    cmesh->AutoBuild();
    cmesh->SetAllCreateFunctionsContinuous();

 
//    std::ofstream out23("Triangle.vtk");
//    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out23, true);
    
   // cmesh->SetDimModel(2);
//    std::ofstream file("MALLA.txt");
//    cmesh->Print(file);
//    cmesh->Element(1)->Connect(0).fSequenceNumber;
//    TPZBlock<STATE> BLOCK_SOL = cmesh->Block();
//    BLOCK_SOL.Print();
   // cmesh->LoadSolution(sol);
    LoadRelativePermeabilities(cmesh);
    
//    std::ofstream file2("MALLA_sol.txt");
//    cmesh->Print(file2);
    
    
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
        an->DefineGraphMesh(dim,scalar_names,vec_names,plotfile);
        an->PostProcess(div,dim);
        std::cout << "Standard post-processing finished." << std::endl;
    }
    
     
    return 0;
}
void LoadRelativePermeabilities(TPZCompMesh *cmesh){
     int nel = cmesh->NElements();
     int gdl = cmesh->Solution().Rows();
    
   
    
    Interpol3RelPerm Test;
//
//    Test.SetData(Krw.GetData(), Interpol3RelPerm::Kralpha::EKrw);
//    Test.SetData(Krow.GetData(), Interpol3RelPerm::Kralpha::EKrow);
//    Test.SetData(Krg.GetData(), Interpol3RelPerm::Kralpha::EKrg);
//    Test.SetData(Krog.GetData(), Interpol3RelPerm::Kralpha::EKrog);
  
    Test.ReadData("krw.txt", Interpol3RelPerm::Kralpha::EKrw);
    Test.ReadData("krow.txt", Interpol3RelPerm::Kralpha::EKrow);
    Test.ReadData("krg.txt", Interpol3RelPerm::Kralpha::EKrg);
    Test.ReadData("krog.txt", Interpol3RelPerm::Kralpha::EKrog);
    
    
    
    
     TPZFMatrix<STATE> sol;
     sol.Resize(gdl, 1);
    
        for (int i=0; i<nel; i++) {
            TPZCompEl *cel = cmesh->Element(i);
            for(int j=0; j<3; j++){
            
                
                int sec_number = cel->Connect(j).fSequenceNumber;
                int fNElConnected = cel->Connect(j).NElConnected();
                //TEST
                if(sec_number > -1)
                {
                    int64_t pos = cmesh->Block().Position(sec_number);
                        sol(pos,0)=calcKro(cel, Test, j);
                    
                }
                
                //
                
                
            //    sol(sec_number,0)=calcKro(cel, j);
       }
    }

     cmesh->LoadSolution(sol);

}
double calcKro(TPZCompEl *cel, Interpol3RelPerm Test,int locPoint){
    
    TPZFMatrix<REAL> cooridnates;
    cel->Reference()->NodesCoordinates(cooridnates);
    
    
    double x = cooridnates(0,locPoint);
    double y = cooridnates(1,locPoint);
    double Sw = x - 0.5*y;
    double kro = Test.Val(Sw, y);
//    double kro = Test.Deriv(Sw, y).GetVal(1, 1);

    return kro;
    
    // calculat Kro
}
