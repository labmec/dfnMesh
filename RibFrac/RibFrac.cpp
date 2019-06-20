
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzintel.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "TPZGeoLinear.h"


#include "TPZMaterial.h"
// #include "pzlog.h"

#include "pzgengrid.h"

#include <stdio.h>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>

#include <set>
#include <map>
#include <vector>
#include <TPZRefPattern.h>
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "TPZExtendGridDimension.h"

#include <opencv2/opencv.hpp>
#include "TRSRibFrac.h"
#include "TRSRibs.h"
#include "TRSFace.h"

//MATERIAL ID MAP
// 1 gmesh
// 4 skeleton
// 5 uncut ribs
// 12 ribs that will be divided
// 20 mid-fracture cut faces
// 35 end-fracture cut faces
// 40 Fracture plane
// 50 divided ribs

using namespace std;

int main(){
  
  // Creating the Geo mesh
	int dimel = 5;
	TPZManVector<REAL, 3> x0(3, 0.), x1(3, 6.0);
	x1[2] = 0.;
	TPZManVector<int, 2> nelx(2, dimel);
	TPZGenGrid gengrid(nelx, x0, x1);
	gengrid.SetElementType(EQuadrilateral);
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	gmesh->SetDimension(2);
	gengrid.Read(gmesh);

	/// Mesh 3D

  TPZExtendGridDimension extend(gmesh, 1);
  extend.SetElType(1);
  TPZGeoMesh *gmesh3d = extend.ExtendedMesh(3);
  gmesh = gmesh3d;
  //std::ofstream out3("3DMESH.vtk");
  //TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out3, true);


	//    Reading coordinates of a plane from txt file
	Matrix plane(3, 4);
	int i = 0;
	int j = 0;
	string value;
	ifstream plane_file("fracture.txt");
	if (!plane_file)
	{
		std::cout << "Error reading file" << endl;
		DebugStop();
	}
	cout << "Fracture plane defined as: \n";
	string line;
	while (getline(plane_file, line))
	{
		std::stringstream ss(line);
		while (getline(ss, value, ' '))
		{
			while (value.length() == 0)
			{
				getline(ss, value, ' ');
			}
			plane(i, j) = std::stod(value);
			cout << plane(i, j) << " , ";
			j++;
		}
		j = 0;
		i++;
		cout << endl;
	}
	std::cout << endl;
    

    //    Setting a plane
    
    TPZVec<int64_t> nodeindices;
    gmesh->Element(0)->GetNodeIndices(nodeindices);
    

    TPZVec<TPZGeoNode> Node(4);
    TPZVec<REAL> nod(3,0);
    nod[0]=plane(0,0);
    nod[1]=plane(1,0);
    nod[2]=plane(2,0);
    Node[0].SetCoord(nod);
    
    nod[0]=plane(0,1);
    nod[1]=plane(1,1);
    nod[2]=plane(2,1);
    Node[1].SetCoord(nod);
    
    nod[0]=plane(0,2);
    nod[1]=plane(1,2);
    nod[2]=plane(2,2);
    Node[2].SetCoord(nod);
    
    nod[0]=plane(0,3);
    nod[1]=plane(1,3);
    nod[2]=plane(2,3);
    Node[3].SetCoord(nod);
    
    int nNods= gmesh->NNodes();
    
    gmesh->NodeVec().Resize(nNods+4);
    gmesh->NodeVec()[nNods]=Node[0];
    gmesh->NodeVec()[nNods+1]=Node[1];
    gmesh->NodeVec()[nNods+2]=Node[2];
    gmesh->NodeVec()[nNods+3]=Node[3];
    
    TPZVec<int64_t> cords(4);
    cords[0]=nNods;
    cords[1]=nNods+1;
    cords[2]=nNods+2;
    cords[3]=nNods+3;
    

    TRSFracPlane fracplane(plane);
    TRSRibFrac RibV(fracplane,gmesh);
    
    int64_t Nels = gmesh->NElements();
    RibV.CreateSkeletonElements(2, 4);
    RibV.CreateSkeletonElements(1, 4);

    gmesh->CreateGeoElement(EQuadrilateral, cords, 40, Nels);
    gmesh->BuildConnectivity();
    RibV.CreateSkeletonElements(1, 40);
    
    // std::ofstream out3("3DMESH.vtk");
    // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out3, true);




    //search gmesh for cut ribs
    for(int i = 0; i< Nels; i++){
      TPZGeoEl *gel = gmesh->Element(i);
      int nSides = gel->NSides();
      //skip all elements that aren't ribs
      if(gel->Dimension()!=1){continue;}
      for (int side=0; side< nSides; side++){
        if(gel->NSideNodes(side)==2){
          int64_t p1 =  gel->SideNodeIndex(side, 0);
          int64_t p2 =  gel->SideNodeIndex(side, 1);
          TPZVec<REAL> pp1(3);
          TPZVec<REAL> pp2(3);
          gmesh->NodeVec()[p1].GetCoordinates(pp1);
          gmesh->NodeVec()[p2].GetCoordinates(pp2);
          bool resul = RibV.Check_rib(pp1, pp2);
          if(resul==true){
              TRSRibs rib(i, true);
              RibV.AddRib(rib);
              TPZVec<REAL> ipoint =RibV.CalculateIntersection(fracplane, pp1, pp2);
              //5O is the material of children ribs
              RibV.Rib(i)->DivideRib(gmesh, ipoint, 50);
              gel->SetMaterialId(12); //when will we delete these ribs? //gmesh->DeleteElement(gel,i);
              
              // std::cout<<"Element: "<<i<<" Side: "<<side<<" Rib status: "<<resul<<" Fracture : 1"<<"\n";
          }
        // bool resul2 = RibV2.Check_rib(pp1, pp2);
        // if(resul2==true){
        //  gel->SetMaterialId(3);
        //  std::cout<<"Element: "<<i<<" Side: "<<side<<" Rib status: "<<resul<<" Fracture : 2"<<"\n";
        // }
        }
      }
    }
    
    // std::ofstream out2("./TestRibs.vtk");
    // TPZVTKGeoMesh::PrintGMeshVTK(RibV.GetgeoMesh(), out2, true);
    
    RibV.CreateSurfaces(20);  
    std::ofstream out("./TestSurfaces.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(RibV.GetgeoMesh(), out, true);
    
    return 0;
}
















  //  RibV.HasLowerDimensionNeighbour(geoside);
//
//
////    Setting a second plane
//    plane(0,0)=22;
//    plane(0,1)=20;
//    plane(0,2)=-10;
//    plane(1,0)=22;
//    plane(1,1)=20;
//    plane(1,2)=10;
//    plane(2,0)=2;
//    plane(2,1)=2;
//    plane(2,2)=10;
//    plane(3,0)=2;
//    plane(3,1)=2;
//    plane(3,2)=-10;
//
//    TPZVec<TPZGeoNode> Node2(4);
//    TPZVec<REAL> nod2(3,0);
//    nod2[0]=plane(0,0);
//    nod2[1]=plane(0,1);
//    nod2[2]=plane(0,2);
//    Node2[0].SetCoord(nod2);
//
//    nod2[0]=plane(1,0);
//    nod2[1]=plane(1,1);
//    nod2[2]=plane(1,2);
//    Node2[1].SetCoord(nod2);
//
//    nod2[0]=plane(2,0);
//    nod2[1]=plane(2,1);
//    nod2[2]=plane(2,2);
//    Node2[2].SetCoord(nod2);
//
//    nod2[0]=plane(3,0);
//    nod2[1]=plane(3,1);
//    nod2[2]=plane(3,2);
//    Node2[3].SetCoord(nod2);
//    int nNods2= gmesh->NNodes();
//
//    gmesh->NodeVec().Resize(nNods2+4);
//    gmesh->NodeVec()[nNods2]=Node2[0];
//    gmesh->NodeVec()[nNods2+1]=Node2[1];
//    gmesh->NodeVec()[nNods2+2]=Node2[2];
//    gmesh->NodeVec()[nNods2+3]=Node2[3];
//
//    TPZVec<int64_t> cords2(4);
//    cords2[0]=nNods2;
//    cords2[1]=nNods2+1;
//    cords2[2]=nNods2+2;
//    cords2[3]=nNods2+3;
//
//    int64_t Nels2 = gmesh->NElements();
//    gmesh->CreateGeoElement(EQuadrilateral, cords2, 4, Nels2);
//
//
//    TRSRibFrac RibV2(plane);
   
 //    Setting coplanar Tolerance
//    RibV.SetTolerance(0.0001);
////    RibV2.SetTolerance(0.0001);
//
//    for(int i=0; i< Nels; i++){
//         TPZGeoEl *gel = gmesh->Element(i);
//        int nSides = gel->NSides();
//         if(gel->Dimension()!=1){continue;}
//        for (int side=0; side< nSides; side++){
//            if(gel->NSideNodes(side)==2){
//               int64_t p1 =  gel->SideNodeIndex(side, 0);
//               int64_t p2 =  gel->SideNodeIndex(side, 1);
//                TPZVec<REAL> pp1(3);
//                TPZVec<REAL> pp2(3);
//                gmesh->NodeVec()[p1].GetCoordinates(pp1);
//                gmesh->NodeVec()[p2].GetCoordinates(pp2);
//
//
//
//                bool resul = RibV.Check_rib(pp1, pp2);
//                if(resul==true){
//
//
//                    gel->SetMaterialId(2);
//                std::cout<<"Element: "<<i<<" Side: "<<side<<" Rib status: "<<resul<<" Fracture : 1"<<"\n";
//                }
////                bool resul2 = RibV2.Check_rib(pp1, pp2);
////                if(resul2==true){
////                      gel->SetMaterialId(3);
////                    std::cout<<"Element: "<<i<<" Side: "<<side<<" Rib status: "<<resul<<" Fracture : 2"<<"\n";
////                }
//            }
//        }
//    }
//    std::ofstream out("TestSkeleton2.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(TEST.GetgeoMesh(), out, true);
    
   // TPZVTKGeoMesh::PrintGMeshVTK(TEST.GetgeoMesh(), out2, true);
    
    
//    std::ofstream out("TestJorge.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
 
    
//    RibV.GetPlane().Print(std::cout);
    
//    bool verific = RibV.Check_rib(pcom, pcom2);
//    std::cout<<verific<<"\n";
//    RibV.fdata.Print(std::cout);
    
//    int nel = gmesh->NElements();
//    for(int i=0; i<nel; i++){
//        TPZGeoEl *gel = gmesh->Element(i);
//        if(gel->Dimension()==2){continue;}
//
//        TPZFMatrix<REAL> coordinates;
//        gel->NodesCoordinates(coordinates);
//        TPZVec<double> p1(3);
//        p1[0]=coordinates(0,0);
//        p1[1]=coordinates(1,0);
//        p1[2]=coordinates(2,0);
//
//        TPZVec<double> p2(3);
//        p2[0]=coordinates(0,1);
//        p2[1]=coordinates(1,1);
//        p2[2]=coordinates(2,1);
//
//        TEST.CalculateIntersection(p1, p2);
//
//    }