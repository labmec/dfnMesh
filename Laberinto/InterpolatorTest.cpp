//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeopyramid.h"
#include "TPZGeoLinear.h"

#include <time.h>
#include <stdio.h>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include "TPZVTKGeoMesh.h"
#include <set>
#include <map>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/plot.hpp>
#include "TRSLinearInterpolator.h"
#include "pzgengrid.h"

using namespace std;
using namespace cv;

int main()
{
   
    // Creating the Geo mesh
    TPZManVector<REAL,3> x0(3,0.),x1(3,2);
    x1[2] = 0.;
    TPZManVector<int,2> nelx(2,4);
    nelx[0] = 2;
    TPZGenGrid gengrid(nelx,x0,x1);
    gengrid.SetElementType(ETriangle);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gengrid.Read(gmesh);
    
    //gengrid.SetBC(TPZGeoMesh *gr, int side, int bc)
    gengrid.SetBC(gmesh, 4, -1);
    gengrid.SetBC(gmesh, 5, -2);
    gengrid.SetBC(gmesh, 6, -3);
    gengrid.SetBC(gmesh, 7, -4);
    
    
    std::ofstream out("DamnerMalla.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    
    
    
    //gengrid.Read(gmesh,2);
    
    
    
    
    
    
    
//TRSLinearInterpolator Test1;
//Test1.ReadData("Bg.txt");
//Test1.GetData().Print(std::cout);
//   //Test1.GetData().Print(std::cout);
//Test1.SetInterpolationType(TRSLinearInterpolator::InterpType::TLinear);
//Test1.SetLeftExtension(TRSLinearInterpolator::Extension::ESlope,0.001);
//Test1.SetRightExtension(TRSLinearInterpolator::Extension::ELinear);
//auto BgFunction(Test1.GetFunction());
//
//std::cout<<BgFunction(50000.0)<<"\n";
   
return 0;
    
}
