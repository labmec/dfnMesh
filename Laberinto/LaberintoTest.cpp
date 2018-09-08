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


using namespace std;
using namespace cv;

int main()
{
//    Matrix data(2,2);
//    data(0,0)=1;
//    data(1,0)=2;
//    data(0,1)=1;
//    data(1,1)=4;
//    TRSLinearInterpolator jos;
//    jos.SetData(data);
//    jos.SetLeftExtension(TRSLinearInterpolator::Slope,5);
//    std::cout<<jos.Val(0);
//    
//    auto val(jos.GetFunction());
//    std::function<double(double)> Rs(jos.GetFunction());
//    auto valDeriv(jos.GetFunctionDeriv());
//    std::cout<<Rs(40000)<<"\n";
//    std::cout<<val(10);
    
    
    Mat image = imread("Small.png",IMREAD_GRAYSCALE);
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
    
    
        //Creando la malla geometrica
        TPZManVector<REAL,3> x0(3,0.),x1(3,px);
        x1[2] = 0.;
        TPZManVector<int,2> nelx(2,py);
        nelx[0] = px;
        TPZGenGrid gengrid(nelx,x0,x1);
        gengrid.SetElementType(EQuadrilateral);
        TPZGeoMesh *gmesh = new TPZGeoMesh;
        gmesh->SetDimension(2);
        gengrid.Read(gmesh);
        //MatsID
        int nels = gmesh->NElements();
    
        for (int i=0; i<nels; i++) {
          
            gmesh->Element(i)->SetMaterialId(vec[i]);
            
        }
    
    
        {
            std::ofstream out("LaberintoTestt.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
        }
    gmesh->BuildConnectivity();
    
    TPZCompMesh cmesh(gmesh);
    
    
    
    
    return 0;
}



