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

#include <set>
#include <map>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/plot.hpp>
#include "TRSLinearInterpolator.h"


using namespace std;
using namespace cv;

double f(double x);
int main()
{
   
    int nPoints = 1000;
  
    Matrix data(nPoints,2);
    double h = 2*M_PI/nPoints;
    for (int i = 0; i<nPoints; i++) {
        data(i,0)= i*h;
        data(i,1)= f(i*h);
    }
   
    
    TRSLinearInterpolator CosInter(data);
    auto CosInterpolator(CosInter.GetFunction());
    std::cout<<CosInterpolator(4.7)<<"\n";
    std::cout<<f(4.7)<<"\n";
    for (int i = 0; i<nPoints; i++) {
        xData.at<double>(i) = i*h;
        xData.at<double>(i) = CosInterpolator(i*h);
        
    }

    
    
 
    return 0;
}

double f(double x){
    double y = cos(x);
    return y;
}



