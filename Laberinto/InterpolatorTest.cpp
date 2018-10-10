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

int main()
{
   
    TRSLinearInterpolator Test1;
    Test1.ReadData("Bg.txt");
    Test1.GetData().Print(std::cout);
   // Test1.GetData().Print(std::cout);
Test1.SetInterpolationType(TRSLinearInterpolator::InterpType::TLinear);
    Test1.SetLeftExtension(TRSLinearInterpolator::Extension::Slope,0.001);
    Test1.SetRightExtension(TRSLinearInterpolator::Extension::Linear);
    auto BgFunction(Test1.GetFunction());
    std::cout<<BgFunction(50000.0)<<"\n";
   
    return 0;
}


