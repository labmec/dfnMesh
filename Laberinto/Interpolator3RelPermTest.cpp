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

int main(){
    TRSLinearInterpolator krw;
    TRSLinearInterpolator krow;
    TRSLinearInterpolator krg;
    TRSLinearInterpolator krog;
    krw.ReadData("krw.txt");
    krow.ReadData("krow.txt");
    krg.ReadData("krg.txt");
    krg.ReadData("krog.txt");
    
    
    Interpol3RelPerm Kro;
    Kro.SetData(krw.GetFunctionDeriv(), Interpol3RelPerm::Kralpha::EKrw);
    Kro.SetData(krow.GetFunctionDeriv(), Interpol3RelPerm::Kralpha::EKrow);
    Kro.SetData(krg.GetFunctionDeriv(), Interpol3RelPerm::Kralpha::EKrg);
    Kro.SetData(krog.GetFunctionDeriv(), Interpol3RelPerm::Kralpha::EKrog);
    
    auto KroFunc(Kro.GetFunctionDeriv());
    Kro.SetData(KroFunc, Interpol3RelPerm::Kralpha::EKrw);
    Kro.Val(0.4, 0.3);
}
