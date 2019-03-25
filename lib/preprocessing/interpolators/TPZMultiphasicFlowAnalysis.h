
#ifndef TPZMultiphasicDarcyAnalysisH
#define TPZMultiphasicDarcyAnalysisH

#include "pznonlinanalysis.h"
#include "pzcondensedcompel.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzstepsolver.h"
#include "pzl2projection.h"
#include "pzgradientreconstruction.h"

// #include "pzbfilestream.h"


class TPZCompMesh;

/** @author Jose Villegas in 07/03/2019
 * @brief class which implements 2D analysis for multiphasic axisimetric darcy flow
 */
class TPZMultiphasicFlowAnalysis : public TPZNonLinearAnalysis
{
private:
    TPZCompMesh *fcmesh;
    
public:
 
    void SetCmesh(TPZCompMesh *cmesh);
    TPZCompMesh * GetCmesh();
    void Assemble( TPZAnalysis *an);
 
    void LoadSolution();
    void PostProc( TPZStack<std::string> scalar_names, TPZStack<std::string> vec_names);
   
};

#endif
