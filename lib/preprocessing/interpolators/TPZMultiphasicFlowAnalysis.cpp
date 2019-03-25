#include "TPZMultiphasicFlowAnalysis.h"

void TPZMultiphasicFlowAnalysis::SetCmesh(TPZCompMesh *cmesh){
    fcmesh = cmesh;
}

TPZCompMesh * TPZMultiphasicFlowAnalysis::GetCmesh(){
    
    return fcmesh;
}

void TPZMultiphasicFlowAnalysis::Assemble(TPZAnalysis *an){
    
}

void TPZMultiphasicFlowAnalysis::LoadSolution(){
    
}

void TPZMultiphasicFlowAnalysis::PostProc( TPZStack<std::string> scalar_names, TPZStack<std::string> vec_names){
    
    int nscalar = scalar_names.size();
    int nvectors= vec_names.size();
    
}
