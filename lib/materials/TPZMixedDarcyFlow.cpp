/**
 * @file
 * @brief Contains the methods of the TPZMixedPoisson class (multiphysics environment)
 * @author Agnaldo Farias
 * @date 2012/05/28
 */

//#include "mixedpoisson.h"
#include "TPZDarcyFlow.h"
#include "TPZMixedDarcyFlow.h"
#include "pzlog.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"


#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedpoisson.data"));
static LoggerPtr logerror(Logger::getLogger("pz.mixedpoisson.error"));
#endif

TPZMixedDarcyFlow::TPZMixedDarcyFlow(): TPZRegisterClassId(&TPZMixedDarcyFlow::ClassId), TPZDarcyFlow() {
  
}

TPZMixedDarcyFlow::TPZMixedDarcyFlow(int matid, int dim): TPZRegisterClassId(&TPZMixedDarcyFlow::ClassId), TPZDarcyFlow(matid) {
 
    
}

TPZMixedDarcyFlow::~TPZMixedDarcyFlow() {
}

TPZMixedDarcyFlow::TPZMixedDarcyFlow(const TPZMixedDarcyFlow &cp) :TPZRegisterClassId(&TPZMixedDarcyFlow::ClassId), TPZDarcyFlow(cp) {
    
}

TPZMixedDarcyFlow & TPZMixedDarcyFlow::operator=(const TPZMixedDarcyFlow &copy){
   
}

int TPZMixedDarcyFlow::NStateVariables() {
    return 1;
}

void TPZMixedDarcyFlow::Print(std::ostream &out) {
    out << "name of material : " << Name() << "\n";
    out << "Base Class properties :";
    Print(out);
    out << "\n";
}

void TPZMixedDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    
#ifdef PZDEBUG
    int nref = datavec.size();
    if (nref != 2) {
        std::cout << " Error. The size of the datavec is different from 2." << std::endl;
        DebugStop();
    }
#endif
    
    // Setting the phis
    TPZFMatrix<REAL> phiQ = datavec[0].phi;
    TPZFMatrix<REAL> phiP = datavec[1].phi;
    TPZFMatrix<REAL> dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> dphiP = datavec[1].dphix;
    
    phiQ.Print(std::cout);
    //Contribution: Matrix A
    int nphiQ  = datavec[0].fVecShapeIndex.NElements();
    int nphiP = phiP.Rows();
    
   
    for(int i=0; i<nphiQ; i++){
        int ivecind = datavec[0].fVecShapeIndex[i].first;
        int ishapeind = datavec[0].fVecShapeIndex[i].second;
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fNormalVec(id,ivecind);
            //ivec(1,0) = datavec[0].fNormalVec(1,ivecind);
            //ivec(2,0) = datavec[0].fNormalVec(2,ivecind);
        }
        
        for(int j=0; j<nphiQ; j++){
            int jvecind = datavec[0].fVecShapeIndex[j].first;
            int jshapeind = datavec[0].fVecShapeIndex[j].second;
            TPZFNMatrix<3,REAL> jvec(3,1,0.);
            for(int id=0; id<3; id++){
                jvec(id,0) = datavec[0].fNormalVec(id,jvecind);
                //ivec(1,0) = datavec[0].fNormalVec(1,ivecind);
                //ivec(2,0) = datavec[0].fNormalVec(2,ivecind);
            }
            REAL val = (ivec(0,0)*jvec(0,0))+(ivec(1,0)*jvec(1,0))+(ivec(2,0)*jvec(2,0));
            ek(i,j) += weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*val;
            
        }
    }

    
    // Coupling terms between flux and pressure. Matrix B
    for(int iq=0; iq<nphiQ; iq++)
    {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fNormalVec(id,ivecind);
            //ivec(1,0) = datavec[0].fNormalVec(1,ivecind);
            //ivec(2,0) = datavec[0].fNormalVec(2,ivecind);
        }
        TPZFNMatrix<3,REAL> axesvec(3,1,0.);
        datavec[0].axes.Multiply(ivec,axesvec);
        
        REAL divwq = 0.;
        for(int iloc=0; iloc<2; iloc++)
        {
            divwq += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
        }
        for (int jp=0; jp<nphiP; jp++) {
            
            REAL fact = (-1.)*weight*phiP(jp,0)*divwq;
            // Matrix B
            ek(iq, nphiQ+jp) += fact;
            
            // Matrix B^T
            ek(nphiQ+jp,iq) += fact;
        }
            
    }
    
    STATE force = ff;
    if(fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction->Execute(datavec[1].x,res);
        force = res[0];
    }
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<nphiP; ip++){
        ef(nphiQ+ip,0) += (-1.)*weight*force*phiP(ip,0);
    }
    
    
 
    
}



void TPZMixedDarcyFlow::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    

#ifdef PZDEBUG
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
        DebugStop();
    }
    if (bc.Type() > 2 ) {
        std::cout << " Erro.!! Neste material utiliza-se apenas condicoes de Neumann e Dirichlet\n";
        DebugStop();
    }
#endif
    
    TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
    int phrq = phiQ.Rows();
    
    REAL v2;
    if(bc.HasForcingFunction())
    {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9,STATE> gradu(Dimension(),1);
        bc.ForcingFunction()->Execute(datavec[0].x,res,gradu);
        v2 = res[0];
    }else
    {
        v2 = bc.Val2()(0,0);
    }
    
    switch (bc.Type()) {
        case 0 :        // Dirichlet condition
            //primeira equacao
            for(int iq=0; iq<phrq; iq++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq,0) += (-1.)*v2*phiQ(iq,0)*weight;
            }
            break;
            
        case 1 :            // Neumann condition
            //primeira equacao
            for(int iq=0; iq<phrq; iq++)
            {
                ef(iq,0)+= gBigNumber*v2*phiQ(iq,0)*weight;
                for (int jq=0; jq<phrq; jq++) {
                    
                    ek(iq,jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight;
                }
            }
            break;
            
        case 2 :            // mixed condition
            for(int iq = 0; iq < phrq; iq++) {
                
                ef(iq,0) += v2*phiQ(iq,0)*weight;
                for (int jq = 0; jq < phrq; jq++) {
                    ek(iq,jq) += weight*bc.Val1()(0,0)*phiQ(iq,0)*phiQ(jq,0);
                }
            }
            
            break;
    }
    
    
  
}

/** Returns the variable index associated with the name */
int TPZMixedDarcyFlow::VariableIndex(const std::string &name){
    if(!strcmp("Flux",name.c_str()))        return  31;
    if(!strcmp("Pressure",name.c_str()))    return  32;
    if(!strcmp("GradFluxX",name.c_str()))   return  33;
    if(!strcmp("GradFluxY",name.c_str()))   return  34;
    if(!strcmp("DivFlux",name.c_str()))   return  35;
    
    if(!strcmp("ExactPressure",name.c_str()))  return 36;
    if(!strcmp("ExactFlux",name.c_str()))  return 37;
    
    if(!strcmp("POrder",name.c_str()))        return  38;
    if(!strcmp("GradPressure",name.c_str()))        return  39;
    if(!strcmp("Divergence",name.c_str()))      return  40;
    if(!strcmp("ExactDiv",name.c_str()))        return  41;
    if (!strcmp("Derivative",name.c_str())) {
        return 42;
    }
    if (!strcmp("Permeability",name.c_str())) {
        return 43;
    }
    
    DebugStop();
    return -1;
}

int TPZMixedDarcyFlow::NSolutionVariables(int var){
    if(var == 31) return 3;
    if(var == 32 || var==8) return 1;
    if(var == 33) return 3;
    if(var == 34) return 3;
    if(var == 35) return 1;
    if(var == 36) return 1;
    if(var == 37) return 2;
    if(var == 38) return 1;
    if(var == 39) return 2;
    if(var == 40 || var == 41) return 1;
    if(var == 42) return 3;
    if(var == 43) return 1;
    return TPZMaterial::NSolutionVariables(var);
}

void TPZMixedDarcyFlow::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    Solout.Resize( this->NSolutionVariables(var));
    
    TPZVec<STATE> SolP, SolQ;
    

    
    // SolQ = datavec[0].sol[0];
    SolP = datavec[1].sol[0];
    
    if(var == 31){ //function (state variable Q)
        for (int i=0; i<3; i++)
        {
            Solout[i] = datavec[0].sol[0][i];
        }
        return;
    }
    
    if(var == 32){
        Solout[0] = SolP[0];//function (state variable p)
        return;
    }
    
    if(var==33){
        Solout[0]=datavec[0].dsol[0](0,0);
        Solout[1]=datavec[0].dsol[0](1,0);
        Solout[2]=datavec[0].dsol[0](2,0);
        return;
    }
    
    if(var==34){
        Solout[0]=datavec[0].dsol[0](0,1);
        Solout[1]=datavec[0].dsol[0](1,1);
        Solout[2]=datavec[0].dsol[0](2,1);
        return;
    }
    
    if(var==35){
        Solout[0]=datavec[0].dsol[0](0,0)+datavec[0].dsol[0](1,1);
        return;
    }
    
    TPZVec<REAL> ptx(3);
    TPZVec<STATE> solExata(1);
    TPZFMatrix<STATE> flux(fDim+1,1);
    
    //Exact soluion
    if(var == 36){
        fForcingFunctionExact->Execute(datavec[0].x, solExata,flux);
        Solout[0] = solExata[0];
        return;
    }//var6
    
    if(var == 37){
        fForcingFunctionExact->Execute(datavec[0].x, solExata,flux);
        Solout[0] = flux(0,0);
        Solout[1] = flux(1,0);
        return;
    }//var7
    
    if(var==38){
        Solout[0] = datavec[1].p;
        return;
    }
    
    if(var==39){
        TPZFNMatrix<3,REAL> dsoldx;
        TPZFMatrix<REAL> dsoldaxes(fDim,1);
        dsoldaxes(0,0) = datavec[1].dsol[0][0];
        dsoldaxes(1,0) = datavec[1].dsol[0][1];
        TPZAxesTools<REAL>::Axes2XYZ(dsoldaxes, dsoldx, datavec[1].axes);
        Solout[0] = dsoldx(0,0);
        Solout[1] = dsoldx(1,0);
        return;
    }
    
    if(var==40){
        Solout[0]=datavec[0].dsol[0](0,0)+datavec[0].dsol[0](1,1);
        return;
    }
    
    if(var==41){
        fForcingFunctionExact->Execute(datavec[0].x,solExata,flux);
        Solout[0]=flux(2,0);
        return;
    }
    if(var == 42)
    {
        for(int i=0; i<fDim; i++)
        {
            Solout[i] = 0.;
        }
        for (int i=0; i<fDim; i++) {
            for (int j=0; j<fDim; j++) {
          //      Solout[i] -= InvPermTensor(i,j)*datavec[0].sol[0][i];
            }
        }
        return;
    }
    if(var ==43)
    {
   //     Solout[0] = PermTensor(0,0);
        return;
    }
    TPZMaterial::Solution(datavec,var,Solout);
   
}


void TPZMixedDarcyFlow::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++ )
    {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsNeighborSol = false;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = false;
        datavec[i].fNeedsHSize = true;
    }
    
}


void TPZMixedDarcyFlow::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    
    
    //                             TPZVec<REAL> &x,TPZVec<STATE> &u,
    //                             TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &/*flux*/,
    //                             TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
    
   
}

int TPZMixedDarcyFlow::ClassId() const{
    return Hash("TPZMixedPoisson") ^ TPZMixedDarcyFlow::ClassId() << 1;
}

