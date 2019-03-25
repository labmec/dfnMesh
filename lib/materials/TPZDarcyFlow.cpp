/**
 * @file
 * @brief Contains implementations of the TPZDarcyFlow methods.
 */

#include "TPZDarcyFlow.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzmaterialdata.h"
#include <math.h>
#include "pzlog.h"

#include <cmath>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.poisson3d"));
#endif

using namespace std;
/**
 * Empty Constructor
 */
TPZDarcyFlow::TPZDarcyFlow(){
    
}

/** Creates a material object and inserts it in the vector of
 *  material pointers of the mesh.
 */
TPZDarcyFlow::TPZDarcyFlow(int matid): TPZDiscontinuousGalerkin(matid){
    
}


/** Creates a material object based on the referred object and
 *  inserts it in the vector of material pointers of the mesh.
 */
TPZDarcyFlow::TPZDarcyFlow(const TPZDarcyFlow &mat){
    
}

/**
 * Destructor
 */
TPZDarcyFlow::~TPZDarcyFlow(){
    
}

/** Fill material data parameter with necessary requirements for the
 * Contribute method. Here, in base class, all requirements are considered
 * as necessary. Each derived class may optimize performance by selecting
 * only the necessary data.
 * @since April 10, 2007
 */
void TPZDarcyFlow::FillDataRequirements(TPZVec<TPZMaterialData> &datavec){
    
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(true);
        datavec[idata].fNeedsSol = true;
    }
    
    
}

void TPZDarcyFlow::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec){
    
}



/** print out the data associated with the material */
void TPZDarcyFlow::Print(std::ostream &out ){
    
}

/** returns the variable index associated with the name */
int TPZDarcyFlow::VariableIndex(const std::string &name){
    if (!strcmp("Solution", name.c_str())) return 0;
    if (!strcmp("Derivative", name.c_str())) return 1;
}

/** returns the number of variables associated with the variable
 indexed by var.  var is obtained by calling VariableIndex */
int TPZDarcyFlow::NSolutionVariables(int var){
    switch (var) {
        case 0:
            return 1; // Scalar
            break;
        case 1:
            return 3; // Vector
            break;
    }

}

/** Computes the divergence over the parametric space */
void TPZDarcyFlow::ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi){
    
}

/** returns the solution associated with the var index based on
 * the finite element approximation */
void TPZDarcyFlow::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0: //Presure
        {
            Solout = datavec[0].sol[0];
            break;
        }
        case 1:
        {
            Solout[0]=datavec[0].dphix[0];
            Solout[0]=datavec[1].dphix[1];
            Solout[0]=datavec[2].dphix[2];
            
        }
            
    }
}

// Contribute Methods being used

/**
 * It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TPZDarcyFlow::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
//    data.SetAllRequirements(true);
//    data.fNeedsSol=true;
    
    TPZAutoPointer<TPZFunction<STATE>> fo = ForcingFunction();
    TPZFMatrix<REAL> phiuH1 = data.phi;
   
    TPZFMatrix<REAL> dphiuH1 = data.dphix;
    TPZVec<REAL> point = data.x;
    TPZVec<double> f(3,0.0);
    int nphi = data.dphix.Cols();
   
    for (int i=0; i<nphi; i++){
        for (int j=0; j<nphi; j++) {
            ek(i,j)+= weight*(dphiuH1(0,i)*dphiuH1(0,j) + dphiuH1(1,i)*dphiuH1(1,j)) ;
        }
        fo->Execute(point,f);
        ef(i,0) += (-1.0)*weight*f[0]*phiuH1(i,0);
      
    }
    
}



 void TPZDarcyFlow::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){

     TPZMaterial *material = bc.Material();
     int mat_ID = bc.Id();
     TPZFMatrix<REAL> dphiuH1 = data.dphi;
     TPZFMatrix<REAL> phiuH1 = data.phi;
     
     REAL value = bc.Val2()(0,0);
  //   double un = data.sol[0][0];
     int nphi = data.dphix.Cols();
     
//     for (int iq = 0; iq < nPhiHdiv; iq++)
//     {
//         ef(iq) += weight * (gBigNumber * (un - Value)) * phiuH1(iq,0);
//
//         for (int jq = 0; jq < nPhiHdiv; jq++)
//         {
//             ek(iq,jq) += gBigNumber * weight * (phiuH1(jq,0)) * phiuH1(iq,0);
//         }
//
//     }
//
    

     int type = bc.Type();
     switch (type) {
         case 0: //Pressure
         {
             for (int i=0; i<nphi; i++){
                 for (int j=0; j<nphi; j++) {
                     ek(i,j)+= weight * gBigNumber * phiuH1(i,0)*phiuH1(j,0); ;
                     ek.Print(std::cout);
                 }
                 
                 ef(i,0) +=gBigNumber*weight*value*phiuH1(i,0);
                 ef.Print(std::cout);
                 std::cout<<std::endl;
             }
             break;
         }
         case 1: //Flux
         {
             for (int i=0; i<nphi; i++){
                 ef(i,0) += weight*value*phiuH1(i,0);
                 ef.Print(std::cout);
                 std::cout<<std::endl;
             }
             break;
         }
        
     }
     
     
}


void TPZDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    ContributeDarcy(datavec, weight,ek, ef);
    std::cout<<"HERE";
    ek.Print(std::cout);
    
}

/**
 * It computes a contribution to the load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TPZDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    std::cout<<"Pare"<<std::endl;
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
}


/**
 * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
    
    
}


//   Contribute for Darcy system

/**
 * It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
     TPZFMatrix<REAL> phiuH1 = datavec[0].phi;
     TPZFMatrix<REAL> dphiuH1 = datavec[0].dphix;
     int nphi = datavec[0].fVecShapeIndex.size();
    
    
}

/**
 * It computes a contribution to the load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeInterfaceDarcy(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeInterfaceDarcy(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
    
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeBCDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
}


//   Contribute for Transport of alpha phase

/**
 * It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeAlpha(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    
}

/**
 * It computes a contribution to the load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeAlpha(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeBCInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeBCInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
}



/**
 * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZDarcyFlow::ContributeInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
    
    
}


/**
 * Save the element data to a stream
 */
void TPZDarcyFlow::Write(TPZStream &buf, int withclassid) const{
    
}

/**
 * Read the element data from a stream
 */
void TPZDarcyFlow::Read(TPZStream &buf, void *context){
    
    
}
