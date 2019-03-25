//  TPZMixedDarcyFlow_H
//  PZ
//
//  Created by José Villegas on 03/14/19.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef TPZMixedDarcyFlow_H
#define TPZMixedDarcyFlow_H

#include "TPZMaterial.h"
#include "pzdiscgal.h"
#include "pzpoisson3d.h"
#include "TPZDarcyFlow.h"
#include "TPZMaterial.h"
#include "pzfunction.h"

/**
 * @ingroup material
 * @author Agnaldo Farias
 * @since 5/28/2012
 * @brief Material to solve a mixed poisson problem 2d by multiphysics simulation
 * @brief Pressure(p): uses L2 space.  Velocity (Q): uses Hdiv space
 */
/**
 * \f$ Q = -(k/visc)*grad(p)  ==> Int{Q.q}dx - (k/visc)*Int{p*div(q)}dx + (k/visc)*Int{pD(q.n)}ds = 0  (Eq. 1)  \f$
 *
 * \f$ div(Q) = f  ==> Int{div(Q)*v}dx = Int{f*v}dx (Eq. 2) \f$
 *
 * \f$ p = pD in Dirichlet boundary and Q.n = qN in Neumann boundary\f$
 */


class TPZMixedDarcyFlow : public TPZDarcyFlow {
    
protected:
    REAL ff;
    int fDim =2;
    
public:
    TPZMixedDarcyFlow();
    
    TPZMixedDarcyFlow(int matid, int dim);
    
    virtual ~TPZMixedDarcyFlow();
    
    TPZMixedDarcyFlow(const TPZMixedDarcyFlow &cp);
    
    TPZMixedDarcyFlow &operator=(const TPZMixedDarcyFlow &copy);
    
    virtual TPZMaterial * NewMaterial(){
        return new TPZMixedDarcyFlow(*this);
    }
    
    
    virtual void Print(std::ostream & out);
    
    virtual std::string Name() { return "TPZMixedPoisson"; }
    
    virtual int NStateVariables();
    
    void SetPermeability(REAL perm) {
      
    }
    
    //Set the permeability tensor and inverser tensor
    void SetPermeabilityTensor(TPZFMatrix<REAL> K, TPZFMatrix<REAL> invK){
        
       
    }
    
    void SetViscosity(REAL visc) {
       
    }
    
    void GetPermeability(REAL &perm) {
      
    }
    
    void SetInternalFlux(REAL flux) {
     
    }
    
    void SetStabilizedMethod(){
      
    }
    
    void SetHdois(){
       
    }
    
    void SetStabilizationCoeficients(REAL delta1, REAL delta2){
        
    }
    
    void SetPermeabilityFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
       
    }
    
    TPZAutoPointer<TPZFunction<STATE> > PermeabilityFunction()
    {
       
    }
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    
    virtual int VariableIndex(const std::string &name);
    
    virtual int NSolutionVariables(int var);
    
    /**
     * @brief It return a solution to multiphysics simulation.
     * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
    {
        // default is no specific data requirements
        int nref = datavec.size();
        for (int iref = 0; iref <nref; iref++) {
            datavec[iref].SetAllRequirements(false);
        }
        if(type == 50)
        {
            for(int iref = 0; iref<nref; iref++){
                datavec[iref].fNeedsSol = true;
            }
        }
    }
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    
    virtual int NEvalErrors() {return 3;}
    
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors);
    
public:
    virtual int ClassId() const;
};

#endif
