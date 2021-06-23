/*! 
 *	DFNRib.cpp
 *  @authors   Pedro Lima
 *  @authors   Jorge Ordoñez
 *  @date      2018-2019
 */

#include "DFNRib.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "TPZRefPatternDataBase.h"
#include "DFNFracture.h"
#include <math.h>

#if PZ_LOG
    static TPZLogger logger("dfn.mesh");
#endif // PZ_LOG

/// Default constructor takes a pointer to geometric element
DFNRib::DFNRib(TPZGeoEl *gel, DFNFracture *Fracture) :
    fGeoEl(gel),
    fIntersectionSide(2),
    fFracture(Fracture),
    fIntersectionIndex(-1),
    fOffbound(false)
{fCoord.resize(0);};

// Copy constructor
DFNRib::DFNRib(const DFNRib &copy){
    this->operator=(copy);
}

// Assignment operator
DFNRib &DFNRib::operator=(const DFNRib &copy){
    fGeoEl = copy.fGeoEl;
    fIntersectionSide = copy.fIntersectionSide;
    fCoord = copy.fCoord;
    fIntersectionIndex = copy.fIntersectionIndex;
    fFracture = copy.fFracture;
    fOffbound = copy.fOffbound;
    return *this;
}

inline REAL VectorNorm(TPZManVector<REAL,3> &vector);
inline REAL Distance(TPZManVector<REAL,3> &vector1, TPZManVector<REAL,3> &vector2);


/// Get real intersection coordinates
TPZManVector<REAL, 3> DFNRib::RealCoord() const{
    TPZManVector<REAL, 3> coord(3, 0);
    TPZGeoMesh *gmesh = fGeoEl->Mesh();

    if(fIntersectionSide == 2 && fIntersectionIndex < 0)
        {coord = this->fCoord;}
    else
        {gmesh->NodeVec()[fIntersectionIndex].GetCoordinates(coord);}
    
    return coord;
}


void DFNRib::Refine(){
    if(!this->NeedsRefinement()) return;
    if(!fGeoEl){
        std::cout<<"No gel associated to the Rib\n";
        DebugStop();
    }

    // set new node in gmesh if it hasn't been set yet
    TPZGeoMesh *gmesh = fGeoEl->Mesh();
    if(fIntersectionIndex < 0){
        fIntersectionIndex = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[fIntersectionIndex].Initialize(fCoord,*gmesh);
    }
    
    // Set refinement pattern
    // @TODO discuss with your advisor if this is really necessary
    // where is the refinement pattern used?
    this->CreateRefPattern();

    // set children
        TPZManVector<int64_t,2> cornerindices(2,0);
        TPZGeoEl *child;
        int64_t elindex = gmesh->NElements();
        for(int i=0; i<2; i++){
            cornerindices[i%2] = fGeoEl->NodeIndex(i%2);
            cornerindices[(i+1)%2] = fIntersectionIndex;
            child = gmesh->CreateGeoElement(EOned, cornerindices, fGeoEl->MaterialId(), elindex);
            fGeoEl->SetSubElement(i,child);
            child->SetFatherIndex(fGeoEl->Index());
            elindex++;
        }
}







void DFNRib::CreateRefPattern(){
        // refinement mesh
        TPZGeoMesh refPatternMesh;
        int refnnodes = fGeoEl->NCornerNodes() + this->NeedsRefinement();
        TPZGeoMesh *gmesh = fGeoEl->Mesh();

        // set nodes
        refPatternMesh.NodeVec().Resize(refnnodes);
        TPZManVector<REAL,3> coord(3);
        for(int i = 0; i<3; i++){
            if(i < 2) {gmesh->NodeVec()[fGeoEl->NodeIndex(i)].GetCoordinates(coord);}
            else      {gmesh->NodeVec()[fIntersectionIndex].GetCoordinates(coord);}

            refPatternMesh.NodeVec()[i].Initialize(coord,refPatternMesh);
        }
        
        // insert father
        TPZManVector<int64_t,2> cornerindices({0,1});
        int64_t elindex = 0;
        refPatternMesh.CreateGeoElement(EOned, cornerindices, fGeoEl->MaterialId(), elindex);
        
        // insert children
        cornerindices = {0,2};
        elindex = 1;
        refPatternMesh.CreateGeoElement(EOned, cornerindices, fGeoEl->MaterialId(), elindex);
        cornerindices = {2,1};
        elindex = 2;
        refPatternMesh.CreateGeoElement(EOned, cornerindices, fGeoEl->MaterialId(), elindex);
                    
        // define refPattern
        refPatternMesh.BuildConnectivity(); 
        TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(refPatternMesh);
        fGeoEl->SetRefPattern(refpat);
    }





void DFNRib::SnapIntersection_force(int closestnode){
    REAL BIG_NUMBER = 1.e12;
    // this->SnapIntersection_try(BIG_NUMBER);
    if(this->CanBeSnapped() == false) return;
    // If closest node was filled, snapped will be forced to node of local index == closestnode. If left -1, code will compute closest node.
    switch (closestnode){
        case -1: NeedsSnap(closestnode,BIG_NUMBER); break; // If user gave -1, they don't know which node is closest, and need this method to compute it
        case  0:
        case  1: break;
        default:{
            std::stringstream error; 
	        error << "An edge doesn't have a node of local index = " << closestnode << "\n  where():  " << __PRETTY_FUNCTION__ << '\n';
            throw std::invalid_argument(error.str()); break;
        } 
    }

    // Snap intersection
    fIntersectionIndex = fGeoEl->NodeIndex(closestnode);
    fIntersectionSide = closestnode;
#if PZ_LOG
    if(logger.isDebugEnabled())
        LOGPZ_DEBUG(logger,"Snap rib " << Index() << ". Towards LocalNode " << closestnode << "  - GlobalNode " << fGeoEl->NodeIndex(closestnode));
#endif // PZ_LOG
    // When an intersection gets snapped, neighbours should also be snapped to keep refinements consistant
    UpdateNeighbours(closestnode);
}



bool DFNRib::SnapIntersection_try(REAL tolDist){
    // If there's an intersection at the 1D side of this DFNRib, check if 
    // that intersection violates the tolerable distance and, if necessary, snap it down to a lower-dimensional side
    int closestnode = -1;
    // NeedsSnap will modify closestnode to zero or one if Snapping is required
    if(this->NeedsSnap(closestnode,tolDist)){
        SnapIntersection_force(closestnode);
        return true;
    }
    return false;
}

/// @todo extract a CanBeSnapped function from this to rewrite SnapIntersection_force and SnapIntersection_try
bool DFNRib::NeedsSnap(int& closestnode, REAL tolDist){
    if(this->CanBeSnapped() == false){
        closestnode = fGeoEl->NodeIndex(fIntersectionSide);
        return false;
    }
    
    TPZManVector<REAL,3> coord(3,0);
    fGeoEl->NodePtr(0)->GetCoordinates(coord);
    REAL dist0 = Distance(fCoord,coord);

    fGeoEl->NodePtr(1)->GetCoordinates(coord);
    REAL dist1 = Distance(fCoord,coord);

    closestnode = (dist0 < dist1 ? 0 : 1);
    REAL dist = MIN(dist0,dist1);

    return (dist < tolDist);
}



void DFNRib::UpdateNeighbours(int iside){

    TPZGeoElSide gelside(fGeoEl,iside);
    TPZGeoElSide neighbour = gelside.Neighbour();
    TPZGeoEl *gel;
    for(/*void*/; neighbour != gelside; neighbour = neighbour.Neighbour()){
        gel = neighbour.Element();
        if(!gel) break;
        if(gel->HasSubElement()) continue;
        int neig_side = neighbour.Side();

        switch(gel->Dimension()){
            case 0: {continue;}
            case 1:{
                // check if DFNRib exists
                DFNRib *rib_ptr = fFracture->Rib(gel->Index());
                if(!rib_ptr){
                    break;
                    // DFNRib neig_rib(gel,fFracture);
                    // fFracture->AddRib(neig_rib);
                    // rib_ptr = fFracture->Rib(gel->Index());
                }
                
                // A rib needs refinement when its intersection node wasn't snapped to one of the nodes
                if(rib_ptr->NeedsRefinement()) {rib_ptr->SnapIntersection_force();}
                break;
            }
            case 2:{
                // check if DFNFace exists
                DFNFace *neig_face = fFracture->Face(gel->Index());
                if(!neig_face) break;
                if(!neig_face->UpdateStatusVec()) break;
		        neig_face->UpdateRefMesh();
                break;
            }
            case 3: {continue;}
        }
    }
    return;
}

/** @brief Adds a pointer of this rib into the corresponding position of its neighbour faces ribvectors*/
void DFNRib::AppendToNeighbourFaces(){
    // Find a neighbour DFNFace through side 2, and change its fRibs vector.
    TPZGeoElSide gelside(fGeoEl,2);
    TPZGeoElSide neig;
    for(neig=gelside.Neighbour(); neig != gelside; neig = neig.Neighbour()){
        if(neig.Element()->Dimension() != 2) continue;
        if(neig.Element()->HasSubElement()) continue;
        // if(neig.Element()->MaterialId() != DFNMaterial::Efracture && 
        //    neig.Element()->MaterialId() != DFNMaterial::Eskeleton) continue;
        int64_t neigindex = neig.Element()->Index();
        DFNFace* dfnface = fFracture->Face(neigindex);
        // If DFNFace is missing, create one
        if(!dfnface){
            // TPZManVector<DFNRib*,4> rib_vec(neig.Element()->NSides(1),nullptr);
            // rib_vec[neig.Side()-neig.Element()->NSides(0)] = this;
            // DFNFace new_dfnface(neig.Element(),fFracture,rib_vec);
            DFNFace new_dfnface(neig.Element(),fFracture);
            dfnface = fFracture->AddFace(new_dfnface);
        }
        dfnface->AddRib(this,neig.Side());
        dfnface->UpdateStatusVec();
        dfnface->UpdateRefMesh();
    }
}



// // Might be useful at some point... but not yet
// bool DFNRib::UpdateMaterial(){
//     DebugStop();
//     if(fGeoEl->MaterialId() == DFNMaterial::Efracture) {return false;}
//     if(fIntersectionSide[0] && fIntersectionSide[1]){
//         fGeoEl->SetMaterialId(DFNMaterial::Efracture);
//         return true;
//     }
//     return false;
// }




inline REAL VectorNorm(TPZManVector<REAL,3> &vector){
    int n = vector.NElements();
    REAL temp = 0.0;
    for (int i = 0; i < n; i++){
        temp += vector[i]*vector[i];
    }
    return sqrt(temp);
}

inline REAL Distance(TPZManVector<REAL,3> &vector1, TPZManVector<REAL,3> &vector2){
    int n = vector1.NElements();
    TPZManVector<REAL,3> difference(n,0);
    for(int i = 0; i<n;i++){
        difference[i] = vector1[i] - vector2[i];
    }
    return VectorNorm(difference);
}



void DFNRib::Print(std::ostream &out) const
{
    out<<"\nRib GeoEl index # "<<fGeoEl->Index();
	out<<"\nSide intersected : "<<fIntersectionSide;
	out<<"\nIntersection Coord : "<< fCoord;
	out<<"\nIntersection Node index: " << fIntersectionIndex;
    out<<"\nRib Nodes:\n";
    out<<"\t0: index: "<<fGeoEl->NodeIndex(0)<<" "; fGeoEl->Node(0).Print(out);
    out<<"\t1: index: "<<fGeoEl->NodeIndex(1)<<" "; fGeoEl->Node(1).Print(out);
    out<<"Subelements:";
    if(fGeoEl->HasSubElement()){
        out<<"  "<<fGeoEl->SubElement(0)->Index();
        out<<"  "<<fGeoEl->SubElement(1)->Index();
    }else{ out << "xxx";}
    out<<"\n\n";
}
