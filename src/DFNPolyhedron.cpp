#include "pzgmesh.h"
#include "DFNPolyhedron.h"
#include "DFNMesh.h"

/// Assignment operator
DFNPolyhedron &DFNPolyhedron::operator=(const DFNPolyhedron &copy){
    fShell.clear();
    fShell = copy.fShell;
    fIndex = copy.fIndex;
    fDFN = copy.fDFN;
    return *this;
}

/// Copy constructor
DFNPolyhedron::DFNPolyhedron(const DFNPolyhedron &copy){
    this->operator=(copy);
}


void DFNPolyhedron::SwapIndex(const int newid){
    TPZVec<std::array<int,2>>& polyh_per_face = fDFN->GetPolyhedraPerFace();
    fIndex = newid;
    for(auto& orientedface : fShell){
        int64_t faceid = orientedface.first;
        int orient = orientedface.second==1?0:1;
        polyh_per_face[faceid][orient] = newid;
    }
}


/** @brief Print method for logging */
void DFNPolyhedron::Print(std::ostream& out){
    int nelements = this->fDFN->Mesh()->NElements();
    int width = 2 + int(std::log10(nelements)+1);
    out << "\nPolyh#"<<std::setw(width-1)<<fIndex<<":";
    for(auto& oriented_face : fShell){
        // out << std::setw(5) << std::showpos << oriented_face << " ";
        out << std::setw(width) << oriented_face;
    }
    if(IsRefined()) out << "\t [Refined]";
}

void DFNPolyhedron::ListDFNFaces(DFNFracture* fracture, TPZStack<DFNFace*> facelist){
    facelist.clear();
    for(auto& orient_face : fShell){
        DFNFace* dfnface = fracture->Face(orient_face.first);
        if(!dfnface) continue;

        facelist.push_back(dfnface);
    }
}

/** @brief Remove faces from this polyhedron*/
void DFNPolyhedron::RemoveFaces(const TPZVec<std::pair<int64_t,int>>& facestack){
    for(auto& face_orient : facestack){
        int64_t index = face_orient.first;
        fShell.erase(index);
    }
}

/** @brief Checks if this polyhedron was refined*/
bool DFNPolyhedron::IsRefined(){
    std::pair<const int64_t,int>& firstface = *(fShell.begin());
    int firstface_polyh = fDFN->GetPolyhedralIndex(firstface);
    return firstface_polyh != this->fIndex;
}

/** @brief Remove father from shell and add its subelements */
void DFNPolyhedron::SwapForChildren(TPZGeoEl* father){
    // if no children, return
    int nchildren = father->NSubElements();
    if(!nchildren) return;
    // check if father is actually part of this polyhedron
    auto position = fShell.find(father->Index());
    if(position == fShell.end()) DebugStop();

    // get orientation and remove father
    int orientation = (*position).second;
    fShell.erase(position);

    // enter children
    for(int i=0; i<nchildren; i++){
        int64_t childindex = father->SubElement(i)->Index();
        fShell.insert({childindex,orientation});
    }
}

/** @brief Checks if this polyhedron was intersected by a fracture*/
bool DFNPolyhedron::IsIntersected(DFNFracture& fracture){
    for(auto& face_orient : fShell){
        int64_t index = face_orient.first;
        DFNFace* face = fracture.Face(index);
        if(face) return true;
    }
    return false;
}

/** @brief Returns true if any of the faces in this polyhedron's shell contains only one Inbound rib*/
bool DFNPolyhedron::IntersectsFracLimit(DFNFracture& fracture){
    for(auto& face_orient : fShell){
        int64_t index = face_orient.first;
        DFNFace* face = fracture.Face(index);
        if(!face) continue;
        if(face->NInboundRibs() == 1) return true;
    }
    return false;
}

void DFNPolyhedron::GetEdges(TPZVec<TPZGeoEl*>& edgelist){
    // Throw every edge of every face from fShell into an std::set then copy into edgelist
    std::set<int64_t> edge_set;
    TPZGeoMesh* gmesh = fDFN->Mesh();
    for(auto& face_orient : fShell){
        TPZGeoEl* face = gmesh->Element(face_orient.first);
        TPZManVector<int64_t,4> face_edges = DFN::GetEdgeIndices(face);
        for(int64_t iedge : face_edges){
            edge_set.insert(iedge);
        }
    }
    edgelist.resize(edge_set.size());
    int i=0;
    for(int64_t iedge : edge_set){
        edgelist[i] = gmesh->Element(iedge);
        i++;
    }
}