#include "pzgmesh.h"
#include "DFNPolyhedron.h"
#include "DFNMesh.h"
#include "TPZVTKGeoMesh.h"

/// Assignment operator
DFNPolyhedron &DFNPolyhedron::operator=(const DFNPolyhedron &copy){
    fShell.clear();
    fShell = copy.fShell;
    fIndex = copy.fIndex;
    fDFN = copy.fDFN;
    fCoarseIndex = copy.fCoarseIndex;
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
void DFNPolyhedron::Print(std::ostream& out, bool detailed) const{
    int nelements = this->fDFN->Mesh()->NElements();
    int width = 2 + int(std::log10(nelements)+1);
    out << "\nPolyh#"<<std::setw(width-1)<<fIndex<<":";
    out << (IsRefined()?"(R)":"   ");
    for(auto& oriented_face : fShell){
        // out << std::setw(5) << std::showpos << oriented_face << " ";
        out << std::setw(width) << oriented_face;
    }
    if(!detailed) return;
    out << "\n\tCoarseIndex = " << std::setw(width-1) << fCoarseIndex
        << "\n\tConvex      = " << std::setw(width-1) << (fIsConvex?"true":"false");
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
bool DFNPolyhedron::IsRefined()const{
    const auto& firstface = *(fShell.begin());
    int firstface_polyh = fDFN->GetPolyhedralIndex(firstface);
    return firstface_polyh != this->fIndex || firstface_polyh == -1;
}

/** @brief Remove father from shell and add its subelements */
void DFNPolyhedron::SwapForChildren(TPZGeoEl* father){
    /// @todo I don't have time to check now, but I think this function keeps getting called on polyhedra that have stopped being important. So there might be room for optimization fixing this.
    // if no children, return
    int nchildren = father->NSubElements();
    if(!nchildren) return;
    // check if father is actually part of this polyhedron
    auto position = fShell.find(father->Index());
    if(position == fShell.end()) DebugStop();

    // get orientation and remove father
    int fatherorientation = (*position).second;
    fShell.erase(position);

    // enter children
    for(int i=0; i<nchildren; i++){
        int64_t childindex = father->SubElement(i)->Index();
        int orientation = fatherorientation * DFN::SubElOrientation(father,i);
        // fShell.insert({childindex,orientation});
        fShell[childindex] = orientation;
    }
}

/** @brief Checks if this polyhedron was intersected by a fracture*/
bool DFNPolyhedron::IsIntersected(DFNFracture& fracture)const{
    for(auto& face_orient : fShell){
        int64_t index = face_orient.first;
        DFNFace* face = fracture.Face(index);
        if(face) return true;
    }
    return false;
}

/** @brief Returns true if any of the faces in this polyhedron's shell contains only one Inbound rib*/
bool DFNPolyhedron::IntersectsFracLimit(DFNFracture& fracture)const{
    for(auto& face_orient : fShell){
        int64_t index = face_orient.first;
        DFNFace* face = fracture.Face(index);
        if(!face) continue;
        if(face->NInboundRibs() == 1) return true;
    }
    return false;
}

std::set<int64_t> DFNPolyhedron::GetEdges() const{
    // Throw every edge of every face from fShell into an std::set
    std::set<int64_t> edge_set;
    TPZGeoMesh* gmesh = fDFN->Mesh();
    for(auto& face_orient : fShell){
        TPZGeoEl* face = gmesh->Element(face_orient.first);
        TPZManVector<int64_t,4> face_edges = DFN::GetEdgeIndices(face);
        for(int64_t iedge : face_edges){
            edge_set.insert(iedge);
        }
    }
    return edge_set; // move semantics or RVO will handle this. Right?
}

std::set<int64_t> DFNPolyhedron::GetEdges_InSet(const std::set<int64_t>& SuperSet) const{
    std::set<int64_t> result;
    TPZGeoMesh* gmesh = fDFN->Mesh();
    const auto& end = SuperSet.end();
    for(auto& face_orient : fShell){
        TPZGeoEl* face = gmesh->Element(face_orient.first);
        TPZManVector<int64_t,4> face_edges = DFN::GetEdgeIndices(face);
        for(int64_t iedge : face_edges){
            if(SuperSet.find(iedge) != end)
                {result.insert(iedge);}
        }
    }
    return result; // move semantics or RVO will handle this. Right?
}

void DFNPolyhedron::Refine(){
    TPZStack<std::pair<int64_t,int>> Shell(fShell.size(),{-1,0});
    int i=0;
    for(auto oriented_face : fShell){
        Shell[i].first = oriented_face.first;
        Shell[i].second = oriented_face.second;
        i++;
    }
    fDFN->MeshPolyhedron(Shell, fCoarseIndex);
}


bool DFNPolyhedron::IsTetrahedron() const{
    bool condition = (fShell.size() == 4 
                        && fDFN->Mesh()->Element(fShell.begin()->first)->Type() == MElementType::ETriangle
                    );

    return condition;

    // Don't really need these, right?
    // TPZGeoMesh* gmesh = fDFN->Mesh();
    // int ntriangles = 0;
    // for(auto& orientedface : fShell){
    //     TPZGeoEl* face = gmesh->Element(orientedface.first);
    //     ntriangles += (face->Type() == MElementType::ETriangle);
    // }

    // return ntriangles == 4;
}

void DFNPolyhedron::PrintVTK() const {
    
    TPZGeoMesh bogusMesh;
    TPZGeoMesh *gmesh = fDFN->Mesh();
    bogusMesh.ElementVec().Resize(gmesh->NElements());
    for (int i = 0; i < bogusMesh.NElements(); i++) {
        bogusMesh.ElementVec()[i] = nullptr;
    }
    bogusMesh.NodeVec() = gmesh->NodeVec();
    
    std::set<int64_t> elindexes;
    for (auto shell : fShell) {
        TPZGeoEl *gel = gmesh->Element(shell.first);
        TPZGeoEl *copiedgel = gel->Clone(bogusMesh);
        copiedgel->SetMaterialId(shell.second);
    }
    std::string filename = "polyhedron_index_" + std::to_string(fIndex) + ".vtk";
    std::ofstream out(filename);
    TPZVTKGeoMesh::PrintGMeshVTK(&bogusMesh, out);
}
