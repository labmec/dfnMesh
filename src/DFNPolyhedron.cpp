#include "pzgmesh.h"
#include "DFNPolyhedron.h"
#include "DFNMesh.h"
#include "TPZVTKGeoMesh.h"
#include <filesystem>
#include <numeric>

#if PZ_LOG
    static TPZLogger logger("dfn.mesh");
#endif

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
    TPZStack<int64_t> newgels;
    fDFN->MeshPolyhedron(Shell, fCoarseIndex,newgels);
#ifdef PZDEBUG
    CoherentRefinementTest(newgels);
#endif // PZDEBUG
}

#ifdef PZDEBUG
void DFNPolyhedron::CoherentRefinementTest(const TPZVec<int64_t>& newgels){
    constexpr REAL _2deg = 0.2*M_PI/180.;
    constexpr REAL tol = _2deg;
    TPZGeoElSide badangle = DFN::TestInternalAngles(fDFN->Mesh(),newgels,tol);
    if(badangle.Element()){
        const std::string dirpath = "./LOG/SliverInducingVolume";
        std::filesystem::create_directories(dirpath);
        LOGPZ_FATAL(logger,
            "Gmsh created a sliver when refining volume " << fIndex
            << "\nThe first internal angle below the threshold of " << tol
            << " rad, was the side " << badangle.Side() 
            << " of " << badangle.Element()->TypeName() << ' ' << badangle.Element()->Index()
            << "\nThere's a .msh file with the Gmsh model for this volume at ./LOG/gmshAPI_LastVolumeMeshed.msh"
            << "\nI'll plot VTK graphics to " << dirpath << '\n'
        )
        this->PlotVTK(dirpath+"/Volume.vtk");
        this->PlotVTK_NeighbourVolumes(dirpath+"/NeigVolume_");
        fDFN->PlotAllPolygons(dirpath+"/AllPolygons.vtk");
        DFN::PlotVTK_elementList(dirpath+"/SubMesh.vtk",newgels,fDFN->Mesh());
        DFN::PlotVTK_elementList(dirpath+"/CoarseVolume.vtk",{fCoarseIndex},fDFN->Mesh());
        DebugStop();
    }
}
#endif // PZDEBUG

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

void DFNPolyhedron::PlotVTK(const std::string filepath) const {
    
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
        if(copiedgel->HasSubElement()){
            copiedgel->SetSubElement(0,nullptr);
        }
        copiedgel->SetMaterialId(shell.second);
    }
    std::string filename; 
    if(filepath.length() == 0){
        filename = "polyhedron_index_" + std::to_string(fIndex) + ".vtk";
    }else{
        filename = filepath;
    }
    std::ofstream out(filename);
    TPZVTKGeoMesh::PrintGMeshVTK(&bogusMesh, out);
}

void DFNPolyhedron::PlotVTK_NeighbourVolumes(const std::string filepathprefix) const{
    std::string pathprefix = filepathprefix.length() == 0 ? "./LOG/NeigVolume_" : filepathprefix;

    // I SHOULD PROBABLY NOT PLOT MYSELF USING THE SAME PATHPREFIX SO NOT TO INDUCE CONFUSION
    // this->PlotVTK(pathprefix+std::to_string(fIndex));
    
    for(const auto& [faceindex, defaultorient] : fShell){
        int neigVolumeIndex = fDFN->GetPolyhedralIndex(faceindex,defaultorient);
        if(neigVolumeIndex != this->fIndex){
            LOGPZ_WARN(logger, "OrientedFace " << faceindex <<" | "<< defaultorient << " is listed for volume " << this->fIndex << " but thinks it belongs to volume " << neigVolumeIndex);
        }

        for(int i = 0; i < 2; i++){
            const int orient = 2*i-1;
            neigVolumeIndex = fDFN->GetPolyhedralIndex(faceindex,orient);
            if(neigVolumeIndex == this->fIndex) continue;
            if(neigVolumeIndex < 0) continue;
            const DFNPolyhedron& neigVolume = fDFN->Polyhedron(neigVolumeIndex);
            neigVolume.PlotVTK(pathprefix + std::to_string(neigVolumeIndex) + ".vtk");
        }
    }

}

bool DFNPolyhedron::IsBoundedBy(const int64_t faceindex) const{
    return IsBoundedBy(faceindex,1)||IsBoundedBy(faceindex,-1);
}
bool DFNPolyhedron::IsBoundedBy(const std::pair<int64_t, int> orientedFace) const{
    return IsBoundedBy(orientedFace.first,orientedFace.second);
}
bool DFNPolyhedron::IsBoundedBy(const int64_t faceindex, const int orientation) const{
    return fIndex == fDFN->GetPolyhedralIndex(faceindex,orientation);
}
bool DFNPolyhedron::IsBoundedByFather(const int64_t childindex) const{
    TPZGeoEl* father = fDFN->Mesh()->Element(childindex)->Father();
    if(!father){return false;}
    else       {return IsBoundedBy(father->Index());}
}

std::set<int64_t> DFNPolyhedron::GetShellSubset(const std::set<TPZGeoElSide>& delimiter) const{
    TPZGeoMesh* gmesh = fDFN->Mesh();
#ifdef PZDEBUG
    // Coherence
    if(this->IsRefined()) DebugStop();
    if(delimiter.size() < 4) DebugStop();
    for(const auto& delim_side : delimiter){
        const int64_t elindex = delim_side.Element()->Index();
        if(delim_side.Element()->Dimension() != 2) DebugStop();
        if(delim_side.Dimension() != 1) DebugStop();
        if(!this->IsBoundedBy(elindex) && !this->IsBoundedByFather(elindex)) DebugStop();
        // @note: To check if delimiter is definitely consistent, I'd need to have an orientation attributed to the edges (like a Sub-polygon).
        // Better if we trust the user assembled this delimiter from a sub-polygon and therefore everything is consistent
    }
#endif // PZDEBUG

    // Assuming delimiter forms a closed loop of 1D-sides of 2D-elements.
    // From the 2D-elements of the shell of this volume (and their sub-elements)
    // Build the subset of (unrefined) 2D elements all to the same side of this delimiter
    std::set<int64_t> subset;
    std::vector<int64_t> queue;
    queue.reserve(4*fShell.size());
    // Faces in the delimiter are part of the subset, so start by inserting them
    for(auto delim_side : delimiter){
        const auto test = subset.insert(delim_side.Element()->Index());
        // Then queue them for neighbour verification
        if(test.second) queue.push_back(delim_side.Element()->Index());
    }

    const int ndelimiters = queue.size();
    // Loop over faces in queue
    for(int i=0; i < queue.size(); i++){
        TPZGeoEl* currentface = gmesh->Element(queue[i]);
        const int nsides = currentface->NSides();
        // Loop over 1D sides of queued faces
        for(int iside = currentface->FirstSide(1); iside < nsides; iside++){
            TPZGeoElSide gelside {currentface,iside};
            if(gelside.Dimension() != 1) continue;
            // Skip delimiter sides, so the subset is exclusive to the same section isolated by the delimiter
            if(i<ndelimiters && delimiter.find(gelside) != delimiter.end()) continue;
            
            // Loop over unrefined 2D neighbours
            TPZGeoElSide neig = gelside.Neighbour();
            for(/*void*/; neig != gelside; ++neig){
                const TPZGeoEl* neigel = neig.Element();
                if(neig.Element()->Dimension() != 2) continue;
                if(neigel->HasSubElement()) continue;
                const int64_t neigindex = neigel->Index();
                // Add neighbours that bound the same volume
                if(this->IsBoundedBy(neigindex) || this->IsBoundedByFather(neigindex)){
                    const auto test = subset.insert(neigindex);
                    if(test.second) {queue.push_back(neigindex);}
                    break;
                }
            }
        }
    }


#ifdef PZDEBUG
    // If subset contains the same area of fShell, than the input delimiter wasn't a closed loop and this code returned the complete shell instead of a subset
    // const REAL shellArea = this->ComputeArea();
    // REAL subsetArea = 0.;
    // for(int64_t index : subset){ 
    //     TPZGeoEl* gel = fDFN->Mesh()->Element(index);
    //     subsetArea += gel->SideArea(gel->NSides()-1);
    // }
    // if(shellArea - subsetArea < gDFN_SmallNumber) {
        // LOGPZ_ERROR(logger, "at: " << __PRETTY_FUNCTION__ 
        //     <<"\nYou asked me for a subset of the shell of the following volume:"
        //     << *this
        //     <<"\nBut the delimiter of this subset does not form a closed loop of edge-sides. So I can't know what subset you want."
        //     <<"\nDelimiter:\n" << delimiter
        // );
        // std::copy(delimiter.begin(),delimi)
        // DFN::PlotVTK_SideList()
    //     DebugStop();
    // }
#endif // PZDEBUG

    return subset;
}


REAL DFNPolyhedron::ComputeArea() const{
    TPZGeoMesh* gmesh = fDFN->Mesh();

    // auto elArea = [gmesh](const std::pair<int64_t, int> orientedFace){
    //     TPZGeoEl* gel = gmesh->Element(orientedFace.first);
    //     return gel->SideArea(gel->NSides()-1);
    // };
    // REAL area = std::accumulate(fShell.begin(),fShell.end(),0.0,elArea);

    REAL area = 0.;
    for(const auto [index,orient] : fShell){
        TPZGeoEl* gel = gmesh->Element(index);
        area += gel->SideArea(gel->NSides()-1);
    }
    return area;
}
