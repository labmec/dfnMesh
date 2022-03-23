/*! 
 *  @brief     Contains implementation of class DFNFracture
 *  @authors   Pedro Lima
 *  @date      2019
 */

#include "DFNFracture.h"
#include "DFNMesh.h"
#include <math.h>
#include <cstdio>
#include <filesystem>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include "TPZRefPatternDataBase.h"
#include "TPZGeoMeshBuilder.h"
#include "DFNNamespace.h"
#include "DFNGraph.h"

#if PZ_LOG
    // #include "log4cxx/fileappender.h"
    // #include "log4cxx/patternlayout.h"

    static TPZLogger logger("dfn.mesh");
#endif


// Empty Constructor
DFNFracture::DFNFracture(){
}

// Constructor with corner points, a geomesh and material ID
DFNFracture::DFNFracture(DFNPolygon &Polygon, DFNMesh *dfnMesh, FracLimit limithandling, int matid)
    :fPolygon(Polygon),
    fmatid(matid)
{
    fdfnMesh = dfnMesh;
    fLimit = limithandling;
    fIndex = dfnMesh->NFractures();
#ifdef PZDEBUG
    fPolygon.PlotVTK("./LOG/CurrentFracture.vtk",fmatid,fIndex);
#endif // PZDEBUG
}

void DFNFracture::Initialize(DFNPolygon &Polygon, DFNMesh *dfnMesh, FracLimit limithandling, int matid)
{
    fPolygon = Polygon;
    fdfnMesh = dfnMesh;
    fLimit = limithandling;
    fmatid = matid;
    fRibs.clear();
    fFaces.clear();
    fSurfaceFaces.clear();
    fSurfaceEdges.clear();
}

// Copy constructor
DFNFracture::DFNFracture(const DFNFracture &copy){
    this->operator=(copy);
}

// Assignment operator
DFNFracture &DFNFracture::operator=(const DFNFracture &copy){
    fdfnMesh = copy.fdfnMesh;
    fRibs = copy.fRibs;
	fFaces = copy.fFaces;
    fPolygon = copy.fPolygon;
    fmatid = copy.fmatid;
    fLimit = copy.fLimit;
    fSurfaceFaces = copy.fSurfaceFaces;
    fSurfaceEdges = copy.fSurfaceEdges;
    return *this;
}

/**
 * @brief Get the fracture plane
 * @return The plane corner coordinates
 */

DFNPolygon &DFNFracture::Polygon() {
    return fPolygon;
}











DFNFace* DFNFracture::AddFace(DFNFace &face){
    int index= face.Index();
    auto res = fFaces.emplace(index,std::move(face));
    if(res.second == false) return nullptr;
#if PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "[Adding face]\n";
        face.Print(sout,true);
        LOGPZ_DEBUG(logger,sout.str());
    }
#endif
    return &(res.first->second);
}
DFNRib* DFNFracture::AddRib(DFNRib &rib){
    int index= rib.Index();
    auto res = fRibs.insert({index,rib});
    // Check if the rib is already included
    if(res.second == false) return nullptr;
    #if PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "[Adding rib]\n";
        res.first->second.Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
    #endif // PZ_LOG
    return &(res.first->second);
}


DFNRib * DFNFracture::Rib(int64_t index){
    auto candidate = fRibs.find(index);
    if(candidate != fRibs.end()){
        return &candidate->second;
    }
    return nullptr;
}
DFNFace * DFNFracture::Face(int64_t index){
    auto candidate = fFaces.find(index);
    if(candidate != fFaces.end()){
        return &candidate->second;
    }
    return nullptr;
}


void DFNFracture::CreateFaces(){
#if PZ_LOG
    LOGPZ_INFO(logger, "[Start][Searching faces]")
#endif // PZ_LOG
    std::cout<<"\r#Faces intersected = 0";
    TPZGeoMesh *gmesh = fdfnMesh->Mesh();
    TPZGeoEl *gel;
 
    // iterate over 2D elements and check their 1D neighbours for intersections
    int64_t nel = gmesh->NElements();
    for(int iel=0; iel<nel; iel++){
        gel = gmesh->Element(iel);
        if(!gel)continue;
        if(gel->Dimension() != 2) continue;
        if(gel->HasSubElement() ) continue;
 
        int nnodes = gel->NCornerNodes();
        int nedges = nnodes;
 
        // gather ribs
        TPZManVector<int64_t,4>  rib_index(nedges,-1);
        TPZManVector<DFNRib*,4> rib_vec(nedges,nullptr);
        bool is_intersected = false;
        for(int iedge = 0; iedge < nedges; iedge++){
            TPZGeoElSide gelside(gel,iedge+nnodes);
            TPZGeoElSide neig = gelside.Neighbour();
            for(/*void*/;neig != gelside; neig = neig.Neighbour()){
                if(neig.Element()->Dimension() != 1) continue;
                // @TODO shouldn't you check on the material id? It has to be a skeleton matid
                rib_index[iedge] = neig.Element()->Index();
            }
            // rib_index contains the 1D geometric element index (of the skeleton mesh)
            if(rib_index[iedge] == -1) fdfnMesh->DFN_DebugStop(); //Missing 1D skeleton
            // rib_index contains the geometric element index of the edge elements
            // Rib returns the pointer to the DFNRib object IF IT EXISTS
            rib_vec[iedge] = Rib(rib_index[iedge]);
            // the face has at least one rib intersecting
            if(rib_vec[iedge]) is_intersected = true;
        }
        if(!is_intersected) continue;
        
        // build face. Status vector is initialized in the constructor when given a rib_vec
        DFNFace face(gel,this,rib_vec);
        // create a mesh representing the division of the face
        face.UpdateRefMesh();
        // Setup a refinement mesh whose quality measures are checked in DFNFace::NeedsSnap()
        // insert the face in the DFNFracture data structure. A copy of DFNFace is created
        AddFace(face);
//        std::cout<<"\r#Faces intersected = "<<fFaces.size()<<std::flush;
    }

    // Depending on the fracture limit directive (DFNFracture::fLimit) we may have to extend 
    // the fracture surface inside the polyhedra intersected by fracture limits
    if(gmesh->Dimension() == 3 && fLimit != Etruncated) IsolateFractureLimits();
    
    std::cout<<std::endl;
#if PZ_LOG
    LOGPZ_INFO(logger, "[End][Searching faces]")
#endif // PZ_LOG
}























void DFNFracture::CreateRibs(){

    LOGPZ_INFO(logger, "[Start][Searching ribs]");

    std::cout<<"\r\n#Ribs intersected = 0";
    //search gmesh for intersected ribs
    int64_t Nels = fdfnMesh->Mesh()->NElements();
    TPZManVector<int64_t, 2> inode(2,0);
    for (int iel = 0; iel < Nels; iel++){
        TPZGeoEl *gel = fdfnMesh->Mesh()->Element(iel);
        if(!gel) continue;
        //skip all elements that aren't ribs
        if (gel->Dimension() != 1){continue;}
        // skip all elements that have been cut by a previous fracture
        if(gel->HasSubElement()){continue;}

        // Check rib and return the coordinate of the intersection point
        // DFNPolygon::IsCutByPolygon checks if a rib intersects the FracturePolygon by looking if its nodes are on opposite sides of the plane and also checking if intersection point is within bounds of the polygon.
        // There is a SmallNumber tolerance check to handle machine precision, but no 'geometrical tolerance' as its defined in DFNMesh::fTolDist for example.
        TPZManVector<REAL,3> intpoint(3,0);
        bool result = fPolygon.IsCutByPolygon(gel, intpoint);

        // Add rib
        if (result == true){
            DFNRib rib(gel, this);
            rib.SetIntersectionCoord(intpoint);
            AddRib(rib);
//            std::cout<<"\r#Ribs intersected = "<<fRibs.size()<<std::flush;
        }
    }
    std::cout<<std::endl;

    LOGPZ_INFO(logger, "[End][Searching ribs]")

}


void DFNFracture::SetFracMaterial_2D(){
    if(fdfnMesh->Dimension() != 2) return;
    for(auto& itr : fFaces){
        DFNFace& face = itr.second;
        int64_t iline = face.LineInFace();
        if(iline < 0) continue;
        fdfnMesh->Mesh()->Element(iline)->SetMaterialId(fmatid);
    }
}

void DFNFracture::RefineRibs(){
    std::cout << " -Refining ribs\r" << std::flush;
    for(auto itr = fRibs.begin(); itr!=fRibs.end(); itr++){
        DFNRib &rib = itr->second;
        rib.Refine();
    }
    std::cout << "               \r" << std::flush;
}


void DFNFracture::RefineFaces(){
    std::cout << " -Refining faces\r" << std::flush;
    for(auto itr = fFaces.begin(); itr!=fFaces.end(); itr++){
        DFNFace *face = &itr->second;
#if PZ_LOG
        if(logger.isDebugEnabled())
        {
            std::stringstream sout;
            face->Print(sout,true);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif // PZ_LOG
        face->Refine();
    }
    std::cout << "                \r" << std::flush;
    fdfnMesh->CreateSkeletonElements(1);
    if(fdfnMesh->Dimension() < 3){SetFracMaterial_2D();}

}





void DFNFracture::SnapIntersections_ribs(REAL tolDist)
{
    if(fIndex < 0) {LOGPZ_INFO(logger,"[Start][Snapping intersections Ribs]" << "[Orthogonal plane]");}
    else           {LOGPZ_INFO(logger,"[Start][Snapping intersections Ribs]" << "[Fracture " << fIndex << "]");}

    if(tolDist < 0.) tolDist = fdfnMesh->TolDist();
    for(auto itr = fRibs.begin(); itr!=fRibs.end(); itr++){
        DFNRib* rib = &itr->second;
        rib->SnapIntersection_try(tolDist);
    }


    if(fIndex < 0) {LOGPZ_INFO(logger,"[End][Snapping intersections Ribs]" << "[Orthogonal plane]");}
    else           {LOGPZ_INFO(logger,"[End][Snapping intersections Ribs]" << "[Fracture " << fIndex << "]");}
}







void DFNFracture::SnapIntersections_faces(REAL tolDist, REAL tolAngle){
#if PZ_LOG
    if(logger.isInfoEnabled()){
        std::stringstream stream;
        stream <<"[Start][Snapping intersections Faces]";
        if(fIndex < 0)
            {stream << "[Orthogonal plane]";}
        else
            {stream << "[Fracture " << fIndex << "]";}
        LOGPZ_INFO(logger,stream.str());
    }
#endif // PZ_LOG
    
    // Negative entries default to fdfnMesh tolerances
    if(tolDist < 0.) tolDist = fdfnMesh->TolDist();
    if(tolAngle < 0.) tolAngle = fdfnMesh->TolAngle();
    tolAngle = std::cos(tolAngle);
    for(auto &itr : fFaces){
        DFNFace* face = &itr.second;
        face->SnapIntersection_try(tolDist, tolAngle);
    }


#if PZ_LOG
    if(logger.isInfoEnabled()){
        std::stringstream stream;
        stream <<"[End][Snapping intersections Faces]";
        if(fIndex < 0)
            {stream << "[Orthogonal plane]";}
        else
            {stream << "[Fracture " << fIndex << "]";}
        LOGPZ_INFO(logger,stream.str());
    }
#endif // PZ_LOG
}

















/**
 * @brief Get an oriented curve loop in gmsh fashion for a 2D element that has 1D neighbours for all its 1D sides
 * @param shift: inform constant shift if you want to shift element indices (+1 is usually the case for gmsh)
*/
void GetCurveLoop(TPZGeoEl* el, std::vector<int> &loop, const int shift=0){
    if(el->Dimension() != 2) DebugStop();
    int nedges = el->NCornerNodes();
    loop.resize(nedges,-1);
    for(int iside=nedges; iside<2*nedges; iside++){
        TPZGeoElSide gelside(el,iside);
        TPZGeoElSide neig = gelside.Neighbour();
        while(neig.Element()->Dimension() != 1){neig = neig.Neighbour();}
        int orientation = (neig.Element()->NodeIndex(0)==el->NodeIndex(iside-nedges)?1:-1);
        loop[iside-nedges] = orientation*(neig.Element()->Index()+shift);
    }
}








void DFNFracture::SetPolygonIndex(std::pair<int64_t,int> face_orient, int polyg_index,TPZVec<std::array<int, 2>>& Polygon_per_face){
	switch(face_orient.second){
		case  1: Polygon_per_face[face_orient.first][0] = polyg_index; break;
		case -1: Polygon_per_face[face_orient.first][1] = polyg_index; break;
		default: DebugStop();
	}
}
int DFNFracture::GetPolygonIndex(std::pair<int64_t,int> face_orient,const TPZVec<std::array<int, 2>>& Polygon_per_face){
	int polyg_index = -1;
	switch(face_orient.second){
		case  1: polyg_index = Polygon_per_face[face_orient.first][0]; break;
		case -1: polyg_index = Polygon_per_face[face_orient.first][1]; break;
		default: DebugStop();
	}
	return polyg_index;
}
std::pair<int64_t,int> DFNFracture::PolyhNeighbour(std::pair<int64_t,int>& currentface_orient, int current_side, int& neig_side){
    int polyh_index = fdfnMesh->GetPolyhedralIndex(currentface_orient);
    if(polyh_index < 0) DebugStop();
    // Loop through 2D neighbours and find the one that shares the same polyhedron
    TPZGeoEl* gel = fdfnMesh->Mesh()->Element(currentface_orient.first);
    TPZGeoElSide gelside(gel,current_side);
    if(gelside.Dimension() != 1) DebugStop();
    TPZGeoElSide neig = gelside.Neighbour();
    for(/*void*/; neig != gelside; neig = neig.Neighbour()){
        if(neig.Element()->Dimension() != 2) continue;
        int64_t neigindex = neig.Element()->Index();
        // if(!Face(neigindex)) continue;
        neig_side = neig.Side();
        // To share the same subpolygon, they have to share the same polyhedron
        if(      polyh_index == fdfnMesh->GetPolyhedralIndex({neigindex,+1})){
            return {neigindex,+1};
        }else if(polyh_index == fdfnMesh->GetPolyhedralIndex({neigindex,-1})){
            return {neigindex,-1};
        }
    }
    // fdfnMesh->DFN_DebugStop();
    return {-1,0};
}

void DFNFracture::SetLoopOrientation(TPZStack<int64_t>& edgelist){
    int nedges = edgelist.size();
    if(nedges < 3) DebugStop();

    TPZGeoEl* prev_gel = fdfnMesh->Mesh()->Element(edgelist[0]);
    TPZGeoEl* gel = fdfnMesh->Mesh()->Element(edgelist[1]);

    int initialedge_orientation=0;
    if(prev_gel->NodeIndex(1) == gel->NodeIndex(0) ||
       prev_gel->NodeIndex(1) == gel->NodeIndex(1)){
            initialedge_orientation =  1;
    }else{
            initialedge_orientation = -1;
    }

    edgelist[0] *= initialedge_orientation;
    for(int i=1; i<nedges; i++){
        int prev_orientation = (edgelist[i-1]>0?1:-1);
        prev_gel = fdfnMesh->Mesh()->Element(edgelist[i-1]*prev_orientation);
        gel = fdfnMesh->Mesh()->Element(edgelist[i]);
        int prev_node = prev_orientation > 0;
        int orientation = (prev_gel->NodeIndex(prev_node) == gel->NodeIndex(0)?1:-1);
        edgelist[i] *= orientation;
    }
}





TPZGeoEl* DFNFracture::FindPolygon(const TPZStack<int64_t>& polygon) const{
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    
    // Changed on Oct 1st, 2021. See Section 3.3.2 Effects of intersection coalescing on the fracture surface
    // of Pedro's dissertation for discussion and explanations
    TPZGeoEl* loopedface = nullptr;
    
    std::set<int64_t> reducedSubPolygon;
    
        std::set<int64_t> nodeSet; // set with all nodes in polygon. Since it is set it will be of unique entries.
        for(const int64_t rib_index : polygon){
            TPZGeoEl* rib = gmesh->Element(std::abs(rib_index));
            if (rib->NNodes() != 2) DebugStop();
            for (int i = 0; i < 2; i++)
                nodeSet.insert(rib->NodeIndex(i));
        }
        
        for(const int64_t rib_index : polygon){
            TPZGeoEl* rib = gmesh->Element(std::abs(rib_index));
            TPZGeoEl* father = rib->Father();
            int64_t elindextoadd = -1;
            while (father) {
                int count = 0;
                if (father->NNodes() != 2) DebugStop();
                for (int i = 0; i < 2; i++){
                    if (nodeSet.count(father->NodeIndex(i)) != 0)
                        count++;
                }
                if (count == 2) {
                    elindextoadd = father->Index();
                }
                
                father = father->Father();
            }
            if (elindextoadd > -1) {
                // reducedSubPolygon is a set so only unique entries go in
                reducedSubPolygon.insert(elindextoadd);
            }
            else{
                reducedSubPolygon.insert(std::abs(rib_index));
            }
        }
        if (reducedSubPolygon.size() > 4) {
            return nullptr; // there is no 2D element with more than 4 edges
        }
                
        loopedface = DFN::GetLoopedFace(reducedSubPolygon,gmesh);

    return loopedface;
    
}

const bool DFNFracture::CheckSubPolygonPlanarity(TPZStack<int64_t>& subpolygon, const int polyhindex) const {
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    const DFNPolyhedron& volume = fdfnMesh->Polyhedron(polyhindex);
    const int nelements = subpolygon.size();
    if(nelements < 4) return true;
    // if(volume.IsTetrahedron()) return true; // @todo leaving this commented for debug
    TPZManVector<int64_t,4> cornerindices(nelements,-1);
    int i=0;
    for(int64_t edge : subpolygon){
        int64_t index = std::abs(edge);
        if(edge < 0){
            cornerindices[i] = gmesh->Element(index)->NodeIndex(1);
        }else{
            cornerindices[i] = gmesh->Element(index)->NodeIndex(0);
        }
        i++;
    }
    
    TPZFMatrix<REAL> nodeMat(3,nelements);
    i = 0;
    for(auto inod : cornerindices){
        TPZGeoNode& nod = gmesh->NodeVec()[inod];
        for (int j = 0; j < 3; j++) {
            nodeMat(j,i) = nod.Coord(j);
        }
        i++;
    }
    
    DFNPolygon polyg(nodeMat);
    const REAL cosangle = polyg.GetWorstAngleCos();
    const auto limit = M_SQRT2/2.;
    
    if (cosangle < limit) {
//        polyg.Print();
//        PlotVTK_SubPolygon(subpolygon,polyhindex,"FailedSubPolygon");
//        const REAL cosangle = polyg.GetWorstAngleCos();
#if PZ_LOG
        std::stringstream sout;
        sout << "\nSubpolygon with bad planarity found!"
        << "\nCosine of worst angle = " << cosangle << std::endl;
        polyg.Print(sout);
        LOGPZ_INFO(logger,sout.str());
#endif // PZ_LOG
        if(volume.IsTetrahedron() ) DebugStop();
        if(nelements < 4) DebugStop(); 
        return false;
    }
    return true;
}

void DFNFracture::MeshFractureSurface(TPZStack<int> &badVolumes){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    std::cout<<"\r#SubPolygons meshed = 0";
    #if PZ_LOG
        LOGPZ_INFO(logger,"[Start][Meshing Fracture Surface] Frac# " << fIndex);
    #endif // PZ_LOG
    // fdfnMesh->CreateSkeletonElements(1);
    fdfnMesh->ExpandPolyhPerFace();
    // SubPolygons are subsets of the fracture surface contained by a polyhedral volume
    // A subpolygon is formed whenever (at least) 2 DFNFaces, are refined and part of the same polyhedron
    TPZVec<std::array<int, 2>> Polygon_per_face(fdfnMesh->Mesh()->NElements(),{-1,-1});
    TPZStack<int64_t> subpolygon(10,gDFN_NoIndex);
    int polygon_counter = 0;
    // Loop over DFNFaces
    for(auto& itr : fFaces){
        TPZStack<int64_t> newelements;
        DFNFace& initial_face = itr.second;
        if(initial_face.NInboundRibs() < 1) continue;

        for(int i=0; i<2; i++){
            // Setup initial face
            int orientation = 1 - 2*i; // branchless way to do (i==0?1:-1)
            std::pair<int64_t,int> initialface_orient = {initial_face.Index(),orientation};

            // skip 'boundary polyhedron'
            int polyhindex = fdfnMesh->GetPolyhedralIndex(initialface_orient);
            if(polyhindex==0) continue;
            if(polyhindex< 0) fdfnMesh->DFN_DebugStop();
            
            // Skip polyhedra on Fracture limits if limit directive is Etruncated
            if(fLimit==Etruncated && fdfnMesh->Polyhedron(polyhindex).IntersectsFracLimit(*this)) continue; // this can be used to truncate the fracture
            
            // Skip if subpolygon has already been built
            if(GetPolygonIndex(initialface_orient,Polygon_per_face) > -1) continue;

            // Clear container
            #ifdef PZDEBUG
                subpolygon.Fill(gDFN_NoIndex);
            #endif // PZDEBUG 
            subpolygon.clear();
            
            // Get entry side
            int inletside = initial_face.FirstRibSide();
            SetPolygonIndex(initialface_orient,polygon_counter,Polygon_per_face);
            
            // Recursively search for next face until subpolygon loop is closed
            BuildSubPolygon(Polygon_per_face,initialface_orient,inletside,subpolygon);
            polygon_counter++;
            
            // A subpolygon of area zero is not a valid subpolygon and should simply be skipped
            if(DFN::IsValidPolygon(subpolygon,gmesh) == false) continue;

            // Clean up sub-polygon
            SetLoopOrientation(subpolygon);
            std::set<int> locDuplicateIndices;
            const bool hasDuplicateEdges = CheckForDuplicateEdges(subpolygon,locDuplicateIndices);            
            if (hasDuplicateEdges) {
                ClearDuplicateEdges(locDuplicateIndices, subpolygon);
                const bool isClosedLoop = CheckIfPolygonIsClosedLoop(subpolygon);
                if (!isClosedLoop) {
                    badVolumes.Push(polyhindex);
                    continue;
                }
            }

            LOGPZ_DEBUG(logger,"SubPolyg " << polygon_counter <<" : [" << subpolygon << "] in Polyh# " << polyhindex);
            
            bool incorporatedElement = TryFaceIncorporate_Topology(subpolygon,polyhindex,newelements);
            if(incorporatedElement) continue;

            const bool isSubPolPlanarEnough = CheckSubPolygonPlanarity(subpolygon,polyhindex);
            if (!isSubPolPlanarEnough) {
                badVolumes.push_back(polyhindex);
                continue;
            }
            
            MeshPolygon(subpolygon,polyhindex,newelements);
            TryFaceIncorporate_Geometry(subpolygon,polyhindex,newelements,badVolumes);
            // CheckVolumeAngles(subpolygon,polyhindex,newelements,badVolumes);
//            std::cout<<"\r#SubPolygons meshed = "<<polygon_counter<<std::flush;
        }
    }
    if(badVolumes.size() > 0) return; // no need for skeleton elements or buildconnectivity if the mesh is getting rolled back
    // Update connectivity and skeleton of new surface elements, after surface mesh is complete
    std::cout << std::endl;
    std::cout << " -Building connectivity\r" << std::flush;
    fdfnMesh->Mesh()->BuildConnectivity();
    std::cout << "                       \r" << std::flush;
    fdfnMesh->CreateSkeletonElements(1);

    #if PZ_LOG
        LOGPZ_INFO(logger,"[End][Meshing Fracture Surface] Frac# " << fIndex);
    #endif // PZ_LOG
}

const bool DFNFracture::CheckForDuplicateEdges(const TPZStack<int64_t>& subpolygon,
                                               std::set<int>& locDuplicateIndices) const {
    bool hasduplicate = false;
    for (int64_t i = 0; i < subpolygon.size(); i++) {
        for (int64_t j = i + 1; j < subpolygon.size(); j++) {
            if (subpolygon[i] == subpolygon[j]){
                locDuplicateIndices.insert(i);
                locDuplicateIndices.insert(j);
                hasduplicate = true;
            }
        }
    }
    return hasduplicate;
}

void DFNFracture::ClearDuplicateEdges(const std::set<int>& locDuplicateIndices,
                                      TPZStack<int64_t>& subpolygon) const {
    
    if(locDuplicateIndices.size()%2 != 0) DebugStop();
    TPZStack<int> subpolWithoutRepeated;
    for (int i = 0; i < subpolygon.size(); i++) {
        const bool is_in = locDuplicateIndices.find(i) != locDuplicateIndices.end();
        if (!is_in) subpolWithoutRepeated.push_back(subpolygon[i]);
    }
    subpolygon.clear();
    for(auto i : subpolWithoutRepeated) subpolygon.push_back(i);
}

const bool DFNFracture::CheckIfPolygonIsClosedLoop(const TPZStack<int64_t>& subpolygon) const {
    TPZGeoMesh *gmesh = fdfnMesh->Mesh();
    const int nedges = subpolygon.size();
    
    // Two adjacent edges should always share a node
    for(int i=0; i<nedges; i++){
        TPZGeoEl* edge = gmesh->Element(std::abs(subpolygon[i]));
        TPZGeoEl* nextedge = gmesh->Element(std::abs(subpolygon[(i+1)%nedges]));
        const int localnodeID = int(subpolygon[i] > 0);
        const int nextlocalnodeID = int(subpolygon[(i+1)%nedges] < 0);

        const int64_t node = edge->NodeIndex(localnodeID);
        const int64_t nextnode = edge->NodeIndex(nextlocalnodeID);
        // they should share a node
        if(node != nextnode) return false;
    }
    if(nedges < 6) return true;
    // Should also be a non-self-intersecting polygon
    std::unordered_set<int64_t> uniqueNodes;
    for(auto index : subpolygon){
        TPZGeoEl* edge = gmesh->Element(std::abs(index));
        uniqueNodes.insert(edge->NodeIndex(0));
        uniqueNodes.insert(edge->NodeIndex(1));
    }
    const int nnodes = uniqueNodes.size();
    if(nnodes != nedges) return false;

    return true;
}

void DFNFracture::BuildSubPolygon(TPZVec<std::array<int, 2>>& Polygon_per_face,
                                    std::pair<int64_t,int> currentface_orient,
                                    int inlet_side,
                                    TPZStack<int64_t>& subpolygon)
{
    int polyg_index = GetPolygonIndex(currentface_orient,Polygon_per_face);
    if(polyg_index < 0) DebugStop();

    DFNFace* current_dfnface = Face(currentface_orient.first);
    if(!current_dfnface) DebugStop();
    // add line in face to polygon
    int64_t nextedge = current_dfnface->LineInFace();
    subpolygon.push_back(nextedge);

    // Get next face
    int outlet_side = current_dfnface->OtherRibSide(inlet_side);
    int nextinlet_side = -1;
    std::pair<int64_t,int> nextface_orient = PolyhNeighbour(currentface_orient, outlet_side, nextinlet_side);
    if(nextface_orient.second == 0){
        int volumeindex = fdfnMesh->GetPolyhedralIndex(currentface_orient);
        PlotVTK_SubPolygon(subpolygon,volumeindex,"FailedSubPolygon");
        fdfnMesh->Polyhedron(volumeindex).PlotVTK_NeighbourVolumes("./LOG/FailedSubPolygon/NeigVolume_");
        DebugStop();
    }
#ifdef PZDEBUG
    if(!Face(nextface_orient.first)){ std::cout<<"\nPolyhNeighbour returned a next face that was not intersected by the fracture\n"; fdfnMesh->DFN_DebugStop();}
#endif // PZDEBUG

    // Check if its set to polygon
    int nextface_polyg_index = GetPolygonIndex(nextface_orient,Polygon_per_face);
    if(nextface_polyg_index < 0){
        SetPolygonIndex(nextface_orient,polyg_index,Polygon_per_face);
        BuildSubPolygon(Polygon_per_face,nextface_orient,nextinlet_side,subpolygon);
    }
    else if(nextface_polyg_index != polyg_index) DebugStop();
}

// void DFNFracture::GetSubPolygons2(){
//     // SubPolygons are subsets of the fracture surface contained by a polyhedral volume
//     // A subpolygon is formed whenever (at least) 2 DFNFaces, are refined and part of the same polyhedron

//     TPZStack<DFNPolyhedron,20>& polyhedra = fdfnMesh->Polyhedra();
//     if(polyhedra.size() < 2){
//         PZError << "\nUninitialized polyhedra stack in DFNMesh\n";
//         DebugStop();
//     }
//     TPZStack<DFNFace*> dfnfaces(10,nullptr);
//     TPZStack<int64_t> polygon;
//     // @todo: std::vector<int> polygon;
//     for(DFNPolyhedron& polyh : polyhedra){
//         // gather DFNFaces that share the same polyhedron
//         dfnfaces.clear();
//         polyh.ListDFNFaces(this,dfnfaces);

//         if(dfnfaces.size() == 2) continue; // TODO review this condition when we start dealing with fracture boundary

//         int nrefined=0;
//         for(auto face : dfnfaces){
//             nrefined += face->NeedsRefinement();
//         }
//     }
    
// }








/** @brief Projects a non-planar polygon onto its best fitting plane and uses Gmsh to mesh it
 * @param orientedpolygon an oriented loop of edges that don't necessarily occupy the same plane
*/
void DFNFracture::MeshPolygon_GMSH(TPZStack<int64_t>& orientedpolygon, std::set<int64_t>& nodes, TPZStack<int64_t>& newelements, bool isplane){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    const int nnodes = nodes.size();
    const int nedges = nnodes;
    newelements.clear();
    // Project nodes onto best fitting plane
    TPZManVector<REAL,3> centroid(3,0.);
    TPZManVector<REAL,3> normal(3,0.);
    if(!isplane){
        TPZFMatrix<REAL> nodecloud(3,nnodes);
        int j=0;
        for(int64_t inode : nodes){
            TPZManVector<REAL,3> coord(3,0.);
            gmesh->NodeVec()[inode].GetCoordinates(coord);
            for(int i=0; i<3; i++){
                nodecloud(i,j) = coord[i];
            }
            j++;
        }
        DFN::BestFitPlane(nodecloud,centroid,normal);
    }
    
    // gmsh::initialize();
	std::string modelname = "model_polyg";
	gmsh::model::add(modelname);
	gmsh::model::setCurrent(modelname);
	std::string mshfilename = "LOG/gmshAPI_polyg.msh";
    gmsh::option::setNumber("Mesh.Algorithm", 5); // (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
	// Insert nodes ____________________________________
	{TPZManVector<REAL,3> coord(3,0.);
    REAL meshsize = 0.;
    if(isplane){
        for(int64_t inode : nodes){
            gmesh->NodeVec()[inode].GetCoordinates(coord);
            gmsh::model::geo::addPoint(coord[0],coord[1],coord[2],meshsize,inode+gmshshift);
        }
    }else{
        TPZManVector<REAL,3> projcoord(3,0.);
        for(int64_t inode : nodes){
            gmesh->NodeVec()[inode].GetCoordinates(coord);
            projcoord = isplane? coord : DFN::GetProjectedX(coord,centroid,normal);
            gmsh::model::geo::addPoint(projcoord[0],projcoord[1],projcoord[2],meshsize,inode+gmshshift);
        }
	}}
	// Insert lines ____________________________________
    std::vector<int> lineloop;
    lineloop.resize(nedges);
	for(int i=0; i<nedges; i++){
        int64_t iline = abs(orientedpolygon[i]);
		TPZGeoEl *gel = gmesh->Element(iline);
		int64_t node0 = gel->NodeIndex(0)+gmshshift;
		int64_t node1 = gel->NodeIndex(1)+gmshshift;
		gmsh::model::geo::addLine(node0,node1,iline+gmshshift);
		gmsh::model::geo::mesh::setTransfiniteCurve(iline+gmshshift,2);
        int orientation = (orientedpolygon[i] > 0 ? 1 : -1);
        lineloop[i] = (iline+gmshshift)*orientation;
	}
	// Insert faces ____________________________________
    // wiretag is a dummy vector with the shifted index of the face/curve-loop
    std::vector<int> wiretag(1,-1);
    wiretag[0] = gmsh::model::geo::addCurveLoop(lineloop);
    if(lineloop.size() < 5){///< To make it more resistant to new nodes
        gmsh::model::geo::addSurfaceFilling(wiretag,wiretag[0]);
        gmsh::model::geo::mesh::setTransfiniteSurface(wiretag[0]); 
        // gmsh::model::geo::mesh::setRecombine(2,wiretag[0]);
    }else{
        gmsh::model::geo::addPlaneSurface(wiretag,wiretag[0]);
    }
    
    // This line is important so the code can conserve the material ids set by the user in the input coarse mesh. It's important for example that the material id inserted here is different than this->fmatid, otherwise we'll fail the recovery of fracture boundary conditions. There's a more robust way to do it, I know. But it's veeery inneficient.
    gmsh::model::addPhysicalGroup(2,wiretag,DFNMaterial::Eintact);
    // int x = <whatever you want as long as it's not zero>
    // gmsh::model::addPhysicalGroup(2,wiretag, this->fmatid + x);

	
	// synchronize before meshing
	gmsh::model::geo::synchronize();
	// mesh
	gmsh::model::mesh::generate(2);
	#ifdef PZDEBUG
		gmsh::write(mshfilename);
	#endif //PZDEBUG
	// import meshed volume back into PZ geoMesh
	std::set<int64_t>& old_nodes = nodes;
	DFN::ImportElementsFromGMSH(gmesh,2,old_nodes,newelements);
	gmsh::model::remove();
	gmsh::clear();
	// gmsh::finalize();
}

/** @brief Mesh a convex polygon from a list of sequentialy connected edges. If not simple, calls on Gmsh
 * @param polygon a loop of edges that don't necessarily occupy the same plane
 * * @param polyhindex passed only for debugging purposes in case one needs to plot the polyhedron
*/
void DFNFracture::MeshPolygon(TPZStack<int64_t>& polygon, const int polyhindex, TPZStack<int64_t>& newelements){
    
    // New elements to be created
//    TPZStack<int64_t> newelements(1,-1);
    newelements.resize(1);
    newelements.Fill(-1);

    // Get set of nodes
    std::set<int64_t> nodes;
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    for(int64_t line : polygon){
		TPZGeoEl *gel = gmesh->Element(std::abs(line));
		nodes.insert(gel->NodeIndex(0));
		nodes.insert(gel->NodeIndex(1));
	}
    
//    SetLoopOrientation(polygon); // Now done in MeshFractureSurface()
    int nedges = polygon.size();
//    std::cout << "Polygon coords\n";
//    for(auto no : nodes) gmesh->NodeVec()[no].Print();
    // std::cout<<"SubPolygon# "<<polygon_counter<<": "<<polygon<<std::endl;
    // If polygon is planar quadrilateral or triangle, we can skip gmsh
    bool isplane = false;
    try{
        isplane = DFN::AreCoPlanar(gmesh,nodes,1e-5);
        switch(nedges){
            case 0: 
            case 1: 
            case 2: DebugStop();
            case 3:             
            case 4: if(isplane){
                        // If you create elements in the surface with matid = this->fmatid, Fracture Boundary condition recovery may fail
                        DFN::MeshSimplePolygon(gmesh,polygon,DFNMaterial::Eintact,newelements);
                        break;
                    }
            default: MeshPolygon_GMSH(polygon,nodes,newelements,isplane);
        }
    }
    catch(...){
        PlotVTK_SubPolygon(polygon,polyhindex,"FailedSubPolygon");
        PZError << "Failed to mesh a SubPolygon\n"
                << "Oriented edge indices:\n";
        for(int64_t index : polygon)
            PZError << (index>0?' ':'-') << std::abs(index) << '\n';
        PZError << "Plotted SubPolygon to ./LOG/FailedSubPolygon/";
        fdfnMesh->DFN_DebugStop();
    }
#ifdef PZDEBUG
    for(const int64_t index : newelements){
        TPZGeoEl *gel = gmesh->Element(index);
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++) {
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            int count = 0;
            while (neighbour != gelside && count < 1000) {
                neighbour = neighbour.Neighbour();
                count++;
            }
            if(count == 1000)
            {
                DebugStop();
            }
        }
    }
#endif

    InsertFaceInSurface(newelements);

    #if PZ_LOG
        if(logger.isDebugEnabled())
            {LOGPZ_DEBUG(logger,"\tNew elements indices : [" << newelements << "]");}
    #endif // PZ_LOG

}


// void DFNFracture::GetSubPolygons_old(){
//     TPZGeoMesh* gmesh = fdfnMesh->Mesh();
//     std::map<int, TPZAutoPointer<std::vector<int>>> subpolygons_map; // @todo maybe change this to a vector of autopointers...
//     // initialize a data structure to track which subpolygons have included which lines
//     std::map<int64_t, std::pair<int,int>> LineTracker;
//     for(auto iterator=fFaces.begin(); iterator!=fFaces.end(); iterator++){
//         DFNFace* face = &iterator->second;
//         int64_t line = face->LineInFace();
//         if(line == -1) continue;
//         LineTracker[line] = {0,0};
//     }
//     TPZManVector<REAL,3> frac_normal(3,0);
//     fPolygon.GetNormal(frac_normal);
//     for(auto iterator=fFaces.begin(); iterator!=fFaces.end(); iterator++){
//         DFNFace* initial_face = &iterator->second;
//         //Check if face intersection is an actual line or has been coalesced down to a point
//         int64_t firstline = initial_face->LineInFace();
//         int current_line = firstline;
//         if(firstline == -1) continue;
//         //Check if the line in face has already been attributed to all its possible subpolygons
//         std::pair<int,int>* tracker = &LineTracker[firstline];
//         int nloops = int(tracker->first>0) + int(tracker->second>0);
//         if(nloops == 2) continue;
//         int npolyhedra = this->fdfnMesh->FaceTracker[initial_face->Index()];
//         // Decide a direction to follow
//         int direction = (tracker->first>0?0:1);
//         // Track if a direction was tried but failed
//         bool DirectionFailed = false;
//         while(nloops < npolyhedra){
//             // TPZAutoPointer<std::vector<int>> subpolygon = new std::vector<int>;
//             // int debugsize = subpolygons_map.size();
//             // std::printf("#%i ---------\n",debugsize);
//             std::vector<int>* subpolygon = new std::vector<int>;

//             // Find next side based in decided direction
//             int nedges = initial_face->GeoEl()->NCornerNodes();
//             int edge;
//             for(edge = 0; edge<nedges; edge++){
//                 DFNRib* edgerib = initial_face->Rib(edge);
//                 if(!edgerib) continue;
//                 if(edgerib->GeoEl()->NodeIndex(0)==gmesh->Element(firstline)->NodeIndex(direction)) break;
//                 if(edgerib->GeoEl()->NodeIndex(1)==gmesh->Element(firstline)->NodeIndex(direction)) break;
//                 if(edgerib->IntersectionIndex()  ==gmesh->Element(firstline)->NodeIndex(direction)) break;
//             }
//             int nextside = edge + initial_face->GeoEl()->NCornerNodes();

//             // Follow direction by getting the next neighbour with the smallest dihedral angle to close the subpolygon
//             TPZGeoEl* current_face = initial_face->GeoEl();
//             do{ //while(current_face != initial_face->GeoEl());
//                 float angle = DFN::_2PI;
//                 if(current_line != -1){
//                     (*subpolygon).push_back(current_line);
//                     // std::cout<<current_line<<std::endl;
//                     if(current_line > 0)    LineTracker[abs(current_line)].first = 1;
//                     else                    LineTracker[abs(current_line)].second = 1;
//                 }
//                 TPZGeoElSide gelside(current_face,nextside);
//                 TPZGeoElSide neig = gelside.Neighbour();
//                 // Check if gelside's side orientation agree's with fracture normal vector
//                 // @todo there's a better way to do this. One that doesn't depend on the fracture polygon. I've done it in GetPolyhedra, so just adapt it here. The orientation of the in-side is the same of the out-side.
//                 TPZManVector<REAL,3> gelside_vec(3,0);
//                 DFN::GetSideVector(gelside,gelside_vec);
//                 int sideorientation;
//                 if(DFN::DotProduct_f(gelside_vec,frac_normal)>0.){
//                     sideorientation = -1;}
//                 else{
//                     sideorientation = 1;
//                 }
//                 TPZGeoEl* next_face;
//                 int current_side = -1;
//                 for(/*void*/; neig!=gelside; neig = neig.Neighbour()){
//                     if(neig.Element()->Dimension() != 2) continue;
//                     if(!Face(neig.Element()->Index())) continue;
//                     float temp_angle = DFN::DihedralAngle<REAL>(gelside,neig,sideorientation);
//                     if(temp_angle < angle){
//                         angle = temp_angle;
//                         next_face = neig.Element();
//                         current_side = neig.Side();
//                     }
//                 }
//                 if(angle > M_PI+gDFN_SmallNumber){
//                     if(npolyhedra == 1 && !DirectionFailed) {
//                         // Might be an unluckly bad oriented line in initial_face. So try going the other way before DebugStop.
//                         DirectionFailed = true;
//                         direction = (direction+1)%2;
//                         current_line = -firstline;
//                         delete &*subpolygon;
//                         break;
//                     }
//                     PZError << "\nNon-convex regions shouldn't exist at this point\n" << __PRETTY_FUNCTION__ << std::endl;
//                     DebugStop();
//                 }
//                 DirectionFailed = false;
//                 // Get line in face, its orientation and next side
//                 DFNFace* next_dfnface = Face(next_face->Index());
//                 current_line = next_dfnface->LineInFace();
//                 int orientation=0;
//                 for(int iedge=0; iedge < next_face->NCornerNodes(); iedge++){
//                     if(iedge + next_face->NCornerNodes() == current_side) continue;
//                     DFNRib* edge_rib = next_dfnface->Rib(iedge);
//                     if(!edge_rib) continue;
//                     nextside = iedge + next_face->NCornerNodes();
//                     if(current_line != -1){
//                         if(edge_rib->GeoEl()->NodeIndex(0)     ==gmesh->Element(current_line)->NodeIndex(1)) orientation = 1;
//                         else if(edge_rib->GeoEl()->NodeIndex(1)==gmesh->Element(current_line)->NodeIndex(1)) orientation = 1;
//                         else if(edge_rib->IntersectionIndex()  ==gmesh->Element(current_line)->NodeIndex(1)) orientation = 1;
//                         else orientation = -1;
//                     }else{orientation=1;}
//                     break;
//                 }
//                 current_line = orientation*current_line;
//                 current_face = next_face;
//             }while(current_face != initial_face->GeoEl());

//             if(!DirectionFailed){
//                 nloops++;
//                 direction = (direction+1)%2;
//                 // subpolygons_map.insert({subpolygons_map.size(),subpolygon});
//                 if(subpolygon->size() >= 3){ // exception for fractures incorporating an edge of a volume
//                     subpolygons_map[subpolygons_map.size()] = subpolygon;
//                 } 
//             }
//         }
//     }
//     int i_polyg = 0;
//     for(auto iterator=subpolygons_map.begin(); iterator!=subpolygons_map.end(); iterator++){
//         std::vector<int> &polygonloop = *iterator->second;
//         // std::cout<<"\n\nSubPolygon #"<<i_polyg<<"\n";
//         int size = polygonloop.size();
//         for(int iline=0; iline<size; iline++){
//             if(polygonloop[iline] > 0){
//                 std::cout<<" ";
//             }
//             // std::cout<<polygonloop[iline]<<"\n";
//         }
//         i_polyg++;
//     }
// }









/** @brief Identify Ribs, Faces and Polyhedra that are affected by the limits of the fracture*/
void DFNFracture::IsolateFractureLimits(){
    FindOffboundRibs();
    // FindOffboundFaces();
}


void DFNFracture::FindOffboundRibs(){
    // Maybe some consistency checks?
    if(fdfnMesh->Polyhedra().size() < 2) {PZError<<"\nError: Uninitialized polyhedra\n"; DebugStop();}

    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    // TPZStack<TPZGeoEl*> edgelist(10,nullptr);
    // std::set<int64_t> edgelist;
    TPZManVector<REAL,3> intersection(3,0.);
    int npolyh = fdfnMesh->Polyhedra().size();
    // Loop over polyhedra that intersect fracture limits (start at 1 to skip boundary)
    for(int ipoly=1; ipoly<npolyh; ipoly++){
        DFNPolyhedron& polyhedron = fdfnMesh->Polyhedron(ipoly);
        if(polyhedron.IsRefined()) continue;
        if(!polyhedron.IntersectsFracLimit(*this)) continue;
        // Check each edge of this polyhedron for an intersection with the fracture plane
        std::set<int64_t> edgelist = polyhedron.GetEdges();
        for(int64_t index : edgelist){
            TPZGeoEl* edge = gmesh->Element(index);
            if(!fPolygon.IsCutByPlane(edge,intersection)) continue;
            if(Rib(edge->Index())) continue;
            DFNRib rib(edge,this);
            rib.SetIntersectionCoord(intersection);
            rib.FlagOffbound(true);
            // if(edge->MaterialId() != DFNMaterial::Efracture) {edge->SetMaterialId(DFNMaterial::Erefined);}
            // edge->SetMaterialId(-5);
            DFNRib* ribptr = AddRib(rib);
            if(ribptr) ribptr->AppendToNeighbourFaces();
        }
    }
    // SnapIntersections_ribs();

}
void DFNFracture::FindOffboundFaces(){
    DebugStop(); //@todo - Since 2020/nov/06 this is being done by DFNRib::AppendToNeighbourFaces(), so this method will probably be unnecessary
}












void DFNFracture::ExportFractureBC(int matid, std::ofstream& out){
    // DebugStop(); // this is an unreviewed draft
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    int fracdim = gmesh->Dimension()-1;
    int bcdim = fracdim-1;
    std::stringstream stream;
    stream << "\nBCfrac" << this->fIndex << "[] = { ";
    // for(TPZGeoEl* gel : gmesh->ElementVec()){
    TPZGeoEl* gel = nullptr;
    for(int64_t index : this->fSurfaceFaces){
        gel = gmesh->Element(index);
        if(!gel) continue;
        if(gel->Dimension() != fracdim) continue;
        if(gel->HasSubElement()) PZError << "\nYou should call DFNFracture::CleanUp before calling this method.\n";
        int nsides = gel->NSides();
        for(int iside = gel->FirstSide(bcdim); iside < nsides-1; iside++){
            TPZGeoElSide gelside(gel,iside);
            bool IsBoundarySide = true;
            for(TPZGeoElSide neig = gelside.Neighbour(); neig != gelside; ++neig){
                if(neig.Element()->Dimension() != 2) continue;
                if(neig.Element()->HasSubElement()) continue;
                // if(fSurfaceFaces.find(neig.Element()->Index()) != fSurfaceFaces.end()){
                if(neig.Element()->MaterialId() == gelside.Element()->MaterialId()){
                    IsBoundarySide = false;
                    break;
                }
                // TODO: if the matid of the neighbor is equal to one of the fractures that intersect with *this, then it is not a boundary. Something like:
//                this->setOfInterfrac.find()
            }
            if(!IsBoundarySide) continue;
            TPZGeoEl* gelbc = DFN::GetSkeletonNeighbour(gel,iside);
            gelbc->SetMaterialId(matid); // @todo I'll probably want to remove this, right? fracture bc materialid is a concern of the resulting .geo file
            stream << gelbc->Index()+gmshshift << ",";
            /** @todo maybe add to a set?
              * @todo maybe the user would want the material id to match a chosen neighbour? 
              *       this way we would have multiple subsets of the boundary for a fracture, which is likely to come handy
              */
        }
    }
    stream.seekp(stream.str().length()-1);
    stream << "};";
    out << stream.str() << std::endl;
}











// void              ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension, std::set<int64_t> &oldnodes, TPZVec<TPZGeoEl*>& newgels){
// void DFNFracture::ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension, std::set<int64_t> &oldnodes, TPZStack<int64_t> &newelements){
//     // GMsh does not accept zero index entities
//     const int shift = 1;

//     // First, if GMsh has created new nodes, these should be inserted in PZGeoMesh
//     // create a map <node,point>
//     std::map<int,int> mapGMshToPZ;

//     for(int64_t pznode : oldnodes){
// 		std::vector<size_t> node_identifiers;
//         std::vector<double> coord;
//         std::vector<double> parametricCoord;
//         gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord,0,pznode+shift,true);
//         int gmshnode = (int) node_identifiers[0];
// 		// insert with hint (since oldnodes is an already sorted set, these nodes will all go in the end)
//         mapGMshToPZ.insert(mapGMshToPZ.end(),{gmshnode,pznode+shift});
// 	}

//     // add new nodes into PZGeoMesh
//     {
//         // get all nodes from GMsh
//             std::vector<size_t> node_identifiers;
//             std::vector<double> coord;
//             std::vector<double> parametricCoord;
//             gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord);
//         // iterate over node_identifiers
//         int nnodes = node_identifiers.size();
//         for(int i = 0; i < nnodes; i++){
//             int gmshnode = node_identifiers[i];
//             // check if it is contained in the map
//             if(mapGMshToPZ.find(gmshnode) == mapGMshToPZ.end()){
//                 // New node -> add to PZGeoMesh
//                 int pznode = (int) gmesh->NodeVec().AllocateNewElement();
//                 TPZManVector<REAL,3> newnodeX(3);
//                 newnodeX[0] = coord[3*i];
//                 newnodeX[1] = coord[3*i+1];
//                 newnodeX[2] = coord[3*i+2];
//                 gmesh->NodeVec()[pznode].Initialize(newnodeX,*gmesh);
//                 // int pznode = (int) gmesh->NNodes();
//                 // gmesh->NodeVec().resize(pznode+1);
//                 // insert it in map
//                 mapGMshToPZ.insert({gmshnode,pznode+shift});
//             }

//         }
//     }
    

    
//     int64_t nels = gmesh->NElements();
//     std::vector<std::pair<int, int> > dim_to_physical_groups;
//     gmsh::model::getPhysicalGroups(dim_to_physical_groups,dimension);
    
//     /// inserting the elements
//     for (auto group: dim_to_physical_groups) {
       
//         int dim = group.first;
//         // only want elements of a given dimension
//         if(dim != dimension) continue;
//         int physical_identifier = group.second;
       
//         std::vector< int > entities;
//         gmsh::model::getEntitiesForPhysicalGroup(dim, physical_identifier, entities);

// 		for (auto tag: entities) {
// 		// std::cout<<"______________________test - tag = "<<tag;
           
//             std::vector<int> group_element_types;
//             std::vector<std::vector<std::size_t> > group_element_identifiers;
//             std::vector<std::vector<std::size_t> > group_node_identifiers;
//             gmsh::model::mesh::getElements(group_element_types,group_element_identifiers,group_node_identifiers, dim, tag);
//             int n_types = group_element_types.size();
//             for (int itype = 0; itype < n_types; itype++){
//                 int el_type = group_element_types[itype];
//                 int n_nodes = TPZGeoMeshBuilder::GetNumberofNodes(el_type);
//                 std::vector<int> node_identifiers(n_nodes);
//                 int n_elements = group_element_identifiers[itype].size();
//                 for (int iel = 0; iel < n_elements; iel++) {
//                     // int el_identifier = group_element_identifiers[itype][iel]+nels;
//                     int el_identifier = gmesh->CreateUniqueElementId()+gmshshift;
// 					// std::cout<<"\n"<<el_identifier<<"\n";

//                     for (int inode = 0; inode < n_nodes; inode++) {
//                         // Use mapGMshToPZ to translate from GMsh node index to PZ nodeindex
//                         node_identifiers[inode] = mapGMshToPZ[group_node_identifiers[itype][iel*n_nodes+inode]];
//                     }
//                     TPZGeoEl* newel = TPZGeoMeshBuilder::InsertElement(gmesh, physical_identifier, el_type, el_identifier, node_identifiers);
//                     newelements.push_back(newel->Index());
//                     #ifdef PZDEBUG
//                         if(newel->Dimension() != dimension){ 
//                             PZError << "\nGmsh tried to group a " << newel->Dimension() << "D element as " << dimension << "D.\n";
//                             gmesh->Element(el_identifier)->Print(PZError);
//                             DebugStop();
//                         }
//                     #endif // PZDEBUG
// 					// int64_t ntest = gmesh->NElements();
// 					// std::cout<<"nelements = "<<ntest<<"\n";
//                 }
//             }
//         }
//     }
//     gmesh->BuildConnectivity();
// }







void DFNFracture::GetEdgesInSurface(std::set<int64_t>& edges){
    edges.clear();
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    TPZManVector<int64_t,4> edgeindices;
    for(int64_t faceindex : fSurfaceFaces){
        TPZGeoEl* face = gmesh->Element(faceindex);
        edgeindices = DFN::GetEdgeIndices(face);
        for(int64_t edgeindex : edgeindices){
            edges.insert(edgeindex);
        }
    }
    LOGPZ_DEBUG(logger,"Edges in surface [frac " << fIndex << "]:\n" << fSurfaceEdges);
}




void Move(TPZFMatrix<REAL>& cornercoord, const TPZManVector<REAL,3>& centroid, REAL percent){
    int nrows = cornercoord.Rows();
    int ncols = cornercoord.Cols();
    if(nrows != 3) DebugStop();

    percent /= 100.0;

    for(int icoord=0; icoord < 3; icoord++){
        for(int jpoint=0; jpoint < ncols; jpoint++){
            // REAL delta = (cornercoord(icoord,jpoint)-centroid[icoord])*percent;
            // cornercoord(icoord,jpoint) += MIN(tolDist,delta);
            cornercoord(icoord,jpoint) += (cornercoord(icoord,jpoint)-centroid[icoord])*percent;
        }
    }
}






void DFNFracture::CreateOrthogonalFracture(DFNFracture& orthfracture, const int edgeindex){
    

#if PZ_LOG
    if(logger.isDebugEnabled()){
        std::stringstream stream;
        stream << "[Start][Recover limit "<<edgeindex<<"][Fracture " << fIndex << "]";
        LOGPZ_DEBUG(logger,stream.str());
    }
#endif // PZ_LOG

    
    // consistency checks
    if(edgeindex < 0) DebugStop();
    if(edgeindex >= this->fPolygon.NEdges()) DebugStop();
    
    TPZManVector<REAL,3> realnormal(3,0.);

    const TPZFMatrix<REAL>& referencecorners = fPolygon.GetCornersX();

    fPolygon.GetNormal(realnormal);
    REAL edgelength = fPolygon.EdgeLength(edgeindex);
    realnormal[0] *= edgelength;
    realnormal[1] *= edgelength;
    realnormal[2] *= edgelength;

    TPZFMatrix<REAL> cornercoord(3,3,0.);
    // node 0
    cornercoord(0,0) = referencecorners.GetVal(0,(edgeindex+1)%fPolygon.NCornerNodes());
    cornercoord(1,0) = referencecorners.GetVal(1,(edgeindex+1)%fPolygon.NCornerNodes());
    cornercoord(2,0) = referencecorners.GetVal(2,(edgeindex+1)%fPolygon.NCornerNodes());
    // node 1
    cornercoord(0,1) = referencecorners.GetVal(0,edgeindex);
    cornercoord(1,1) = referencecorners.GetVal(1,edgeindex);
    cornercoord(2,1) = referencecorners.GetVal(2,edgeindex);
    // node 2
    cornercoord(0,2) = cornercoord(0,1)+realnormal[0];
    cornercoord(1,2) = cornercoord(1,1)+realnormal[1];
    cornercoord(2,2) = cornercoord(2,1)+realnormal[2];

    // BugFix. When corners of polygon coincide perfectly with nodes in the mesh, orthogonal fractures may fail to perfectly recover a limit. Slightly expanding the fracture (by, say, 0.1%), then letting snap algorithms handle it, worked on all tests I did.
    TPZManVector<REAL,3> centroid(3,0.);
    fPolygon.ComputeCentroid(centroid);
    Move(cornercoord,centroid,0.1);

    TPZGeoMesh* gmesh = this->fdfnMesh->Mesh();
    DFNPolygon dummypolygon(cornercoord,gmesh);

    orthfracture.Initialize(dummypolygon,fdfnMesh,fLimit);
}

/// @todo There's a good room for improvement on this method. I feel like I'm doing 1 or more unnecessary binary searches.
void DFNFracture::UpdateFractureSurface(){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    TPZStack<TPZGeoEl*> to_add;
    TPZStack<TPZGeoEl*> to_remove;
    for(int64_t index : fSurfaceFaces){
        TPZGeoEl* father = gmesh->Element(index);
        if(!father){
            fSurfaceFaces.erase(index);
            continue;
        }
        if(father->Dimension() != 2) DebugStop();
        if(!father->HasSubElement()) continue;
        father->YoungestChildren(to_add);
        to_remove.push_back(father);
    }
    for(TPZGeoEl* gel : to_remove){
        RemoveFromSurface(gel);
    }
    for(TPZGeoEl* gel : to_add){
        AddToSurface(gel);
    }
}

void DFNFracture::RecoverFractureLimits(){
    // fLimit directive decides if this code should run
    if(this->fLimit != FracLimit::Erecovered) return;
    // Nothing to do for fractures that haven't intersected the mesh
    if(fRibs.size() == 0) return;

    LOGPZ_INFO(logger,"[Start][Recover fracture limits][Fracture " << fIndex << "]");

    // int buggylimit = -1;
    DFNFracture& realfracture = *this;
    realfracture.UpdateFractureSurface();
    realfracture.GetEdgesInSurface(fSurfaceEdges);

    // Number of limit edges in this fracture's DFNPolygon
    int nlimits = fPolygon.NCornerNodes();

    for(int ilimit=0; ilimit<nlimits; ++ilimit){
    // for(int ilimit=nlimits-1; ilimit>=0; --ilimit){
    // try{
        DFNFracture orthfracture;
        CreateOrthogonalFracture(orthfracture,ilimit);
        realfracture.GetEdgesInSurface(realfracture.fSurfaceEdges);
        orthfracture.FindRibs(realfracture.fSurfaceEdges);
        orthfracture.SnapIntersections_ribs(fdfnMesh->TolDist());
        orthfracture.SnapIntersections_faces(fdfnMesh->TolDist(),fdfnMesh->TolAngle());
        orthfracture.RefineRibs();
        orthfracture.RefineFaces();
        orthfracture.SortFacesAboveBelow(fmatid,DFNMaterial::Eintact,realfracture);
    // }catch(...){
    //     if(ilimit == buggylimit) DebugStop();
    //     else buggylimit = ilimit;
    //     fdfnMesh->UpdatePolyhedra();
    //     ilimit--;
    //     // Try again with updated polyhedra
    // }

    LOGPZ_INFO(logger,"[ End ][Recover limit "<<ilimit<<"][Fracture " << fIndex << "]");


    }

    // @todo I'm not sure if Updating Polyhedra after all limits is a robust decision. Maybe it should be called after each limit. Although, it hasn't broken yet, so maybe we'll just leave it here until we find an example of this breaking. I've left a TryCatch approach commented above.
    fdfnMesh->UpdatePolyhedra();
    fdfnMesh->Mesh()->BuildConnectivity();


#if PZ_LOG
    LOGPZ_INFO(logger,"[End][Recover fracture limits][Fracture " << fIndex << "]");
#endif // PZ_LOG


}





void DFNFracture::FindRibs(const std::set<int64_t>& ribset){
    if(ribset.size() == 0) return;

    TPZGeoMesh* gmesh = fdfnMesh->Mesh();

    for(int64_t index : ribset){
        TPZGeoEl* gel = gmesh->Element(index);
        if(!gel) continue;
        if(gel->Dimension() != 1) DebugStop();
        if(gel->HasSubElement()) DebugStop();

        TPZManVector<REAL,3> intpoint(3,0.);
        if(fPolygon.IsCutByPlane(gel,intpoint)){
            DFNRib rib(gel,this);
            rib.SetIntersectionCoord(intpoint);
            DFNRib* newrib = AddRib(rib);
            if(newrib) newrib->AppendToNeighbourFaces();
        }
    }
}


void DFNFracture::RemoveFromSurface(const TPZVec<int64_t>& indices){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    for(const int64_t index : indices){
        switch (gmesh->Element(index)->Dimension()){
            case 1: fSurfaceEdges.erase(index); break;
            case 2: fSurfaceFaces.erase(index); break;
            default: DebugStop();
        }
    }
}
void DFNFracture::AddToSurface(const std::set<int64_t>& indices){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    for(const int64_t index : indices){
        switch (gmesh->Element(index)->Dimension()){
            case 1: fSurfaceEdges.insert(index); break;
            case 2: fSurfaceFaces.insert(index); break;
            default: DebugStop();
        }
    }
}
void DFNFracture::AddOrRemoveFromSurface(const std::set<int64_t>& indices){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    for(const int64_t index : indices){
        std::pair<std::set<int64_t>::iterator, bool> test;
        switch (gmesh->Element(index)->Dimension()){
            case 1: {
                test = fSurfaceEdges.insert(index);
                if(!test.second) fSurfaceEdges.erase(test.first);
                break;
            }
            case 2: {
                test = fSurfaceFaces.insert(index);
                if(!test.second) fSurfaceFaces.erase(test.first);
                break;
            }
            default: DebugStop();
        }
    }
}
void DFNFracture::RemoveFromSurface(TPZGeoEl* gel){
    switch (gel->Dimension()){
        case 1: fSurfaceEdges.erase(gel->Index()); break;
        case 2: fSurfaceFaces.erase(gel->Index()); break;
        default: DebugStop();
    }
}
void DFNFracture::AddToSurface(TPZGeoEl* gel){
    switch (gel->Dimension()){
        case 1: fSurfaceEdges.insert(gel->Index()); break;
        case 2: fSurfaceFaces.insert(gel->Index()); break;
        default: DebugStop();
    }
}
void DFNFracture::InsertFaceInSurface(int64_t elindex){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    TPZGeoEl* gel = gmesh->Element(elindex);
    if(!gel) DebugStop();
    TPZStack<TPZGeoEl*> children;
    if(gel->HasSubElement()){
        gel->YoungestChildren(children);
    }else{
        children.push_back(gel);
    }
    for(TPZGeoEl* gel : children){
        // gel->SetMaterialId(fmatid); // Don't set material id here, otherwise fracture BC recovery will fail
        if(gel->Dimension() !=2) DebugStop();
        fSurfaceFaces.insert(gel->Index());
        #if PZDEBUG
            // if(!CheckIsLegalSurfaceElement(gel->Index())){
            //     fdfnMesh->DFN_DebugStop();
            // }
        #endif
    }
}

bool DFNFracture::CheckIsLegalSurfaceElement(const int64_t elindex) const{
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    TPZGeoEl* gel = gmesh->Element(elindex);
    if(!gel) return false;

    for(int iside = gel->FirstSide(1); iside < gel->NSides()-1; iside++){
        TPZGeoElSide gelside(gel,iside);
        TPZGeoElSide neig;
        int nneigh_in_surface = 0;
        for(neig = gelside.Neighbour(); neig != gelside; ++neig){
            nneigh_in_surface += (fSurfaceFaces.find(neig.Element()->Index()) != fSurfaceFaces.end());
        }
        if(nneigh_in_surface > 1){
            PZError << "\n\n[FATAL] Found inconsistent surface\n" << std::endl;
            fdfnMesh->DFN_DebugStop();
            return false;
        }
    }
    return true;
}


void DFNFracture::SortFacesAboveBelow(int id_above, int id_below, DFNFracture& realfracture){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    fPolygon.SortNodes(gmesh);  //< Sort all nodes of the mesh as either being above or below the current orthogonal plane
    realfracture.UpdateFractureSurface(); //< Swap any refined surface element for its children

    std::vector<int64_t> to_remove;
    // Has to be done for all faces on the surface of RealFracture, because we're calling DFNFracture::RecoverFractureLimits after all fractures were inserted into the DFN
    // If you wanna try something different, call RecoverFractureLimts() after each individual fracture has been inserted, then the range of this loop will be:
    // for(auto dfnface : this->fFaces){ int64_t index = dfnface->fGeoEl->Index(); // And you can check the subelements of each dfnface to_add or to_remove
    for(const int64_t index : realfracture.fSurfaceFaces){
        TPZGeoEl* face = gmesh->Element(index);
        if(this->CheckFaceAbove(face,false)) {continue;}
        else {to_remove.push_back(index);}
    }
    for(const int64_t index : to_remove){
        TPZGeoEl* face = gmesh->Element(index);
        realfracture.RemoveFromSurface(face);
    }
}

void DFNFracture::ResetSurfaceMaterial(const int matid){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    for(int64_t index : fSurfaceFaces){
        TPZGeoEl* gel = gmesh->Element(index);
        if(!gel){
            fSurfaceFaces.erase(index);
            continue;
        }
        // if(gel->Dimension() != 2) DebugStop();
        gel->SetMaterialId(matid);
        DFN::SetEdgesMaterialId(gel,matid);
    }
}




void DFNFracture::CleanUp(int surface_matid){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    TPZStack<TPZGeoEl*> to_add;
    TPZStack<TPZGeoEl*> to_remove;
    for(int64_t index : fSurfaceFaces){
        TPZGeoEl* father = gmesh->Element(index);
        if(!father){
            fSurfaceFaces.erase(index);
            continue;
        }
        if(father->Dimension() != 2) DebugStop();
        if(!father->HasSubElement()){
            father->SetMaterialId(surface_matid);
            DFN::SetEdgesMaterialId(father,surface_matid);
            continue;
        }
        father->YoungestChildren(to_add);
        to_remove.push_back(father);
    }
    for(TPZGeoEl* gel : to_add){
        AddToSurface(gel);
        gel->SetMaterialId(surface_matid);
        DFN::SetEdgesMaterialId(gel,surface_matid);
    }
    for(TPZGeoEl* gel : to_remove){
        fSurfaceFaces.erase(gel->Index());
    }

    // // Build the set with every 1D element at the surface of this fracture
    // GetEdgesInSurface(fSurfaceEdges);
}


bool DFNFracture::CheckFaceAbove(TPZGeoEl* face, bool use_face_centroid){
    /// @note If you can guarantee use_face_centroid = true won't lead to errors, go ahead and use it. Else, keep use_face_centroid = false

    TPZManVector<REAL,3> centroid(3,0.);
    TPZGeoMesh* gmesh = face->Mesh();
    // This can be done (less robustly but often cheaper) by checking if centroid of face is above the DFNPolygon that originally defined the fracture
    if(use_face_centroid){
        TPZGeoElSide faceside(face,face->NSides()-1);
        faceside.CenterX(centroid);
        return fPolygon.Compute_PointAbove(centroid);
    // Or it can be done (robust but often costly) by checking the centroid of each edge of the face
    }else{// use_edges_centroid
        int count=0;
        const int nsides = face->NSides();
        for(int iside=face->FirstSide(1); iside < nsides; iside++){
            TPZGeoElSide edgeside(face,iside);
            edgeside.CenterX(centroid);
            count += fPolygon.Compute_PointAbove(centroid);
        }
        return count > 1; // 2 edges above the polygon plane is condition enough for a face to be above the fracture surface even after snap
    }
    DebugStop(); // this shouldn't be reached
    return false;
}




void DFNFracture::Print(std::ostream & out) const
{
    out << "\nFracture #" << fIndex << "\n";
    fPolygon.Print(out);
	out << "\n\nDFNRibs:\n";
    for(auto itr = fRibs.begin(); itr!=fRibs.end(); itr++){
        const DFNRib *rib = &itr->second;
        rib->Print(out);
    }
	out << "\n\nDFNFaces:\n";
    for(auto itr = fFaces.begin(); itr!=fFaces.end(); itr++){
        const DFNFace *face = &itr->second;
        face->Print(out);
    }

    // Surface elements
	out << "\n\nSurface Elements:\n";
    if(fSurfaceFaces.size() < 1) 
        {out << "\"No surface was created/incorporated on this fracture\"";}
    int nelements = fdfnMesh->Mesh()->NElements();
    int width = 2 + int(std::log10(nelements)+1);
    for(int64_t index : fSurfaceFaces){
        out << std::setw(width) << std::right << index << "\n";
    }
    // todo?
    // SubPolygons
}



#if PZ_LOG
    // /** @brief Creates a logger for this object and fills pointer to this->fLogger
    //   * @note A method to create a separate logger + appender for each DFNFracture. (as opposed to the main logger which logs everything from the DFNMesh)*/
    // log4cxx::LoggerPtr DFNFracture::CreateLogger(std::string filename, std::string layout_convpattern){
    //     using namespace log4cxx;
    //     // Default parameters
    //     if(filename == "default") filename = "LOG/dfn.fracture" + to_string(fIndex) + ".log";
    //     if(layout_convpattern == "default") layout_convpattern = "%m%n";

    //     // Create logger
    //     std::string loggername = "dfn.mesh.frac" + to_string(fIndex);
    //     LoggerPtr fLogger = log4cxx::Logger::getLogger(loggername);
        
    //     // Create appender
    //     std::string appendername = "appender_frac" + to_string(fIndex);
    //     FileAppender* appender = new FileAppender();
    //     appender->setAppend(false);
    //     appender->setName(appendername);
    //     appender->setFile(filename);

    //     // Setup layout
    //     PatternLayout* layOut = new PatternLayout();
    //     layOut->setConversionPattern(layout_convpattern);
    //             // "%m%n");
    //             // "[%-5p] %c{1}:%L - %m%n");
    //             // "%d{yyyyf-MM-dd HH:mm:ss} %-5p %c{1}:%L - %m%n");
    //     appender->setLayout(LayoutPtr(layOut));

    //     // Activate options and open log file
    //     log4cxx::helpers::Pool p;
    //     appender->activateOptions(p);
        
    //     // Add appender
    //     fLogger->addAppender(log4cxx::AppenderPtr(appender));

    //     // Set level and additivity
    //     fLogger->setAdditivity(false);
    //     fLogger->setLevel(log4cxx::Level::getInfo());

    //     // LOG4CXX_INFO(fLogger,"test1");
    //     // LOG4CXX_INFO(fLogger,"test2");

    //     return fLogger;
    // }
#endif // PZ_LOG



std::set<int64_t> DFNFracture::IdentifySnapRibs(){
    // Loop over DFNFaces, check if it was intersected through 2 consecutive nodes, if true, get the index of the 1D skeleton connecting those consecutive nodes.
    std::set<int64_t> SnapRibs;
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    // Loop over DFNFaces using fancy C++17 syntax
    for(auto const& [faceindex, face] : fFaces){
        // check if it was intersected through 2 consecutive nodes
        int snapside = face.CheckAssimilatedSide();
        if(snapside < 0) continue;
        // if true, get the index of the 1D skeleton connecting those consecutive nodes
        TPZGeoEl* facegel = gmesh->Element(faceindex);
        TPZGeoEl* SnapRibGel = DFN::GetSkeletonNeighbour(facegel,snapside);

        if(!SnapRibGel){PZError << "\n You may have non-built connectivities or haven't created a necessary skeleton on TPZGeoElSide = {" << facegel->Index() << ", " << snapside << "};"; DebugStop();}

        int64_t SnapRibIndex = SnapRibGel->Index();
        SnapRibs.insert(SnapRibIndex);
    }
    return SnapRibs;
}












std::set<int> DFNFracture::IdentifyIntersectedPolyhedra(){
    std::set<int> IntersectedPolyh;
    for(const auto& volume : fdfnMesh->Polyhedra()){
        if(volume->Index() == 0) continue;           // skip boundary polyhedron
        if(volume->IsRefined()) continue;            // skip refined polyhedra
        if(!volume->IsConvex()) DebugStop();         // they should all be convex at this point
        if(!volume->IsIntersected(*this)) continue;  // get only intersected polyhedra
        IntersectedPolyh.insert(volume->Index());
    }
    return IntersectedPolyh;
}












void DFNFracture::CheckSnapInducedOverlap(){
    {
        TPZGeoMesh* gmesh = fdfnMesh->Mesh();
        std::set<int64_t> SnapRibs = this->IdentifySnapRibs();
        std::set<int> IntersectedPolyh = this->IdentifyIntersectedPolyhedra();
        
        // Gather problematic volumes
        TPZStack<int> problem_volumes;
        for(int vol_index : IntersectedPolyh){
            const DFNPolyhedron& volume = fdfnMesh->Polyhedron(vol_index);
            if(IsProblemVolume(SnapRibs,volume)){
                problem_volumes.push_back(vol_index);
            }
        }

        // Mesh problematic volumes
        if(problem_volumes.size() == 0) return;
        for(int vol_index : problem_volumes){
            DFNPolyhedron& volume = fdfnMesh->Polyhedron(vol_index);
            volume.Refine();
            this->RemoveRefinedDFNFaces(vol_index);
        }

        // Update DFNFracture data
        fdfnMesh->UpdatePolyhedra();
        fPolygon.SortNodes(gmesh);
        this->CreateRibs();
        this->CreateFaces();
        this->SnapIntersections_ribs();
        this->SnapIntersections_faces();
    }
// #ifdef PZDEBUG
//     fdfnMesh->PrintSummary();
// #endif // PZDEBUG
    // for(const int64_t index : SnapRibs){
    //     gmesh->Element(index)->SetMaterialId(DFNMaterial::Eintact);
    // }
    // Continue recursively until there are none problematic volumes
    this->CheckSnapInducedOverlap();
}

bool DFNFracture::IsProblemVolume(const std::set<int64_t>& AllSnapRibs, const DFNPolyhedron& IntersectedVolume) const{
    // Boundary polyhedron can never be a problem volume
    if(IntersectedVolume.Index() == 0) return false;
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();

    /// @return False if is a tetrahedron
    if(IntersectedVolume.IsTetrahedron()) return false;

    /// @return False if volume has N_SnapRibs <= 1
    // @definition: SnapRibs := Set of edges in the surface of this fracture which existed in the mesh before this fracture was inserted.
    // @definition: volumeSnapRibs := Set of SnapRibs that are in the shell of this polyhedral volume
    std::set<int64_t> volumeSnapRibs = IntersectedVolume.GetEdges_InSet(AllSnapRibs);
    int n_snapribs = volumeSnapRibs.size();
    if(n_snapribs <= 1) return false;
    
    /// @return false if the set of volumeSnapRibs form a closed loop of 3 or 4 edges around an existing face in the mesh
    if(DFN::GetLoopedFace(volumeSnapRibs,gmesh)) 
        {return false;}

    /// @return false if the set of rib elders form a closed loop of 3 or 4 edges (which means we've snapped onto a mesh face which will latter be incorporated during DFNFracture::MeshFractureSurface)
    // @definition: rib elders = the set of EldestAncestors of the snap ribs in this polyhedron
    std::set<int64_t> RibElders;
    for(const int64_t rib_index : volumeSnapRibs){
        TPZGeoEl* rib = gmesh->Element(rib_index);
        TPZGeoEl* elder = (rib->Father()?rib->EldestAncestor():rib);
        RibElders.insert(elder->Index());
    }
    switch (RibElders.size()){
        case 1: return false;
        case 2: break; // case 2 is inconclusive, see next test
        case 3:
        case 4: {
            // @return False if VolumeSnapRibs perfectly cover (no more, no less) their eldest ancestors AND the set of rib elders form a closed loop of 3 or 4 edges (which means we've snapped onto a mesh face which will latter be incorporated during DFNFracture::MeshFractureSurface)
            if(DFN::GetLoopedFace(RibElders,gmesh) != nullptr){
                std::set<int64_t> children = DFN::YoungestChildren(RibElders,gmesh);
                if(children == volumeSnapRibs) // @note: set comparison is O(1) if sets have different size
                    {return false;}
            }
            break;
        }
        default: break;
    }

    
    
    /// Count number of neighbour pairs and colinear neighbour pairs of SnapRibs within the volume
    int n_neig_pairs = 0; ///< Number of pairs of neighbours within the set 'volumeSnapRibs'
    int n_colinear_pairs = 0; ///< Number of pairs of COLINEAR neighbours within the set 'volumeSnapRibs'
    for(const int64_t rib_index : volumeSnapRibs){
        TPZGeoEl* rib = gmesh->Element(rib_index);
        for(int inode = 0; inode < 2; inode ++){
            TPZGeoElSide ribside(rib,inode);
            for(TPZGeoElSide neig = ribside.Neighbour(); neig!=ribside; ++neig){
                if(volumeSnapRibs.find(neig.Element()->Index()) == volumeSnapRibs.end()) continue;
                if(neig.Element()->Index() < rib_index) continue;
                n_neig_pairs++;
                TPZGeoEl* elderA = (rib->Father()? rib->EldestAncestor() : rib);
                TPZGeoEl* elderB = (neig.Element()->Father() ? neig.Element()->EldestAncestor() : neig.Element());
                n_colinear_pairs += (int)(elderA == elderB);
            }
        }
    }

    /// @return False if there are no pairs of neighbour SnapRibs
    if(!n_neig_pairs) return false;
    /// @return False if all pairs of neighbour SnapRibs are co-linear (same EldestAncestor)
    if(n_neig_pairs == n_colinear_pairs) return false;

    return true;
}


void DFNFracture::RemoveRefinedDFNFaces(const int vol_index){

    // Refined faces get removed from polyhedra when the polyhedra is meshed and those faces get refined. They're swapped for their children, so we should reach then through these
    // Gather father elements of faces in this polyhedron's shell
    std::set<int64_t> refinedfaces;
    DFNPolyhedron& polyh = fdfnMesh->Polyhedron(vol_index);
    for(const auto& oriented_face : polyh.Shell()){

        int64_t faceindex = oriented_face.first;
        // Get refined faces
        TPZGeoEl* facegel = fdfnMesh->Mesh()->Element(faceindex);
        if(!facegel->Father()) continue;
        auto result = refinedfaces.insert(facegel->FatherIndex());

        // Less readable, but more efficient way to implement this method
        // if(result.second) continue;
        // result.second == false, means we found a father element with at least 2 children in this polyhedron, which makes it a true candidate to be removed
        // truecandidates.insert(facegel->FatherIndex());
    }
    for(int64_t index : refinedfaces){

        // Check if it's a DFNFace
        auto itr = fFaces.find(index);
        if(itr == fFaces.end()) continue;
        #if PZ_LOG
            LOGPZ_DEBUG(logger,"DFNFace # " << index << " \"(Removed from Frac# " << fIndex << ")\"");
        #endif // PZ_LOG
        // Remove it from DFNFaces
        fFaces.erase(itr++); // righthand-increment the iterator to avoid problems with destructor. First copies the iterator to std::map::erase(iterator), then increments, then executes map::erase()
        // fFaces.erase(faceindex);
    }
}

bool DFNFracture::FindFractureIntersection(DFNFracture& OtherFrac, 
                                           TPZStack<int64_t>& EdgeList)
{
    using namespace DFN;
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    EdgeList.clear();

    // Check for user defined real geometrical intersection
    Segment Segment;
    const DFNPolygon& jpolygon = this->Polygon();
    const DFNPolygon& kpolygon = OtherFrac.Polygon();
    bool geom_intersection_Q = jpolygon.ComputePolygonIntersection(kpolygon,Segment);
    if(!geom_intersection_Q){return false;}

    LOGPZ_DEBUG(logger, "[Intersection search] Frac# " << this->Index() << " vs Frac# " << OtherFrac.Index())

    // Get common edges
    const std::set<int64_t> common_edges = DFN::set_intersection(this->fSurfaceEdges,OtherFrac.fSurfaceEdges);
    if(common_edges.size() == 0){
        LOGPZ_DEBUG(logger, "\"This fracture pair has no common edge in their surface. Intersection may have been coalesced into a single node.\"");
        return false;
    }

    // If we assume the set of common edges is always contiguous, the following can be called
    int64_t start, end;
//    ComputeStartAndEndOfSetOfEdges(common_edges, Segment, start, end);
    
    // Because of recoverlimits, the intersection between two fractures
    // can be discontiguous. In this case, we need to apply the shortest
    // path algorithm for each contiguous set of intersectionsd
    TPZManVector<std::set<int64_t>,2> contiguousEdges(1);
    CreateSetsOfContiguousEdges(common_edges,contiguousEdges);
        
    // Build a graph and solve
    for (auto contigCommonEdges : contiguousEdges) {
        DFNGraph graph(fdfnMesh,contigCommonEdges);
        ComputeStartAndEndOfSetOfEdges(contigCommonEdges,Segment,start,end);
        /* If initial and final nodes are the same, then mesh size is too big to catch this intersection
         * If you want to force it, you may try one or more of the following:
         * 1. Pre-refine the mesh (see DFNMesh::PreRefine())
         * 2. Change tolerances (see DFNMesh::fTolDist & DFNMesh::fTolAngle)
         * 3. Design the coarse mesh to have a convenient node that can help to represent this intersection
         */
        if(start == end){
            LOGPZ_DEBUG(logger, "Fracture intersection has been coalesced into a single node. Node index == " << start);
            return false;
        }
        const bool isShortestPathAvailable = graph.ComputeShortestPath(start,end,EdgeList);
        if (!isShortestPathAvailable) {
            std::string filepath = "./LOG/FailedIntersection/";
            PlotVTK_SharedSurface(filepath,OtherFrac,Segment);
        }
    }
    
    return EdgeList.size();
}

void DFNFracture::ComputeStartAndEndOfSetOfEdges(const std::set<int64_t>& common_edges,
                                                 const Segment& Segment,
                                                 int64_t& start, int64_t& end) const {
    using namespace DFN;
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    
    if (common_edges.size() == 1) {
        TPZGeoEl* edge = gmesh->Element(*common_edges.begin());
        start = edge->NodeIndex(0);
        end = edge->NodeIndex(1);
        return;
    }
    
    std::set<int64_t> nodes; // Nodes of the graph
    for(auto index : common_edges){
        TPZGeoEl* edge = gmesh->Element(index);
        nodes.insert(edge->NodeIndex(0));
        nodes.insert(edge->NodeIndex(1));
    }
    
    
    
    TPZStack<int64_t,2> PathNodes(2,-1);
    for(int i=0; i<2; i++){
        TPZManVector<REAL,3> jcoord(3,0.);
        int64_t closestnode = -1;
        REAL closestdist = __FLT_MAX__;

        for(int64_t nodeindex : nodes){
            gmesh->NodeVec()[nodeindex].GetCoordinates(jcoord);
            // Compute distance
            TPZManVector<REAL,3> dif = jcoord - Segment[i];
            REAL normdist = Norm<REAL>(dif);
            if(normdist > closestdist) continue;
            closestdist = normdist;
            closestnode = nodeindex;
        }
        PathNodes[i] = closestnode;
    }
    start = PathNodes[0];
    end = PathNodes[1];
}

void DFNFracture::CreateSetsOfContiguousEdges(std::set<int64_t> common_edges,
                                              TPZManVector<std::set<int64_t>,2> &contiguousEdges) const {
    
    if (contiguousEdges.size() != 1) DebugStop();
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    
    int i = 0; // number of different sets with contiguous edges
    bool keepTrying = true;
    
    while (common_edges.size() != 0) {
        std::set<int64_t>::iterator it = common_edges.begin();
        contiguousEdges[i].insert(*it);
        common_edges.erase(it);
        std::set<int64_t>::iterator itcontig = contiguousEdges[i].begin();

        bool foundNeigh = false;
        for (; itcontig != contiguousEdges[0].end() ; ) {
            TPZGeoEl* currentgel = gmesh->Element(*itcontig);
            if(!currentgel || currentgel->NNodes() != 2) DebugStop();
            TPZGeoElSide currentgelside0(currentgel,0), currentgelside1(currentgel,1);
            it = common_edges.begin();
            for( ; it != common_edges.end() ; it++) {
                TPZGeoEl* gel = gmesh->Element(*it);
                if (!gel || gel->NNodes() != 2) DebugStop();
                TPZGeoElSide gelside0(gel,0), gelside1(gel,1);
                const bool is0 = currentgelside0.NeighbourExists(gelside0);
                const bool is1 = currentgelside0.NeighbourExists(gelside1);
                const bool is2 = currentgelside1.NeighbourExists(gelside0);
                const bool is3 = currentgelside1.NeighbourExists(gelside1);
                if (is0 or is1 or is2 or is3) {
                    contiguousEdges[0].insert(*it);
                    common_edges.erase(it);
                    itcontig = contiguousEdges[0].begin();
                    foundNeigh = true;
                    break;
                }
            }
            if (!foundNeigh)
                itcontig++;
            else{
                foundNeigh = false;
            }
            if (common_edges.size() == 0) {
                break;
            }
        }
        
        if (common_edges.size() != 0) {
            i++;
            contiguousEdges.resize(i+1);
        }
    }
    
    
    
}

void DFNFracture::PlotVTK_SharedSurface(const std::string& filepath, DFNFracture& otherfrac, const Segment& seg) {
    std::filesystem::create_directories(filepath);
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    
    std::set<int64_t> commonFaces = DFN::set_intersection(fSurfaceFaces, otherfrac.Surface());
    std::set<int64_t> commonEdges = DFN::set_intersection(fSurfaceEdges, otherfrac.SurfaceEdges());
    
    // ===> Creation of printable gmesh for the segment between two polygons
    TPZGeoMesh segmentMesh;
    const int64_t npoints = seg.size();
    segmentMesh.ElementVec().Resize(0);
    segmentMesh.NodeVec().Resize(npoints);
    segmentMesh.SetMaxNodeId(npoints-1);
    int eltype = 1; // line
    int bogusmatid = 1; // bogus matid
    TPZGeoNode node_obj;
    for (int i = 0 ; i < npoints ; i++) {
        node_obj.SetCoord(seg[i]);
        node_obj.SetNodeId(i);
        segmentMesh.NodeVec()[i] = node_obj;
        if (i != 0) {
            std::vector<int> nodeident = {i,i+1};
            int index = i-1;
            TPZGeoMeshBuilder::InsertElement(&segmentMesh, bogusmatid,
                                             eltype, index, nodeident);
        }
    }
    std::string outpath = filepath + "Segments.vtk";
    std::ofstream out(outpath);
    TPZVTKGeoMesh::PrintGMeshVTK(&segmentMesh, out, true, true);
    
    // ===> Creation of printable gmesh common edges and facets
    TPZGeoMesh bogusMesh;
    bogusMesh.ElementVec().Resize(gmesh->NElements());
    for (int i = 0; i < bogusMesh.NElements(); i++) {
        bogusMesh.ElementVec()[i] = nullptr;
    }
    bogusMesh.NodeVec() = gmesh->NodeVec();
    
    for (auto index : commonEdges) {
        TPZGeoEl* commonedge = gmesh->Element(index);
        TPZGeoEl* copiededge = commonedge->Clone(bogusMesh);
        copiededge->SetMaterialId(2);
    }
    for (auto index : commonFaces) {
        TPZGeoEl* commonface = gmesh->Element(index);
        TPZGeoEl* copiedface = commonface->Clone(bogusMesh);
        copiedface->SetMaterialId(3);
    }
    
    outpath = filepath + "CommonElements.vtk";
    std::ofstream out2(outpath);
    TPZVTKGeoMesh::PrintGMeshVTK(&bogusMesh, out2, true, true);
    
    // ===> Printing fracture polygons
    outpath = filepath + "Polygon1.vtk";
    this->Polygon().PlotVTK(outpath);
    outpath = filepath + "Polygon2.vtk";
    otherfrac.Polygon().PlotVTK(outpath);
        
    // ===> Printing fracture meshes
    TPZVec<int> matid_backup;
    fdfnMesh->ClearMaterials(GMESHNOMATERIAL, matid_backup);
    outpath = filepath + "FracMesh1.vtk";
    PlotVTK(Index(), outpath, true, false);
    outpath = filepath + "FracMesh2.vtk";
    otherfrac.PlotVTK(otherfrac.Index(), outpath, true, false);

    DebugStop();
}


void DFNFracture::FindFractureIntersection_Trivial(const DFNFracture& OtherFrac, TPZStack<int64_t>& EdgeList){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    std::set<int64_t> EdgeList_set;

    LOGPZ_DEBUG(logger, "[Intersection search][Neighbour MatID] Frac# " << this->Index() << " vs Frac# " << OtherFrac.Index())

    const std::set<int64_t>& othersurface = OtherFrac.fSurfaceFaces;
    
    for(int64_t faceindex : this->fSurfaceFaces){
        TPZGeoEl* face = gmesh->Element(faceindex);
        for(int iside = face->FirstSide(1); iside < face->FirstSide(2); iside++){
            TPZGeoElSide geoside(face,iside);
            for(auto neig = geoside.Neighbour(); neig!=geoside; ++neig){
                if(neig.Element()->Dimension() != 2) continue;
                /* This can be done more efficiently by checking material id, but you'd have to make sure no fractures with different material ids overlap.
                What I mean is, if there's a third fracture that overlaps OtherFrac, then it may have changed the material id of some elements in their shared surface.
                I'm not very enthusiastic on leaving the robustness dependent on details based in material ids, though. Feels hacky.
                So here's the expensive but robust solution: Multiple binary searches. As in, the worst case (which happens whenever 2 fractures do not intersect) being
                twice for every 2D neighbour through each 1D side of each 2D gel in the surface of this fracture by each other fracture tested.*/
                if(othersurface.find(neig.Element()->Index()) == othersurface.end()) continue;
                // if(neig.Element()->MaterialId() != OtherFrac.MaterialId()) continue;
                
                TPZGeoEl* skeleton = DFN::GetSkeletonNeighbour(face,iside);
                EdgeList_set.insert(skeleton->Index());
                break;
            }
        }
    }

    int nedges = EdgeList_set.size();
    if(nedges == 0) return; // No intersection

    EdgeList.clear();
    EdgeList.resize(nedges);
    int i=0;
    for(int64_t index : EdgeList_set){
        EdgeList[i] = index;
        i++;
    }
}

void DFNFracture::SetupGraphicsFractureIntersections(TPZStack<int>& fracfrac_int){
	/// A geometrical intersection between 2 bounded planes in R^3 is a line segment, 
	/// so we represent it by the coordinates of its nodes
	Segment int_segment;
	const int nfrac = fdfnMesh->NFractures();
    TPZGeoMesh* gmesh = this->fdfnMesh->Mesh();

    // Build the set with every 1D element at the surface of this fracture
    GetEdgesInSurface(fSurfaceEdges);

	// Test every pair of fractures for intersection
    int jfrac = this->fIndex;
    for(int kfrac = 0; kfrac<nfrac; kfrac++){
        if(kfrac == jfrac) continue;

        TPZStack<int64_t> intersection_edges;
        this->FindFractureIntersection(*fdfnMesh->FractureList()[kfrac],intersection_edges);
        if(intersection_edges.size() == 0) continue;
        
        int min = std::min(jfrac,kfrac);
        int max = std::max(jfrac,kfrac);
        // You can think of this as the indexing of elements in a symmetric matrix (intersection(i,j)==intersection(j,i)), but shifted 2*N (because N fractures + N boundaries)
        int intersection_matid = nfrac*min - ((min+1)*min)/2 + max + 2*nfrac;

        fracfrac_int.push_back(intersection_matid);
        for(auto index : intersection_edges){
            gmesh->Element(index)->SetMaterialId(intersection_matid);
        }
    }
}

void DFNFracture::SetupGraphicsFractureBC(){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    int fracdim = gmesh->Dimension()-1;
    int bcdim = fracdim-1;
    TPZGeoEl* gel = nullptr;
    int bc_matid = fIndex+fdfnMesh->NFractures();
    for(int64_t index : this->fSurfaceFaces){
        gel = gmesh->Element(index);
        if(!gel) continue;
        if(gel->Dimension() != fracdim) continue;
        if(gel->HasSubElement()) PZError << "\nYou should call DFNFracture::CleanUp before calling this method.\n";
        int nsides = gel->NSides();
        for(int iside = gel->FirstSide(bcdim); iside < nsides-1; iside++){
            TPZGeoElSide gelside(gel,iside);
            bool IsBoundarySide = true;
            for(TPZGeoElSide neig = gelside.Neighbour(); neig != gelside; ++neig){
                if(neig.Element()->Dimension() != 2) continue;
                if(neig.Element()->HasSubElement()) continue;
                // if(fSurfaceFaces.find(neig.Element()->Index()) != fSurfaceFaces.end()){
                if(neig.Element()->MaterialId() == gelside.Element()->MaterialId()){
                    IsBoundarySide = false;
                    break;
                }
            }
            if(!IsBoundarySide) continue;
            TPZGeoEl* gelbc = DFN::GetSkeletonNeighbour(gel,iside);
            // gelbc->SetMaterialId(fmatid_BC);
            gelbc->SetMaterialId(bc_matid);
            /** @todo maybe add to a set?
              * @todo maybe the user would want the material id to match a chosen neighbour? 
              *       this way we would have multiple subsets of the boundary for a fracture, which is likely to come handy
              */
        }
    }
}

void DFNFracture::PlotVTK(const int surface_matid, const std::string exportname,
                          bool putGraphicalElements, bool plotIntersections){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    this->CleanUp(surface_matid);
    TPZStack<int> fracfrac_int;
    // this->ResetSurfaceMaterial(fmatid);
    SetupGraphicsFractureBC();
    if(plotIntersections) SetupGraphicsFractureIntersections(fracfrac_int);

    int mat_BC = fIndex + fdfnMesh->NFractures();
    int mat_frac = fIndex;
    int mat_polygon = -fIndex-1;
    std::set<int> myMaterial {mat_frac, mat_BC, mat_polygon};
    for(int mat_intersection : fracfrac_int) myMaterial.insert(mat_intersection);

    std::ofstream file(exportname);
    TPZVTKGeoMesh::PrintGMeshVTKmy_material(gmesh, file, myMaterial, true, true);
}

void DFNFracture::RollBack(TPZGeoMesh *gmeshBackup) {
    for(auto& ribit : fRibs){
        const int geoindex = ribit.second.GeoEl()->Index();
        ribit.second.SetGeoEl(gmeshBackup->Element(geoindex));
    }
    
    for(auto& faceit : fFaces){
        const int geoindex = faceit.second.GeoEl()->Index();
        faceit.second.SetGeoEl(gmeshBackup->Element(geoindex));
    }
    
    
}
bool DFNFracture::TryFaceIncorporate_Topology(const TPZStack<int64_t>& subpolygon, 
                                            const int polyhindex,
                                            TPZStack<int64_t>& subpolygMesh)
{
    const int nedges = subpolygon.size();
    if(nedges > 4) return false;

    // Check if would-be polygon already exists in the mesh before trying to mesh it
    TPZGeoEl* ExistingGel = FindPolygon(subpolygon);
    if(ExistingGel){
        bool splinterFlag = false;
        if(ExistingGel->HasSubElement()){
            splinterFlag = fSurfaceFaces.find(ExistingGel->SubElement(0)->Index()) != fSurfaceFaces.end();
        }else{
            splinterFlag = fSurfaceFaces.find(ExistingGel->Index()) != fSurfaceFaces.end();
        }
        if(splinterFlag){
            // If this is the second time we're trying to add ExistingGel to the surface, then it's a splinter and it should be removed from the surface
            // Splinter = an element attributed to the fracture surface but kind of loose from it. Like a splinter peeling off of wood
            TPZStack<TPZGeoEl*> children;
            if(ExistingGel->HasSubElement()) {DebugStop(); ExistingGel->YoungestChildren(children);}
            else{children.push_back(ExistingGel);}
            
            for(auto subel : children)
                {RemoveFromSurface(subel);}
            LOGPZ_DEBUG(logger,"\tSplinter removed. Element index : " << ExistingGel->Index());
        }else{
            LOGPZ_DEBUG(logger,"\tIncorporated element index : " << ExistingGel->Index());
            InsertFaceInSurface(ExistingGel->Index());
        }
        return true;
    }
    return false;
}
// void DFNFracture::CheckVolumeAngles(const TPZStack<int64_t>& subpolygon, const int polyhindex, const TPZStack<int64_t>& newelements,TPZStack<int>& badVolumes){
//     // Consistency
// #ifdef PZDEBUG
//     if(newelements.size() < 1) {DebugStop();}
//     if(subpolygon.size() < 3) {DebugStop();}
//     if(polyhindex < 0 || polyhindex >= fdfnMesh->NPolyhedra()) {DebugStop();}
// #endif

//     TPZGeoMesh* gmesh = fdfnMesh->Mesh();
//     DFNPolyhedron& pol = fdfnMesh->Polyhedron(polyhindex);
    
//     // If N Elements == 1, Gmsh was not used to mesh this SubPolygon and we can skip this function
//     // if (newelements.size() == 1) {return;}
    
//     // For each edge in subpolygon, there is a neighbour on the fracture surface. 
//     // Check surfel internal dihedral angle with the 2 faces in the volume shell which are its neighbours
//     // through the 1D-side occupied by the edge in the subpolygon
//     for (int64_t isubpol :subpolygon) {
//         TPZGeoEl* edge = gmesh->Element(abs(isubpol));
//         TPZGeoElSide edgeside(edge,edge->NSides()-1);
//         TPZGeoElSide neig = edgeside.Neighbour();
//         // Get 2D element on the fracture surface
//         TPZGeoElSide surfelside;
//         for (; neig != edgeside ; ++neig) {
//             if (neig.Element()->Dimension() != 2 || neig.Element()->HasSubElement())
//                 {continue;}

//             const int64_t indexneig = neig.Element()->Index();
//             int64_t *p = std::find(newelements.begin(), newelements.end(), indexneig);
//             if (p != newelements.end()){
//                 surfelside = neig;
//                 break;
//             }
//         }
//         if (!surfelside.Element()) {
//             DebugStop();
//         }
        
//         // Search for 2D gel in shell neighbour of surfel through edgeside
//         neig = edgeside.Neighbour();
//         for (; neig != edgeside ; ++neig) {
//             if (neig.Element()->Dimension() != 2 || neig.Element()->HasSubElement())
//                 continue;
            
//             const int64_t neigindex = neig.Element()->Index();
//             // The face we're looking for is either on the volume shell, or is a subelement of an element that is
//             int64_t shellfaceindex = -1;
//             TPZGeoEl* father = neig.Element()->Father();
//             if(pol.IsBoundedBy(neigindex)){
//                 shellfaceindex = neigindex;
//             }else{
//                 if(!father) {continue;}
//                 if(!pol.IsBoundedBy(father->Index())){continue;}
//                 shellfaceindex = father->Index();
//             }

//             // Get orientation of internal angle
//             // If face_in_shell was refined, we'll compute the internal angle to the child element. So we need its orientation.
//             int fatherChildOrientation = 1;
//             if(father && shellfaceindex == father->Index()){
//                 fatherChildOrientation = DFN::SubElOrientation(father,neig.Element()->WhichSubel());
//             }
//             // Orientation of face_in_shell relative to the polyhedral volume
//             const int shellFaceOrientation = pol.IsBoundedBy(shellfaceindex,1) ? 1 : -1;
//             // Orientation of element relative to the edge in the subpolygon
//             // const int faceEdgeOrientation = DFN::OrientationMatch(neig, edgeside) ? 1 : -1;
//             // Proper orientation of the internal angle
//             const int internalAngleOrientation = shellFaceOrientation*fatherChildOrientation;
            
//             // Compute internal angle, and check against a tolerance
//             const REAL internalAngle = DFN::DihedralAngle(neig, surfelside, internalAngleOrientation);
//             constexpr float _5deg = 5.*(M_PI/180.);
//             // constexpr float _355deg = 2.*M_PI - 5.*(M_PI/180.);
//             // if (internalAngle < _5deg || internalAngle > _355deg) {
//             if (internalAngle < _5deg) {
//                 badVolumes.push_back(pol.Index());
//                 return;
//             }
            
//         }
//     }
// }

bool DFNFracture::TryFaceIncorporate_Geometry(const TPZStack<int64_t>& subpolygon,
                                    const int polyhindex,
                                    const TPZStack<int64_t>& newelements,
                                        TPZStack<int>& badVolumes)
{
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();

#ifdef PZDEBUG
    for(const int64_t index : newelements){
        TPZGeoEl *gel = gmesh->Element(index);
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++) {
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            int count = 0;
            while (neighbour != gelside && count < 1000) {
                neighbour = neighbour.Neighbour();
                count++;
            }
            if(count == 1000)
            {
                DebugStop();
            }
        }
    }
#endif
    // Tolerance for incorporation. We should maybe give this as option to the user
    constexpr REAL tol0 = 5.*(M_PI/180.);
    constexpr REAL tol2pi = 2*M_PI - tol0;

    const int nedges = subpolygon.size();
    const std::string plotpath = "./LOG/FailedSubPolygon";
    // Consistency
#ifdef PZDEBUG
    if(newelements.size() < 1) {DebugStop();}
    if(nedges < 3) {DebugStop();}
    if(polyhindex < 0 || polyhindex >= fdfnMesh->NPolyhedra()) {DebugStop();}
#endif

    const DFNPolyhedron& vol = fdfnMesh->Polyhedron(polyhindex);
    
    // For each edge in subpolygon, there is a neighbour on the fracture surface and 2 neighbours on the volume shell, gather those on 3 groups
    TPZManVector<TPZGeoElSide> SurfEl(nedges,{nullptr,-1});          // 2D neighbour of the i-th edge within newelements
    TPZManVector<TPZGeoElSide> PositiveShell(nedges,{nullptr,-1});   // 2D neighbour of the i-th edge on the positive side of the volume shell
    TPZManVector<TPZGeoElSide> NegativeShell(nedges,{nullptr,-1});   // 2D neighbour of the i-th edge on the negative side of the volume shell
    TPZManVector<REAL> PositiveShell_angles(nedges,-999.);
    TPZManVector<REAL> NegativeShell_angles(nedges,-999.);

    int NBadAngles_pos = 0;
    int NBadAngles_neg = 0;
    int NGoodAngles_pos = 0;
    int NGoodAngles_neg = 0;

    try{
    // Loop over edges of the sub-polygon
    for(int iedge=0; iedge < nedges; iedge++){
        int groupcounter = 0;
        const int64_t edgeindex = abs(subpolygon[iedge]);
        const int edgeOrient = DFN::sgn(subpolygon[iedge]);
        TPZGeoEl* edge = gmesh->Element(edgeindex);
        const TPZGeoElSide edgeside(edge,edge->NSides()-1);
        
        // Group neighbours. There should only be one neighbour that matches each of the 3 criteria.
        TPZGeoElSide neig = edgeside.Neighbour();
        for (; neig != edgeside ; ++neig) {
            const TPZGeoEl* neigel = neig.Element();
            if(neigel->Dimension() != 2) continue;
            if(neigel->HasSubElement()) continue;

            const int64_t indexneig = neigel->Index();
            // Check if it's on the shell
            if(vol.IsBoundedBy(indexneig) || vol.IsBoundedByFather(indexneig)){
                groupcounter++;
                // shellFace is the face that is actually on the volume shell. Can either be the neighbour or its father element
                const int64_t shellfaceindex = vol.IsBoundedByFather(indexneig) ? neigel->Father()->Index() : indexneig;
                // If shellFace was refined, we'll compute the internal angle to the child element. So we need its orientation.
                int fatherChildOrientation = 1;
                if(shellfaceindex != indexneig){
                    fatherChildOrientation = DFN::SubElOrientation(neigel->Father(),neigel->WhichSubel());
                }
                // Orientation of shellFace relative to the polyhedral volume
                const int shellFaceOrientation = vol.IsBoundedBy(shellfaceindex,1) ? 1 : -1;
                // Orientation of element relative to the edge in the subpolygon
                const int faceEdgeOrientation = DFN::OrientationMatch(neig, edgeside) ? 1 : -1;
                // Proper orientation of the internal angle
                const int internalAngleOrientation = shellFaceOrientation*fatherChildOrientation;
                // Shell group (positive or negative) depends on the sub-polygon orientation
                const int shellgroup = internalAngleOrientation * faceEdgeOrientation*edgeOrient;

                switch(shellgroup){
                    case  1: PositiveShell[iedge] = neig; break;
                    case -1: NegativeShell[iedge] = neig; break;
                    default: DebugStop();
                }
                continue;
            }

            // Check if it's on the fracture surface
            const int64_t *p = std::find(newelements.begin(), newelements.end(), indexneig);
            if (p != newelements.end()){
                groupcounter++;
                SurfEl[iedge] = neig;
            }
            if(groupcounter > 3) DebugStop();
        }
        if(groupcounter != 3) DebugStop();

        // Compute internal angle, and check against a tolerance
        
        {
            const int faceEdgeOrientation = DFN::OrientationMatch(PositiveShell[iedge], edgeside) ? 1 : -1;
            const int orient = faceEdgeOrientation*edgeOrient;
            PositiveShell_angles[iedge] = DFN::DihedralAngle(PositiveShell[iedge],SurfEl[iedge], orient);
            if(PositiveShell_angles[iedge] < tol0 || PositiveShell_angles[iedge] > tol2pi)
                {NBadAngles_pos++;}
            else
                {NGoodAngles_pos++;}
        }

        {
            const int faceEdgeOrientation = DFN::OrientationMatch(NegativeShell[iedge], edgeside) ? 1 : -1;
            const int orient = -faceEdgeOrientation*edgeOrient;
            NegativeShell_angles[iedge] = DFN::DihedralAngle(NegativeShell[iedge],SurfEl[iedge], orient);
            if(NegativeShell_angles[iedge] < tol0 || NegativeShell_angles[iedge] > tol2pi)
                {NBadAngles_neg++;}
            else
                {NGoodAngles_neg++;}
        }
        if(PositiveShell_angles[iedge] > 3.*M_PI/2.) PositiveShell_angles[iedge] -= 2.*M_PI;
        if(NegativeShell_angles[iedge] > 3.*M_PI/2.) NegativeShell_angles[iedge] -= 2.*M_PI;
        
#ifdef PZDEBUG
        {
            REAL sum = PositiveShell_angles[iedge]+NegativeShell_angles[iedge];
            // the total opening angle should be less than PI, otherwise the volume is not convex
            if(M_PI-sum < -1.e-6) DebugStop();
        }
#endif
        
//        if(PositiveShell_angles[iedge] > M_PI && PositiveShell_angles[iedge] < tol2pi) DebugStop();
//        if(NegativeShell_angles[iedge] > M_PI && NegativeShell_angles[iedge] < tol2pi) DebugStop();
    }
    // If one, but not all, angle violates the threshold, tag this volume as a badVolume
    if((NBadAngles_neg && NGoodAngles_neg && !NBadAngles_pos)
        || (NBadAngles_pos && NGoodAngles_pos && !NBadAngles_neg) ) 
    {
        badVolumes.push_back(vol.Index());
        return false;
    }
    }catch(...){
        PlotVTK_SubPolygon(subpolygon,polyhindex,"FailedSubPolygon"); 
        DFN::PlotVTK_elementList(plotpath+"/SurfaceMesh.vtk",newelements,gmesh);
        DFN::PlotVTK_SideList(plotpath+"/PositiveDelimiter.vtk",PositiveShell);
        DFN::PlotVTK_SideList(plotpath+"/NegativeDelimiter.vtk",NegativeShell);
        DFN::PlotVTK_SideList(plotpath+"/SurfaceDelimiter.vtk",SurfEl);
        DebugStop();
    }
    int closestSubset = 0;
    if(NBadAngles_pos == nedges) closestSubset = +1;
    if(NBadAngles_neg == nedges) closestSubset = -1;
    if(NBadAngles_pos == nedges && NBadAngles_neg == nedges){
        // Accumulate angles
        // closestSubset gets minimum
        // const REAL totalPos = std::reduce( PositiveShell_angles.begin(), PositiveShell_angles.end()
        const REAL totalPos = std::accumulate(  PositiveShell_angles.begin(), PositiveShell_angles.end()
                                            ,0.
                                            ,[](REAL total, REAL angle){
                                                return total + (angle < tol2pi ? angle : 2.*M_PI - angle);
                                            }
                                        );
        const REAL totalNeg = std::accumulate(  NegativeShell_angles.begin(), NegativeShell_angles.end()
                                            ,0.
                                            ,[](REAL total, REAL angle){
                                                return total + (angle < tol2pi ? angle : 2.*M_PI - angle);
                                            }
                                        );
        closestSubset = (totalPos < totalNeg ? +1 : -1);
    }

    if(closestSubset != 0){
        // Get all shell elements that are on the same group (closest subset) delimited by the sub-polygon, and incorporate to surface
        TPZVec<TPZGeoElSide>& delimiter_vector = closestSubset == 1 ? PositiveShell : NegativeShell;
        std::set<TPZGeoElSide> delimiter;
        for(const auto& gelside : delimiter_vector){ delimiter.insert(gelside);}
        std::set<int64_t> shellsubset = vol.GetShellSubset(delimiter);
        AddOrRemoveFromSurface(shellsubset);
        // Unused elements should be discarded
        RemoveFromSurface(newelements);
        for(const int64_t index : newelements){
            gmesh->DeleteElement(gmesh->Element(index));
        }
        return true;
    }


    return false;
}



void DFNFracture::PlotVTK_SubPolygon(const TPZVec<int64_t>& subpolygon, const int volumeindex, const std::string relativefolderpath) const{
    const std::string dirpath = "./LOG/"+relativefolderpath;
    std::filesystem::create_directories(dirpath);

    TPZGeoMesh* gmesh = fdfnMesh->Mesh();

    // Plot subpolygon edges
    {
        std::ofstream outstream(dirpath+"/SubPolygon.vtk");
        std::set<int64_t> subpolygon_set;
        for(int64_t index : subpolygon) subpolygon_set.insert(std::abs(index));
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,subpolygon_set,outstream);
    }

    // Plot polyhedral volume
    fdfnMesh->Polyhedron(volumeindex).PlotVTK(dirpath+"/Volume.vtk");

    // Plot fracture plane
    fPolygon.PlotVTK(dirpath+"/FracturePlane.vtk",fmatid,fIndex);
}
