/*! 
 *  @brief     Contains implementation of class DFNFracture
 *  @authors   Pedro Lima
 *  @date      2019
 */

#include "DFNFracture.h"
#include "DFNMesh.h"
#include <math.h>
#include <cstdio>
// #include <unordered_set>
#include <algorithm>
#include "TPZRefPatternDataBase.h"
#include "TPZGeoMeshBuilder.h"

// const float _2PI = 6.2831853071795865;

// Empty Constructor
DFNFracture::DFNFracture(){
}

// Constructor with corner points, a geomesh and material ID
DFNFracture::DFNFracture(DFNPolygon *Polygon, DFNMesh *dfnMesh){
    fPolygon = Polygon;
    fdfnMesh = dfnMesh;

    // Set corner nodes of fracture into mesh
    if(fdfnMesh->Dimension() == 3){
        TPZManVector<int64_t,4> nodeindices = fPolygon->SetPointsInGeomesh(fdfnMesh->Mesh());
    }
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
    return *this;
}

/**
 * @brief Get the fracture plane
 * @return The plane corner coordinates
 */

DFNPolygon *DFNFracture::Polygon() {
    return fPolygon;
}











void DFNFracture::AddFace(DFNFace &face){
    int index= face.Index();
    fFaces.emplace(index,face);
}
void DFNFracture::AddRib(DFNRib &rib){
    int index= rib.Index();
    fRibs.emplace(index,rib);
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


void DFNFracture::FindFaces(){
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
        for(int iedge = 0; iedge < nedges; iedge++){
            TPZGeoElSide gelside(gel,iedge+nnodes);
            TPZGeoElSide neig = gelside.Neighbour();
            for(/*void*/;neig != gelside; neig = neig.Neighbour()){
                if(neig.Element()->Dimension() != 1) continue;
                rib_index[iedge] = neig.Element()->Index();
            }
        }

        // build face
        DFNFace candidate(gel,this);
        candidate.SetRibs(rib_index);
        candidate.UpdateStatusVec();
        if(!candidate.IsIntersected()) continue;
        if(candidate.IsOnBoundary()){
            candidate.FindInPlanePoint();
        }
        // gel->SetMaterialId(DFNMaterial::Erefined);
        candidate.UpdateRefMesh();
        AddFace(candidate);
    }
}























void DFNFracture::FindRibs(){
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

        // Check rib
        TPZManVector<REAL,3> intpoint(3,0);
        bool resul = fPolygon->Check_rib(gel, &intpoint);

        // Add rib
        if (resul == true){
            DFNRib rib(gel, this);
            rib.StatusVec()[2] = 1; //StatusVec = {0,0,1}
            rib.SetIntersectionCoord(intpoint);
            if(gel->MaterialId() != DFNMaterial::Efracture) {gel->SetMaterialId(DFNMaterial::Erefined);}
            AddRib(rib);
        }
    }
}


void DFNFracture::RefineRibs(){
    for(auto itr = fRibs.begin(); itr!=fRibs.end(); itr++){
        DFNRib* rib = &itr->second;
        rib->Refine();
    }
}
void DFNFracture::RefineFaces(){
    for(auto itr = fFaces.begin(); itr!=fFaces.end(); itr++){
        DFNFace *face = &itr->second;
        face->Refine();
    }
}





void DFNFracture::OptimizeRibs(REAL tolDist){
    for(auto itr = fRibs.begin(); itr!=fRibs.end(); itr++){
        DFNRib* rib = &itr->second;
        rib->Optimize(tolDist);
    }
}
void DFNFracture::OptimizeFaces(REAL tolDist, REAL tolAngle){
    for(auto itr = fFaces.begin(); itr!=fFaces.end(); itr++){
        DFNFace* face = &itr->second;
        face->Optimize(tolDist, tolAngle);
    }
}













void DFNFracture::GetOuterLoop(std::vector<int> &loop){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    // Compute fracture edges' lenghts
    int nedges = fPolygon->GetCornersX().Cols();
    TPZVec<REAL> edgelength(nedges,0);
    Matrix fraccorners(3,nedges);
    fraccorners = fPolygon->GetCornersX();
    for(int i = 0; i < nedges; i++){
        edgelength[i] = sqrtl(pow(fraccorners(0,i)-fraccorners(0,(i+1)%nedges),2)
                              +pow(fraccorners(1,i)-fraccorners(1,(i+1)%nedges),2)
                              +pow(fraccorners(2,i)-fraccorners(2,(i+1)%nedges),2));
    }
    
	// vector of pointers to maps
    // @todo refactor this to tpzautopointer<map> to prevent memory leak
	// TPZManVector<std::map<REAL, int64_t>* > edgemap(nedges);
	TPZManVector<TPZAutoPointer<std::map<REAL, int64_t>> > edgemap(nedges);
    for(int i = 0; i < nedges; i++){
        edgemap[i] = new std::map<REAL, int64_t>;
    }
    // The set of points that have already been checked
    std::set<int64_t> checked;

    // iterate over all endfaces and map it to the fracture-edge that cuts it
    for (auto it = fFaces.begin(); it != fFaces.end(); it++){
        // get intersection node index and coordinates
        DFNFace *iface = &it->second;
        // @todo iface->IsOnBoundary2
        TPZManVector<REAL,3> ipointcoord = iface->IntersectionCoord();
        if(ipointcoord.size()<2) continue;
        int64_t ipointindex = iface->IntersectionIndex();

        // check if point was already checked
        auto aux = checked.insert(ipointindex);
        bool already_checked = !aux.second;
        if(already_checked) continue;

        // @todo if(ipointindex == any of the polygon.CornerIndex) continue;
        // iterate over edges to check if ipoint belongs to it
        for(int iedge = 0; iedge < nedges; iedge++){
			//vectors from ipoint to iedge's nodes
			TPZManVector<REAL, 3> v1(3);
				v1[0] = fraccorners(0,iedge) - ipointcoord[0];
				v1[1] = fraccorners(1,iedge) - ipointcoord[1];
				v1[2] = fraccorners(2,iedge) - ipointcoord[2];
			TPZManVector<REAL, 3> v2(3);
				v2[0] = fraccorners(0,(iedge+1)%nedges) - ipointcoord[0];
				v2[1] = fraccorners(1,(iedge+1)%nedges) - ipointcoord[1];
				v2[2] = fraccorners(2,(iedge+1)%nedges) - ipointcoord[2];
			// square of cross product
			REAL temp = pow(v1[1]*v2[2] - v1[2]*v2[1],2);
				temp += pow(v1[2]*v2[0] - v1[0]*v2[2],2);
				temp += pow(v1[0]*v2[1] - v1[1]*v2[0],2);
				temp = sqrtl(temp);
            // check if point is in edge by calculating if it's normal distance to the edge is zero
            REAL dist = temp/edgelength[iedge];
			if(dist<gDFN_SmallNumber){
				// compute local 1D coordinate (alpha)
				REAL norm = sqrtl(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
				REAL alpha = norm/edgelength[iedge];
				if(alpha > 1+gDFN_SmallNumber){std::cout<<"\nMisattribution of point to edge\n";DebugStop();}
				// map intersection indexes from smallest alpha to biggest
                // map <alpha, index>
				edgemap[iedge]->insert({alpha, ipointindex});
			}
        }
    }
    bool warning_message_Q = true;
    //Once intersections on fracture-edges have been properly ordered and mapped by edge
	//iterate over edges to split them
	for (int iedge = 0; iedge < nedges; iedge++)
	{
		int64_t nels = gmesh->NElements();
		TPZManVector<int64_t,2> inodes(2);     //index of nodes to be connected
        int icorner = iedge; //for readability

		// connect first end-face intersection to iedge's first node
        inodes[0] = fPolygon->CornerIndex(icorner);
		auto it = edgemap[iedge]->begin();
        if(edgemap[iedge]->size() == 0){
            #ifdef PZDEBUG
                if(warning_message_Q){
                    PZError<<"\n Warning: Is there an edge of a fracture that doesn't cut any element? \n\n";
                    warning_message_Q = false;
                }
            #endif //PZDEBUG
            // DebugStop();
        }
        // iterate over iedge's map
        while(it != edgemap[iedge]->end()){
            inodes[1] = it->second;
            gmesh->CreateGeoElement(EOned, inodes, DFNMaterial::Efracture, nels);
            loop.push_back((int) nels);
            inodes[0] = inodes[1];
            it++;
        }

		// connect last end-intersection to edge last node
        inodes[1] = fPolygon->CornerIndex((icorner+1)%nedges);
		gmesh->CreateGeoElement(EOned, inodes, DFNMaterial::Efracture, nels);
        loop.push_back((int) nels);
    }
	
    gmesh->BuildConnectivity();
    // correct duplicates
    int nlines = loop.size();
    TPZGeoEl* gel;
    for(int iline=0; iline<nlines; iline++){
        int iedge = loop[iline];
        gel = gmesh->Element(iedge);
        if(!gel || gel->Dimension() != 1) DebugStop();
        TPZGeoElSide gelside(gel,2);
        TPZGeoElSide neig = gelside.Neighbour();
        bool newedge_Q = true;
        for(/*void*/; neig!=gelside; neig = neig.Neighbour()){
            if(neig.Element()->Dimension() != 1) continue;
            // if duplicate, replace in loop and delete duplicate
            newedge_Q = false;
            loop[iline] = neig.Element()->Index(); 
            gmesh->DeleteElement(gel,iedge);
        }
        // if it's new, add it to fOutline
        if(newedge_Q){fOutline.insert({iedge,gel});}
    }

    // fix orientation
    for(int iline=1; iline<nlines; iline++){
        int iline_index = loop[iline];
        gel = gmesh->Element(iline_index);
        int anterior_index = loop[(iline+nlines-1)%nlines];
        TPZGeoEl *anterior = gmesh->Element(anterior_index);
        if(gel->NodeIndex(0) != anterior->NodeIndex(1)){
            loop[iline] *= -1;
        }
    }
}













void DFNFracture::GetFacesInSurface(std::vector<TPZGeoEl*> &faces){
    faces.reserve(20);
    std::map<int64_t, int> candidates;
    TPZGeoEl *cand_face;
    int nedges;
    TPZGeoEl *edge;
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    // loop through edges in outline and check their 2D neighbours
    for(auto itr = fOutline.begin(); itr != fOutline.end(); itr++){
        edge = itr->second;
        TPZGeoElSide gelside(edge,2);
        TPZGeoElSide neig = gelside.Neighbour();
        for(/*void*/;neig!=gelside; neig = neig.Neighbour()){
            cand_face = neig.Element();
            if(cand_face->Dimension() != 2) continue;
            nedges = cand_face->NCornerNodes();
            int64_t index = cand_face->Index();
            if(++candidates[index] < nedges) continue;
            // @todo{ if 2 or more edges of cand_face belong to fracture outerloop
                // cand_face will not enter;
                // the matching edges should be swaped by the remaining edges of cand_face
                // pay attention to the orientation  }
            faces.push_back(gmesh->Element(index));
            SetMaterialIDChildren(DFNMaterial::Efracture,cand_face);
            fSurface.insert({index,cand_face});
        }
    }
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







void DFNFracture::MeshFractureSurface(){
    fdfnMesh->Mesh()->BuildConnectivity();
    // GMsh does not accept zero index entities
    const int shift = 1;
    // First construct the edges of the fracture surface
    std::vector<int> outerLoop;
    GetOuterLoop(outerLoop);
    std::vector<TPZGeoEl*> facesInSurface;
    GetFacesInSurface(facesInSurface);
    
    // initialize GMsh
    // gmsh::initialize();
	std::string modelname = "modelsurface";
	gmsh::model::add(modelname);
    gmsh::model::setCurrent(modelname);
    gmsh::option::setNumber("Mesh.Algorithm", 5); // (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
    // INSERT POINTS
        // iterate over fOutline and get points
    std::set<int64_t> pointset;
    {
        int64_t index = -1;
        TPZManVector<REAL,3> projcoord(3,0);
        TPZManVector<REAL,3> realcoord(3,0);
        for(auto itr : fOutline){
            TPZGeoEl* gel = itr.second;
            for(int inode=0; inode<gel->NCornerNodes(); inode++){
                index = gel->NodeIndex(inode);
                bool newpoint = pointset.insert(index).second;
                if(!newpoint){continue;}
                gel->NodePtr(inode)->GetCoordinates(realcoord);
                projcoord = fPolygon->GetProjectedX(realcoord);
                gmsh::model::geo::addPoint(projcoord[0],projcoord[1],projcoord[2],0.,index+shift);

            }
        }
    }
    // INSERT LINES
    std::vector<int> curvesInSurface;
    {
        for(auto iter = fOutline.begin(); iter != fOutline.end(); iter++){
            int64_t iel = iter->first+shift;
            TPZGeoEl *gel = iter->second;
            if(gel->Dimension() != 1) continue;
            int64_t node0 = gel->NodeIndex(0)+shift;
            int64_t node1 = gel->NodeIndex(1)+shift;

            gmsh::model::geo::addLine(node0,node1,iel);
            gmsh::model::geo::mesh::setTransfiniteCurve(iel,2); // to reduce number of nodes created by GMsh
            bool lineIsInEdge = (std::find(outerLoop.begin(),outerLoop.end(),iel-shift) != outerLoop.end());
            if(lineIsInEdge == false){
                curvesInSurface.push_back(iel);
            }
        }
    }
    // CURVE LOOPS
        for(auto itr = outerLoop.begin(); itr != outerLoop.end(); itr++){
            *itr += shift;
        }
        int surfaceIndex = 0 + shift;
        std::vector<int> wiretags(facesInSurface.size()+1,-1);
        wiretags[0] = gmsh::model::geo::addCurveLoop(outerLoop,surfaceIndex);
    // Holes in surface
        for(int iface=0;iface<facesInSurface.size();iface++){
            TPZGeoEl* face = facesInSurface[iface];
            std::vector<int> loop;
            GetCurveLoop(face,loop,shift);
            wiretags[iface+1] = gmsh::model::geo::addCurveLoop(loop);
        }

    // SURFACE + HOLES
        gmsh::model::geo::addPlaneSurface(wiretags,surfaceIndex);
    
    // lines in surface
        // @comment GMsh requires synchronize before embedding geometric entities
        gmsh::model::geo::synchronize();
        gmsh::model::mesh::embed(1,curvesInSurface,2,surfaceIndex);
    // PHYSICAL GROUPS
        // physical curve
        int nlines = curvesInSurface.size() + outerLoop.size();
        if(curvesInSurface.size() > outerLoop.size()){
            curvesInSurface.reserve(nlines);
            curvesInSurface.insert(curvesInSurface.end(), outerLoop.begin(), outerLoop.end() );

            gmsh::model::addPhysicalGroup(1,curvesInSurface,DFNMaterial::Efracture);
        }else{
            outerLoop.reserve(nlines);
            outerLoop.insert(outerLoop.end(), curvesInSurface.begin(), curvesInSurface.end() );

            gmsh::model::addPhysicalGroup(1,outerLoop,DFNMaterial::Efracture);
        }
        // physical surface
        gmsh::model::addPhysicalGroup(2,{surfaceIndex},DFNMaterial::Efracture);

    // synchronize before meshing
        gmsh::model::geo::synchronize();
    // mesh
        gmsh::model::mesh::generate(2);
        // gmsh::model::mesh::optimize("Netgen");
    // write (for testing)
        // gmsh::write("testAPI.msh");
    // import meshed plane back into PZ geoMesh
        TPZVec<int64_t> newelements;
        ImportElementsFromGMSH(fdfnMesh->Mesh(),2,pointset,newelements);
    // close GMsh
    gmsh::model::remove();
	// gmsh::clear();
    // gmsh::finalize();
    
    InsertElementsInSurface(newelements);
    fdfnMesh->CreateSkeletonElements(1, DFNMaterial::Efracture);
}






namespace DFN{
	const float _2PI = 6.2831853071795865;
    template<typename Ttype>
    float DotProduct_f(TPZManVector<Ttype,3> &vec1, TPZManVector<Ttype,3> &vec2){
        int size1 = vec1.size();
        int size2 = vec2.size();
        if(size1 != size2){throw std::bad_exception();}
        float dot = 0.;
        for(int j=0; j<size1; j++){
            dot += vec1[j]*vec2[j];
        }
        return dot;
    }

    template<typename Ttype>
    float Norm_f(TPZManVector<Ttype, 3> &vec){
        float norm = 0.;
        for(int j=0, size=vec.size(); j<size; j++){
            norm += vec[j]*vec[j];
        }
        return std::sqrt(norm);
    }

    template<typename Ttype>
    TPZManVector<Ttype,3> CrossProduct_f(TPZManVector<Ttype,3> &vec1, TPZManVector<Ttype,3> &vec2){
        if(vec1.size() != 3){throw std::bad_exception();}
        if(vec2.size() != 3){throw std::bad_exception();}

        TPZManVector<REAL,3> result(3,0.);
        result[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
        result[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
        result[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
        return result;
    }
    template <class T, int NumExtAlloc1, int NumExtAlloc2>
    TPZManVector<T,3> operator-(TPZManVector<T,NumExtAlloc1>& v1,TPZManVector<T,NumExtAlloc2>& v2){
        int64_t size1 = v1.size();
        int64_t size2 = v2.size();
        if(size1 != size2) throw std::bad_exception();
        TPZManVector<T,3> result(size1);
        for(int64_t i = 0; i<size1; i++){
            result[i] = v1[i] - v2[i];
        }
        return result;
    }


    /**
     * @brief Returns the oriented dihedral angle between gel and neighbour
     * @note:1 Make sure neighbour is an actual neighbour, otherwise this method will spit nonsense
     * @note:2 Returned angle is in interval [0, 2pi)
     * @param sideorientation decides if method returns angle or 2*pi-angle, as a right-handed axis orientation
    */
    float DihedralAngle(TPZGeoElSide &gelside, TPZGeoElSide &neighbour, int sideorientation = 1){
        // Consistency checks
        if(gelside.Element()->Dimension() != 2)     DebugStop();
        if(gelside.Dimension() != 1)                DebugStop();
        if(neighbour.Element()->Dimension() !=2)    DebugStop();
        if(neighbour.Dimension() != 1)              DebugStop();
        if(!gelside.NeighbourExists(neighbour))     DebugStop();

        TPZGeoEl* gel = gelside.Element();
        TPZGeoMesh* gmesh = gel->Mesh();
        const int side = gelside.Side();
        TPZManVector<double,3> sharednode0(3,0);
        TPZManVector<double,3> sharednode1(3,0);
        gmesh->NodeVec()[gelside.SideNodeIndex(0)].GetCoordinates(sharednode0);
        gmesh->NodeVec()[gelside.SideNodeIndex(1)].GetCoordinates(sharednode1);

        TPZManVector<REAL,3> oppositenode_gel(3,0);
        TPZManVector<REAL,3> oppositenode_neig(3,0);
        gel->Node((gelside.Side()+2)%gel->NNodes()).GetCoordinates(oppositenode_gel);
        neighbour.Element()->Node((neighbour.Side()+2)%neighbour.Element()->NNodes()).GetCoordinates(oppositenode_neig);
        TPZManVector<REAL,3> tangentvec_gel = oppositenode_gel - sharednode0;
        TPZManVector<REAL,3> tangentvec_neig = oppositenode_neig - sharednode0;
        TPZManVector<REAL,3> tangentvec_edge(3);
        switch(sideorientation){
            case -1:{tangentvec_edge = sharednode1 - sharednode0; break;}
            case  1:{tangentvec_edge = sharednode0 - sharednode1; break;}
            default: DebugStop();
        }

        TPZManVector<REAL,3> normalvec_gel = CrossProduct_f(tangentvec_gel,tangentvec_edge);
        TPZManVector<REAL,3> normalvec_neig = CrossProduct_f(tangentvec_neig,tangentvec_edge);;
        TPZManVector<REAL,3> aux = CrossProduct_f(normalvec_neig,normalvec_gel);
        float x = Norm_f(tangentvec_edge)*DotProduct_f(normalvec_neig,normalvec_gel);
        float y = DotProduct_f(tangentvec_edge,aux);
        float angle = atan2f32(y,x);

        return (angle >= 0? angle : angle + _2PI);
    }

    /**
     * @brief Get a vector from node 0 to node 1 of a 1D side
    */
    void GetSideVector(TPZGeoElSide &gelside, TPZManVector<REAL,3>& vector){
        if(gelside.Dimension() != 1) DebugStop();
        int node0 = gelside.SideNodeLocIndex(0);
        int node1 = gelside.SideNodeLocIndex(1);

        TPZManVector<REAL,3> coord0(3,0);
        TPZManVector<REAL,3> coord1(3,0);
        gelside.Element()->Node(node0).GetCoordinates(coord0);
        gelside.Element()->Node(node1).GetCoordinates(coord1);

        vector = coord1 - coord0;
    }
}








void DFNFracture::GetSubPolygons(){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    std::map<int, TPZAutoPointer<std::vector<int>>> subpolygons_map; // @todo maybe change this to a vector of autopointers...
    // initialize a data structure to track which subpolygons have included which lines
    std::map<int64_t, std::pair<int,int>> LineTracker;
    for(auto iterator=fFaces.begin(); iterator!=fFaces.end(); iterator++){
        DFNFace* face = &iterator->second;
        int64_t line = face->LineInFace();
        if(line == -1) continue;
        LineTracker[line] = {0,0};
    }
    TPZManVector<REAL,3> frac_normal(3,0);
    fPolygon->GetNormal(frac_normal);
    for(auto iterator=fFaces.begin(); iterator!=fFaces.end(); iterator++){
        DFNFace* initial_face = &iterator->second;
        //Check if face intersection is an actual line or has been coalesced down to a point
        int64_t firstline = initial_face->LineInFace();
        int current_line = firstline;
        if(firstline == -1) continue;
        //Check if the line in face has already been attributed to all its possible subpolygons
        std::pair<int,int>* tracker = &LineTracker[firstline];
        int nloops = int(tracker->first>0) + int(tracker->second>0);
        if(nloops == 2) continue;
        int npolyhedra = this->fdfnMesh->FaceTracker[initial_face->Index()];
        // Decide a direction to follow
        int direction = (tracker->first>0?0:1);
        // Track if a direction was tried but failed
        bool DirectionFailed = false;
        while(nloops < npolyhedra){
            // TPZAutoPointer<std::vector<int>> subpolygon = new std::vector<int>;
            std::vector<int>* subpolygon = new std::vector<int>;

            // Find next side based in decided direction
            int nedges = initial_face->GeoEl()->NCornerNodes();
            int edge;
            for(edge = 0; edge<nedges; edge++){
                DFNRib* edgerib = initial_face->Rib(edge);
                if(!edgerib) continue;
                if(edgerib->GeoEl()->NodeIndex(0)==gmesh->Element(firstline)->NodeIndex(direction)) break;
                if(edgerib->GeoEl()->NodeIndex(1)==gmesh->Element(firstline)->NodeIndex(direction)) break;
                if(edgerib->IntersectionIndex()  ==gmesh->Element(firstline)->NodeIndex(direction)) break;
            }
            int nextside = edge + initial_face->GeoEl()->NCornerNodes();

            // Follow direction by getting the next neighbour with the smallest dihedral angle to close the subpolygon
            TPZGeoEl* current_face = initial_face->GeoEl();
            do{
                float angle = DFN::_2PI;
                if(current_line != -1){
                    (*subpolygon).push_back(current_line);
                    std::cout<<current_line<<std::endl;
                    if(current_line > 0)    LineTracker[abs(current_line)].first = 1;
                    else                    LineTracker[abs(current_line)].second = 1;
                }
                TPZGeoElSide gelside(current_face,nextside);
                TPZGeoElSide neig = gelside.Neighbour();
                // Check if gelside's side orientation agree's with fracture normal vector
                TPZManVector<REAL,3> gelside_vec(3,0);
                DFN::GetSideVector(gelside,gelside_vec);
                int sideorientation;
                if(DFN::DotProduct_f(gelside_vec,frac_normal)>0.){
                    sideorientation = -1;}
                else{
                    sideorientation = 1;
                }
                TPZGeoEl* next_face;
                int current_side = -1;
                for(/*void*/; neig!=gelside; neig = neig.Neighbour()){
                    if(neig.Element()->Dimension() != 2) continue;
                    if(!Face(neig.Element()->Index())) continue;
                    float temp_angle = DFN::DihedralAngle(gelside,neig,sideorientation);
                    if(temp_angle < angle){
                        angle = temp_angle;
                        next_face = neig.Element();
                        current_side = neig.Side();
                    }
                }
                if(angle > M_PIf32+gDFN_SmallNumber){
                    if(npolyhedra == 1 && !DirectionFailed) {
                        // Might be an unluckly bad oriented line in initial_face. So try going the other way before DebugStop.
                        DirectionFailed = true;
                        direction = (direction+1)%2;
                        current_line = -firstline;
                        delete &*subpolygon;
                        break;
                    }
                    PZError << "\nNon-convex regions shouldn't exist at this point\n" << __PRETTY_FUNCTION__ << std::endl;
                    DebugStop();
                }
                DirectionFailed = false;
                // Get line in face, its orientation and next side
                DFNFace* next_dfnface = Face(next_face->Index());
                current_line = next_dfnface->LineInFace();
                int orientation=0;
                for(int iedge=0; iedge < next_face->NCornerNodes(); iedge++){
                    if(iedge + next_face->NCornerNodes() == current_side) continue;
                    DFNRib* edge_rib = next_dfnface->Rib(iedge);
                    if(!edge_rib) continue;
                    nextside = iedge + next_face->NCornerNodes();
                    if(current_line != -1){
                        if(edge_rib->GeoEl()->NodeIndex(0)     ==gmesh->Element(current_line)->NodeIndex(1)) orientation = 1;
                        else if(edge_rib->GeoEl()->NodeIndex(1)==gmesh->Element(current_line)->NodeIndex(1)) orientation = 1;
                        else if(edge_rib->IntersectionIndex()  ==gmesh->Element(current_line)->NodeIndex(1)) orientation = 1;
                        else orientation = -1;
                    }else{orientation=1;}
                    break;
                }
                current_line = orientation*current_line;
                current_face = next_face;
            }while(current_face != initial_face->GeoEl());

            if(!DirectionFailed){
                nloops++;
                // subpolygons_map.insert({subpolygons_map.size(),subpolygon});
                subpolygons_map[subpolygons_map.size()] = subpolygon;
            }
        }
    }
    int i_polyg = 0;
    for(auto iterator=subpolygons_map.begin(); iterator!=subpolygons_map.end(); iterator++){
        std::vector<int> &polygonloop = *iterator->second;
        std::cout<<"\n\nSubPolygon #"<<i_polyg<<"\n";
        int size = polygonloop.size();
        for(int iline=0; iline<size; iline++){
            if(polygonloop[iline] > 0){
                std::cout<<" ";
            }
            std::cout<<"\n"<<polygonloop[iline];
        }
        i_polyg++;
    }
}







































void DFNFracture::ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension, std::set<int64_t> &oldnodes, TPZVec<int64_t> &newelements){
    // GMsh does not accept zero index entities
    const int shift = 1;

    // First, if GMsh has created new nodes, these should be inserted in PZGeoMesh
    // create a map <node,point>
    std::map<int,int> mapGMshToPZ;

    for(int64_t pznode : oldnodes){
		std::vector<size_t> node_identifiers;
        std::vector<double> coord;
        std::vector<double> parametricCoord;
        gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord,0,pznode+shift,true);
        int gmshnode = (int) node_identifiers[0];
		// insert with hint (since oldnodes is an already sorted set, these nodes will all go in the end)
        mapGMshToPZ.insert(mapGMshToPZ.end(),{gmshnode,pznode+shift});
	}

    // add new nodes into PZGeoMesh
    {
        // get all nodes from GMsh
            std::vector<size_t> node_identifiers;
            std::vector<double> coord;
            std::vector<double> parametricCoord;
            gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord);
        // iterate over node_identifiers
        int nnodes = node_identifiers.size();
        for(int i = 0; i < nnodes; i++){
            int gmshnode = node_identifiers[i];
            // check if it is contained in the map
            if(mapGMshToPZ.find(gmshnode) == mapGMshToPZ.end()){
                // New node -> add to PZGeoMesh
                int pznode = (int) gmesh->NodeVec().AllocateNewElement();
                TPZManVector<REAL,3> newnodeX(3);
                newnodeX[0] = coord[3*i];
                newnodeX[1] = coord[3*i+1];
                newnodeX[2] = coord[3*i+2];
                gmesh->NodeVec()[pznode].Initialize(newnodeX,*gmesh);
                // int pznode = (int) gmesh->NNodes();
                // gmesh->NodeVec().resize(pznode+1);
                // insert it in map
                mapGMshToPZ.insert({gmshnode,pznode+shift});
            }

        }
    }
    

    
    int64_t nels = gmesh->NElements();
    std::vector<std::pair<int, int> > dim_to_physical_groups;
    gmsh::model::getPhysicalGroups(dim_to_physical_groups,dimension);
   
    /// inserting the elements
    for (auto group: dim_to_physical_groups) {
       
        int dim = group.first;
        // only want elements of a given dimension
        if(dim != dimension) continue;
        int physical_identifier = group.second;
       
        std::vector< int > entities;
        gmsh::model::getEntitiesForPhysicalGroup(dim, physical_identifier, entities);

		for (auto tag: entities) {
		// std::cout<<"______________________test - tag = "<<tag;
           
            std::vector<int> group_element_types;
            std::vector<std::vector<std::size_t> > group_element_identifiers;
            std::vector<std::vector<std::size_t> > group_node_identifiers;
            gmsh::model::mesh::getElements(group_element_types,group_element_identifiers,group_node_identifiers, dim, tag);
            int n_types = group_element_types.size();
            for (int itype = 0; itype < n_types; itype++){
                int el_type = group_element_types[itype];
                int n_nodes = TPZGeoMeshBuilder::GetNumberofNodes(el_type);
                std::vector<int> node_identifiers(n_nodes);
                int n_elements = group_element_identifiers[itype].size();
                for (int iel = 0; iel < n_elements; iel++) {
                    int el_identifier = group_element_identifiers[itype][iel]-1+nels;
					// std::cout<<"\n"<<el_identifier<<"\n";

                    for (int inode = 0; inode < n_nodes; inode++) {
                        // node_identifiers[inode] = group_node_identifiers[itype][iel*n_nodes+inode];
                        // Use mapGMshToPZ to translate from GMsh node index to PZ nodeindex
                        node_identifiers[inode] = mapGMshToPZ[group_node_identifiers[itype][iel*n_nodes+inode]];
                    }
                    TPZGeoMeshBuilder::InsertElement(gmesh, physical_identifier, el_type, el_identifier, node_identifiers);
					int64_t ntest = gmesh->NElements();
					// std::cout<<"nelements = "<<ntest<<"\n";
                }
            }
        }
    }
    gmesh->BuildConnectivity();
}


















void DFNFracture::AssembleOutline(){
    // Build skeleton of 1D elements between new subelements
    fdfnMesh->CreateSkeletonElements(1);
    DFNFace *face = nullptr;
    TPZGeoEl *face_gel = nullptr;

    for(auto itr = fFaces.begin(); itr!=fFaces.end(); itr++){
        face = &itr->second;
        face_gel = face->GeoEl();
        int n_intersection_points = 0;
        for(int istate : face->StatusVec()){
            n_intersection_points += istate;
        }
        if(n_intersection_points < 2) continue;

        int nsides = face_gel->NSides();
        int nnodes = face_gel->NCornerNodes();
        TPZStack<int64_t> framenodes;
        for(int i=0; i<nsides; i++){
            if(face->StatusVec()[i]){
                if(i<nnodes){
                    framenodes.push_back(face_gel->NodeIndex(i));
                }else if(i<nsides-1){
                    framenodes.push_back(face->Rib(i-nnodes)->IntersectionIndex());
                }else if(i == nsides-1){
                    framenodes.push_back(face->IntersectionIndex());
                }
            }
        }
        TPZGeoEl *frame_el = nullptr;
        int64_t frame_el_index = -1;
        int nchildren = face_gel->NSubElements();
        if(nchildren == 0){nchildren++;}
        // queue all possible lines by checking 1D neighbours of children
		std::set<TPZGeoEl *> candidate_ribs;
		for(int ichild=0; ichild<nchildren; ichild++){
            TPZGeoEl* child;
            if(face_gel->HasSubElement()){
                child = face_gel->SubElement(ichild);
            }else{
                child = face_gel;
            }
			int nribs = child->NCornerNodes();
			for(int cside = nribs; cside < 2*nribs; cside++){
				TPZGeoElSide childside(child,cside);
				TPZGeoElSide neig = childside.Neighbour();
				for(/*void*/; neig != childside; neig = neig.Neighbour()){
					if(neig.Element()->Dimension() != 1) continue;
					candidate_ribs.insert(neig.Element());
                    break;
				}
			}
		}
        bool enters_Q = false;
        int nframenodes = framenodes.size();
        int counter = 0;
        for(auto candidate : candidate_ribs){
            for(int inode=0; inode<nframenodes; inode++){
                if(framenodes[inode]==candidate->NodeIndex(0)){
                    for(int nextnode=0; nextnode<nframenodes; nextnode++){
                        if(framenodes[nextnode]==candidate->NodeIndex(1)){
                            enters_Q = true;
                            break;
                        }
                    }
                }
                if(enters_Q) break;
            }
            if(!enters_Q) continue;
            frame_el = candidate;
            frame_el_index = candidate->Index();
            fOutline.insert({frame_el_index,frame_el});
            frame_el->SetMaterialId(DFNMaterial::Efracture);
            enters_Q = false;
            counter++;
            if(counter==2) break; //2 ribs is already an exception... 3 should never happen
            // if(counter==framenodes.size()-1) break;
        }
    }
}