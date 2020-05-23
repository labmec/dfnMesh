/*! 
 *  @brief     Compares a geomesh with fracture plane to refine intersections, and mesh the surface of the fracture.
 *  @details   Intersection search is performed after creation of skeleton mesh
 *  with DFNFractureMesh::CreateSkeletonElements. Fracture plane should
 *  be a DFNFracPlane.
 *  @authors   Pedro Lima
 *  @date      2019
 */

#include "DFNFractureMesh.h"
#include <math.h>
#include <cstdio>
// #include <unordered_set>
#include <algorithm>
#include "TPZRefPatternDataBase.h"
#include "TPZGeoMeshBuilder.h"

// Empty Constructor
DFNFractureMesh::DFNFractureMesh(){
}

// Constructor with corner points, a geomesh and material ID
DFNFractureMesh::DFNFractureMesh(DFNFracPlane &FracPlane, TPZGeoMesh *gmesh, int matID){
    fFracplane = FracPlane;
    fGMesh = gmesh;
    fSurfaceMaterial = matID;

    // @ToDo Maybe implement check for previously created skeleton
    // Create skeleton elements
    int materialSkeleton = 4;
    CreateSkeletonElements(2, materialSkeleton);
    CreateSkeletonElements(1, materialSkeleton);

    // Set corner nodes of fracture into mesh
    TPZManVector<int64_t,4> nodeindices = fFracplane.SetPointsInGeomesh(fGMesh);
    // int64_t index;
    // gmesh->CreateGeoElement(EQuadrilateral,nodeindices,40,index);
}

// Copy constructor
DFNFractureMesh::DFNFractureMesh(const DFNFractureMesh &copy){
    fGMesh = copy.fGMesh;
    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
    fMidFaces = copy.fMidFaces;
	fEndFaces = copy.fEndFaces;
    fFracplane = copy.fFracplane;
    fSurfaceMaterial = copy.fSurfaceMaterial;
    fTransitionMaterial = copy.fTransitionMaterial;
}

// Assignment operator
DFNFractureMesh &DFNFractureMesh::operator=(const DFNFractureMesh &copy){
    fGMesh = copy.fGMesh;
    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
    fMidFaces = copy.fMidFaces;
	fEndFaces = copy.fEndFaces;
    fFracplane = copy.fFracplane;
    fSurfaceMaterial = copy.fSurfaceMaterial;
    fTransitionMaterial = copy.fTransitionMaterial;
    return *this;
}

/**
 * @brief Get the fracture plane
 * @return The plane corner coordinates
 */

DFNFracPlane &DFNFractureMesh::GetPlane() {
    return fFracplane;
}

/**
 * @brief Set the tolerance for the distance between a point-plane
 * @param Tolerance
 */

void DFNFractureMesh::SetTolerance(REAL tolerance){
    fTolerance = tolerance;
}

/**
 * @brief Get the tolerance
 * @return The tolerance
 */

REAL DFNFractureMesh::GetTolerance() const{
    return fTolerance;
}









/**
 * @brief Check if the neighbour has equal dimension
 * @param geliside GeoElement side
 */

bool DFNFractureMesh::HasEqualDimensionNeighbour(TPZGeoElSide &gelside){
    
    int dimension = gelside.Dimension();

    if (gelside.Element()->Dimension() == dimension){
        return true;
    }

    TPZGeoElSide neighbour = gelside.Neighbour();

    while (neighbour != gelside){
        if (neighbour.Element()->Dimension()==dimension){
            return true;
        }
        neighbour = neighbour.Neighbour();
    }
    return false;    
}










/**
  * @brief Creates the skeleton mesh
  * @param dimension
  * @param matid material ID number for skeleton elements
  */

void DFNFractureMesh::CreateSkeletonElements(int dimension, int matid)
{
    int nel = fGMesh->NElements();
    for (int iel = 0; iel < nel; iel++)
    {
        TPZGeoEl *gel = fGMesh->Element(iel);
        if(!gel) continue;
        // Elements can't have a skeleton of higher dimension than itself
        if(gel->Dimension() <= dimension) continue;
        
        int nsides = gel->NSides();
        int ncorners = gel->NCornerNodes();
        // iterating from higher-dimensional sides to lower-dimensional should narrow the search
        for (int iside = nsides-2; iside >= ncorners; iside--)
        {
            TPZGeoElSide gelside = gel->Neighbour(iside);

            if (gelside.Dimension() != dimension){continue;}
            bool haskel = HasEqualDimensionNeighbour(gelside);
            if (haskel == false)
            {
                if(matid == -1) matid = gel->MaterialId();
                // TPZGeoElBC(gelside, matid);
                if(gel->MaterialId() >= fSurfaceMaterial) TPZGeoElBC(gelside, fSurfaceMaterial+6);
                else TPZGeoElBC(gelside, matid);
            }
        }
    }
}











/**
 * @brief Sets the rib idexes in the map
 * @param rib Rib to be set
 */

void DFNFractureMesh::AddRib(DFNRibs rib){
    int index= rib.ElementIndex();
    fRibs.emplace(index,rib);
}

/**
 * @brief Add cut faces using indexes
 * @param Face to be set
 */

bool DFNFractureMesh::AddMidFace(DFNFace &face){
    // iterate over ribs to find intersected ones
    int nribscut = 0;
    TPZManVector<int64_t,2> CutRibsIndex(2);
    TPZVec<int64_t> rib_index = face.GetRibs();
    int nribs = rib_index.size();
    for (int irib = 0; irib < nribs; irib++)
    {
        DFNRibs *ribtest = &fRibs[rib_index[irib]];
        if(ribtest->IsCut() == true){
            CutRibsIndex[nribscut]=rib_index[irib];
            nribscut++;
        }
    }
    // std::cout<<"first rib: "<<CutRibsIndex[0]<<std::endl;
    // std::cout<<"second rib: "<<CutRibsIndex[1]<<std::endl;

    // Connect intersection points
    TPZVec<int64_t> ipoints(2);
    ipoints[0] = Rib(CutRibsIndex[0])->IntersectionIndex();
    ipoints[1] = Rib(CutRibsIndex[1])->IntersectionIndex();
    int64_t nels = fGMesh->NElements();
    this->fSurfEl[nels] = fGMesh->CreateGeoElement(EOned, ipoints, fSurfaceMaterial+6, nels);
                                                                        // +6 for graphical debugging
    
    int index= face.ElementIndex();
    fMidFaces[index]=face;
    return true;
}

/**
 * @brief Add faces that are cut at the edges of fracture (using indexes)
 * @param Face to be set
 */

bool DFNFractureMesh::AddEndFace(DFNFace &face){
    // Create geometric element for intersection node beetween EndFace and fracture edge
        TPZManVector<REAL,3> coords(3);
        bool ipoint_exists = FindEndFracturePoint(face,coords);
        if(ipoint_exists == false) return false;
        TPZManVector<int64_t> ipoints(2);
        ipoints[0] = fGMesh->NodeVec().AllocateNewElement();
        fGMesh->NodeVec()[ipoints[0]].Initialize(coords, *fGMesh);
        face.SetIntersectionIndex(ipoints[0]);
    //iterate over ribs to connect intersection points
    TPZVec<int64_t> rib_index = face.GetRibs();
    int nribs = rib_index.size();
    for (int irib = 0; irib < nribs; irib++)
    {
        DFNRibs *ribtest = &fRibs[rib_index[irib]];
        if(ribtest->IsCut() == true){
            // std::cout<<"single rib cut: "<<rib_index[irib]<<std::endl;
            // Connect intersection points
                int64_t nels = fGMesh->NElements();
                ipoints[1] = Rib(rib_index[irib])->IntersectionIndex();
                this->fSurfEl[nels] = fGMesh->CreateGeoElement(EOned, ipoints, fSurfaceMaterial+6, nels);
                break;                                                                  // +6 for graphical debugging
        }
    }
    int index= face.ElementIndex();
    fEndFaces[index]=face;
    return true;
}









DFNFace * DFNFractureMesh::Face(int64_t index){
    auto candidate = fMidFaces.find(index);
    if(candidate != fMidFaces.end()){
        return &candidate->second;
    }
    candidate = fEndFaces.find(index);
    if(candidate != fEndFaces.end()){
        return &candidate->second;
    }
    return NULL;
}




/**
 * @brief Create cut surfaces
 * @param Material id
 * @note This method was getting too long. Moved some of it into AddEndFace and AddMidFace.
 */

void DFNFractureMesh::SplitFaces(int matID){

    fGMesh->BuildConnectivity();

    // iterate over all 2D elements and check their 1D neighbours for intersections
    int64_t nel = fGMesh->NElements();
    for(int iel = 0; iel<nel; iel++){
        TPZGeoEl *gel = fGMesh->Element(iel);
        if(!gel) continue;
        if (gel->Dimension() != 2){continue;}
        if(gel->HasSubElement()) continue;
        //fSurfaceMaterial is MaterialID for fracture plane
        if(gel->MaterialId()==fSurfaceMaterial){continue;}
        
        int nribscut =0;
        int nedges = gel->NCornerNodes();

        // vector with status for each node and rib of face
        TPZManVector<bool,8> sidestatus(nedges*2,false);
        // TPZManVector<int64_t,2> CutRibsIndex(2);
        // vector with indices of father ribs that outline the face
        TPZManVector<int64_t,4> rib_index(nedges,-1);

        // iterate over ribs to check for intersection
        for(int iside = 0; iside < nedges; iside++){
            TPZGeoElSide gelside(gel,iside+nedges);
            TPZGeoElSide neig;
            // TPZGeoElSide neig = gelside.Neighbour();
            // while(neig.Element()->Dimension()!=1 || !neig.Element()->HasSubElement()){ 
            //     neig=neig.Neighbour();
            // }
            // rib_index[iside] = neig.Element()->Index();
            for(neig = gelside.Neighbour(); neig!= gelside; neig = neig.Neighbour()){
                if(neig.Element()->Dimension() == 1 && neig.Element()->HasSubElement()) {
                    rib_index[iside] = neig.Element()->Index();
                    break;
                }
                if(neig.Element()->Dimension() == 1) rib_index[iside] = neig.Element()->Index();
            }
            DFNRibs *ribtest = &fRibs[rib_index[iside]];
            if(ribtest->IsCut()==true){
            //check if ribtest was divided into two ribs, or a rib and a point
                // get node where ribtest is cut
                int64_t cutnode = ribtest->IntersectionIndex();
                TPZVec<int64_t> ribtestNodes(2);
                // fGMesh->Element(ribtest->ElementIndex())->GetNodeIndices(ribtestNodes);
                ribtestNodes[0] = gel->SideNodeIndex(iside+nedges,0);
                ribtestNodes[1] = gel->SideNodeIndex(iside+nedges,1);
                // check if intersection is a corner node and assign status accordingly
                if(cutnode == ribtestNodes[0]){
                    sidestatus[iside] = true;
                }
                else if(cutnode == ribtestNodes[1]){
                    sidestatus[(iside+1)%nedges] = true;
                }
                else{
                    sidestatus[iside+nedges] = true;
                }
                
                nribscut++;
            }
        }

        // if there are ribs cut, create a face object
        if(nribscut == 0){continue;}

        int edgescut = 0, cornerscut = 0;
        for(int i = 0; i<nedges; i++){
            if(sidestatus[i]) cornerscut++;
            if(sidestatus[i+nedges]) edgescut++;
        }
        if(edgescut == 0 && cornerscut == 1) nribscut = 1;
        // Create DFNFace
        DFNFace face(iel, true);
        face.SetRibs(rib_index); 
        // During development, elements at fracture surface have material id over fSurfaceMaterial
        // if(gel->MaterialId() != fSurfaceMaterial) {gel->SetMaterialId(matID);}
        face.SetStatus(sidestatus);
        face.SetFractureMesh(this);
        // if(nribscut == 1) {gel->SetMaterialId(matID+17);} // this is here for graphical debugging only... comment it on release

        // Add face to map
        bool face_ok = false;
        switch (nribscut){
            case  2: face_ok = AddMidFace(face);break;
            case  1: face_ok = AddEndFace(face);break;
            default: std::cout<<"\n\n"<<__PRETTY_FUNCTION__<<"\n Face # "<<iel<<"\nNo more than 2 ribs should've been cut\n\n";DebugStop();
        }
        if(face_ok){ 
            if(gel->MaterialId() < fSurfaceMaterial) {gel->SetMaterialId(matID);}
            Face(iel)->DivideSurface(gel->MaterialId());
        }

        // //Print result
		// 	std::ofstream out1("./TestSurfaces.vtk");
		// 	TPZVTKGeoMesh::PrintGMeshVTK(fGMesh, out1, true);
    }
    fGMesh->BuildConnectivity();
    CreateSkeletonElements(1, fTransitionMaterial+3);
}














bool DFNFractureMesh::FindEndFracturePoint(DFNFace &face,TPZVec<REAL> &ipoint){
    // Convert TPZGeoEl into DFNFracPlane
    TPZGeoEl *gelface = fGMesh->Element(face.ElementIndex());
    TPZFMatrix<REAL> corners(3,4);
    int n;
    // check if face is quadrilateral
    if(gelface->Type() == ETriangle){
        n = 1;
    }else{ //gelface->Type() == EQuadrilateral
        n = 2;
    }

    gelface->NodesCoordinates(corners);
    for(int iplane = 0; iplane < n; iplane++){
        // divide quadrilaterals into 2 triangles in order to account for sets of points which are not coplanar
        TPZFMatrix<REAL> subcorners(3,3,0);
        if(n>1){
            for(int j = 0; j<3; j++){
                subcorners(j,0) = corners(j,2*iplane);
                subcorners(j,1) = corners(j,2*iplane+1);
                subcorners(j,2) = corners(j,(2*iplane+3)%4);
            }
        }else{
            subcorners = corners;
        }
        DFNFracPlane faceplane(subcorners);
        // Check fFracplane's ribs for intersection with faceplane
        int nribs = fFracplane.GetCornersX().Cols();
        for(int irib = 0; irib < nribs; irib++){
            TPZManVector<REAL,3> p1(3);
            TPZManVector<REAL,3> p2(3);
            for(int i = 0; i<3; i++){
                p1[i] = fFracplane.GetCornersX()(i, irib);
                p2[i] = fFracplane.GetCornersX()(i, (irib+1)%nribs);
            }
            if(faceplane.Check_rib(p1, p2)){
                ipoint = faceplane.CalculateIntersection(p1,p2);
                return true;
            }
        }
    }
    // std::cout<<"\n DFNFractureMesh::FindEndFracturePoint\n";
    // std::cout << "\nFailed to find intersection point in end-fracture face index: " << face.ElementIndex() << std::endl;
    // DebugStop();
    return false;
}







void DFNFractureMesh::SplitRibs(int matID){
    // if(!gRefDBase.GetUniformRefPattern(EOned)){
    //     gRefDBase.InitializeUniformRefPattern(EOned);
    // }
    //search gmesh for cut ribs
    int64_t Nels = fGMesh->NElements();
    for (int iel = 0; iel < Nels; iel++){
        TPZGeoEl *gel = fGMesh->Element(iel);
        int nSides = gel->NSides();

        //skip all elements that aren't ribs
        if (gel->Dimension() != 1){continue;}
        // skip all elements that have been cut by a previous fracture
        if(gel->HasSubElement()) continue;

        // Get rib's vertices
        int64_t p1 = gel->SideNodeIndex(2, 0);
        int64_t p2 = gel->SideNodeIndex(2, 1);
        TPZManVector<REAL,3> pp1(3);
        TPZManVector<REAL,3> pp2(3);
        fGMesh->NodeVec()[p1].GetCoordinates(pp1);
        fGMesh->NodeVec()[p2].GetCoordinates(pp2);

        // Check rib
        bool resul = fFracplane.Check_rib(pp1, pp2);

        // Split rib
        if (resul == true){
            // std::cout<<"\nRib # "<<iel;
            DFNRibs rib(iel, true);
            AddRib(rib);
            TPZVec<REAL> ipoint = fFracplane.CalculateIntersection(pp1, pp2);
            // During development, elements at fracture surface have material id bigger than fSurfaceMaterial
            // if(gel->MaterialId() != fSurfaceMaterial) {gel->SetMaterialId(matID);}
            if(gel->MaterialId() < fSurfaceMaterial) {gel->SetMaterialId(matID);}
            Rib(iel)->DivideRib(fGMesh, ipoint, -1);
            if(gel->MaterialId() < fSurfaceMaterial) {gel->SetMaterialId(fTransitionMaterial);}

            // std::cout<<"Element: "<<iel<<" Side: "<<side<<" Rib status: "<<resul<<" Fracture : 1"<<"\n";
        }
    }
}








void DFNFractureMesh::SplitFractureEdge(std::list<int> &fracEdgeLoop){ 
    
    // Compute fracture edges' lenghts
    int nedges = fFracplane.GetCornersX().Cols();
    TPZVec<REAL> edgelength(nedges,0);
    Matrix fraccorners(3,nedges);
    fraccorners = fFracplane.GetCornersX();
    for(int i = 0; i < nedges; i++){
        edgelength[i] = sqrtl(pow(fraccorners(0,i)-fraccorners(0,(i+1)%nedges),2)
                              +pow(fraccorners(1,i)-fraccorners(1,(i+1)%nedges),2)
                              +pow(fraccorners(2,i)-fraccorners(2,(i+1)%nedges),2));
    }
    
	// vector of pointers to maps
	TPZVec<std::map<REAL, int64_t>* > edgemap(nedges);
    for(int i = 0; i < nedges; i++){
        edgemap[i] = new std::map<REAL, int64_t>;
    }
    
    // iterate over all endfaces and map it to the fracture-edge that cuts it
    for (auto it = fEndFaces.begin(); it != fEndFaces.end(); it++){
        // get intersection node index and coordinates
        DFNFace *iface = &it->second;
        int64_t ipointindex = iface->IntersectionIndex();
        TPZVec<REAL> ipointcoord(3);
        fGMesh->NodeVec()[ipointindex].GetCoordinates(ipointcoord);

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
			if(dist<fTolerance){
				// compute local 1D coordinate (alpha)
				REAL norm = sqrtl(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
				REAL alpha = norm/edgelength[iedge];
				if(alpha > 1+fTolerance){std::cout<<"\nProblem with alpha\n";DebugStop();}
				// map intersection indexes from smallest alpha to biggest
                // map <alpha, index>
				edgemap[iedge]->insert({alpha, ipointindex});
			}
        }
    }
    
    // GMsh does not accept zero index entities
    const int shift = 1;

    bool warning_message_Q = true;
    //Once intersections on fracture-edges have been properly ordered and mapped by edge
	//iterate over edges to split them
	for (int iedge = 0; iedge < nedges; iedge++)
	{
		int64_t nels = fGMesh->NElements();
		TPZManVector<int64_t,2> inodes(2);     //index of nodes to be connected
        int icorner = iedge; //for readability
        
        if(edgemap[iedge]->size() == 0){
            #ifdef PZDEBUG
                if(warning_message_Q){
                    std::cout<<"\n Warning: Is there an edge of a fracture that doesn't cut any element? \n\n";
                    warning_message_Q = false;
                }
            #endif //PZDEBUG
            // DebugStop();
            inodes[0] = fFracplane.CornerIndex(icorner);
            inodes[1] = fFracplane.CornerIndex((icorner+1)%nedges);
            this->fSurfEl[nels] = fGMesh->CreateGeoElement(EOned, inodes, fSurfaceMaterial+6, nels);
            fracEdgeLoop.push_back((int) nels + shift);
            continue;
        }


		// connect first end-face intersection to iedge's first node
		auto it = edgemap[iedge]->begin();
        inodes[0] = fFracplane.CornerIndex(icorner);
		inodes[1] = it->second;
		this->fSurfEl[nels] = fGMesh->CreateGeoElement(EOned, inodes, fSurfaceMaterial+6, nels);
        fracEdgeLoop.push_back((int) nels + shift);
        
        // iterate over iedge's map
        while(++it != edgemap[iedge]->end()){
            nels++;
            inodes[0] = inodes[1];
            inodes[1] = it->second;
            this->fSurfEl[nels] = fGMesh->CreateGeoElement(EOned, inodes, fSurfaceMaterial+6, nels);
            fracEdgeLoop.push_back((int) nels + shift);
        }

		// connect last end-intersection to edge last node
		nels++;
		inodes[0] = inodes[1];
        inodes[1] = fFracplane.CornerIndex((icorner+1)%nedges);
		this->fSurfEl[nels] = fGMesh->CreateGeoElement(EOned, inodes, fSurfaceMaterial+6, nels);
        fracEdgeLoop.push_back((int) nels + shift);
    }
	
fGMesh->BuildConnectivity();  
}












void DFNFractureMesh::SplitFracturePlane(){
    // GMsh does not accept zero index entities
    const int shift = 1;
    // First construct the edges of the fracture surface
    std::list<int> fracEdgeLoop;
    SplitFractureEdge(fracEdgeLoop);

    // initialize GMsh
    gmsh::initialize();
    gmsh::option::setNumber("Mesh.Algorithm", 5); // (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
    // INSERT POINTS
        // iterate over ribs and get their ipoints
        std::set<int64_t> pointset;
        
        for(auto itr = fRibs.begin(); itr != fRibs.end(); itr++){
            DFNRibs *irib = &itr->second;
            if (irib->IsCut() == false){continue;}
            int64_t inode = irib->IntersectionIndex();
            pointset.insert(inode);
        }
        for(auto point : pointset){
            TPZManVector<REAL, 3> coord(3);
            fGMesh->NodeVec()[point].GetCoordinates(coord);

            gmsh::model::geo::addPoint(coord[0],coord[1],coord[2],0.,point+shift);

        }
        // iterate over endFaces and get their ipoints
        for(auto itr = fEndFaces.begin(); itr != fEndFaces.end(); itr++){
            DFNFace *iface = &itr->second;
            if (iface->IsCut() == false){continue;}
            int64_t inode = iface->IntersectionIndex();
            TPZManVector<REAL, 3> coord(3);
            fGMesh->NodeVec()[inode].GetCoordinates(coord);

            gmsh::model::geo::addPoint(coord[0],coord[1],coord[2],0.,inode+shift);
        }
        // Corners of fracture plane
        {
            int ncorners = fFracplane.GetCornersX().Cols();
            for(int i = 0; i<ncorners; i++){
                int64_t inode = fFracplane.CornerIndex(i);
                TPZManVector<REAL, 3> coord(3);
                fGMesh->NodeVec()[inode].GetCoordinates(coord);

                gmsh::model::geo::addPoint(coord[0],coord[1],coord[2],0.,inode+shift);
            }
        }
    // INSERT LINES
        // move loop list into a vector
        std::vector<int> edgeloopvector{std::make_move_iterator(std::begin(fracEdgeLoop)),
                                        std::make_move_iterator(std::end(fracEdgeLoop))};
        std::vector<int> curvesInSurface;
        for(auto iter = fSurfEl.begin(); iter != fSurfEl.end(); iter++){
            int64_t iel = iter->first+shift;
            TPZGeoEl *gel = iter->second;
            if(gel->Dimension() != 1) continue;
            int64_t node0 = gel->NodeIndex(0)+shift;
            int64_t node1 = gel->NodeIndex(1)+shift;

            gmsh::model::geo::addLine(node0,node1,iel);
            gmsh::model::geo::mesh::setTransfiniteCurve(iel,2); // to reduce number of nodes created by GMsh
            bool lineIsInEdge = (std::find(edgeloopvector.begin(),edgeloopvector.end(),iel) != edgeloopvector.end());
            if(lineIsInEdge == false){
                curvesInSurface.push_back(iel);
            }
        }
    
    // CURVE LOOP
        int surfaceIndex = 1;
        gmsh::model::geo::addCurveLoop(edgeloopvector,surfaceIndex);

    // SURFACE
        gmsh::model::geo::addPlaneSurface({surfaceIndex},surfaceIndex);
    
    // lines in surface
        // @comment GMsh requires synchronize before embedding geometric entities
        gmsh::model::geo::synchronize();
        gmsh::model::mesh::embed(1,curvesInSurface,2,surfaceIndex);
    // PHYSICAL GROUPS
        // physical curve
        int nlines = curvesInSurface.size() + edgeloopvector.size();
        if(curvesInSurface.size() > edgeloopvector.size()){
            curvesInSurface.reserve(nlines);
            curvesInSurface.insert(curvesInSurface.end(), edgeloopvector.begin(), edgeloopvector.end() );

            gmsh::model::addPhysicalGroup(1,curvesInSurface,fSurfaceMaterial+6);
        }else{
            edgeloopvector.reserve(nlines);
            edgeloopvector.insert(edgeloopvector.end(), curvesInSurface.begin(), curvesInSurface.end() );

            gmsh::model::addPhysicalGroup(1,edgeloopvector,fSurfaceMaterial+6);
        }
        // physical surface
        gmsh::model::addPhysicalGroup(2,{surfaceIndex},fSurfaceMaterial+7);

    // synchronize before meshing
        gmsh::model::geo::synchronize();
    // mesh
        gmsh::model::mesh::generate(2);
        // gmsh::model::mesh::optimize("Netgen");
    // write (for testing)
        // gmsh::write("testAPI.msh");
    // import meshed plane back into PZ geoMesh
        ImportElementsFromGMSH(fGMesh,2);
    // close GMsh
    gmsh::finalize();
    
    CreateSkeletonElements(1, fSurfaceMaterial+6);
}


void DFNFractureMesh::AddVolume(DFNVolume volume){
    int index = volume.ElementIndex();
    fVolumes[index] = volume;
}



/// Sets material for elements at surface of fracture (40 is default)
void DFNFractureMesh::SetSurfaceMaterial(int matID = 40){
    int64_t nels = fGMesh->NElements();
    for (int64_t iel = 0; iel < nels; iel++)
    {
        if(fGMesh->Element(iel)->MaterialId() == fSurfaceMaterial){
            fGMesh->Element(iel)->SetMaterialId(matID);
        }
    }
    fSurfaceMaterial = matID;
}



































void DFNFractureMesh::ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension){
    // GMsh does not accept zero index entities
    const int shift = 1;

    // First, if GMsh has created new nodes, these should be inserted in PZGeoMesh
    // create a map <node,point>
    std::map<int,int> mapGMshToPZ;

    // iterate over ribs and get intersection nodes
    for(auto itr = fRibs.begin(); itr != fRibs.end(); itr++){
        DFNRibs *irib = &itr->second;
        if (irib->IsCut() == false) continue;
        int pznode = (int) irib->IntersectionIndex() +shift;

        std::vector<size_t> node_identifiers;
        std::vector<double> coord;
        std::vector<double> parametricCoord;
        gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord,0,pznode,true);
        int gmshnode = (int) node_identifiers[0];
        mapGMshToPZ.insert({gmshnode,pznode});
    }
    // iterate over endFaces and get their ipoints
    for(auto itr = fEndFaces.begin(); itr != fEndFaces.end(); itr++){
        DFNFace *iface = &itr->second;
        if (iface->IsCut() == false) continue;
        int pznode = (int) iface->IntersectionIndex() +shift;
        
        std::vector<size_t> node_identifiers;
        std::vector<double> coord;
        std::vector<double> parametricCoord;
        gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord,0,pznode,true);
        int gmshnode = (int) node_identifiers[0];
        mapGMshToPZ.insert({gmshnode,pznode});
    }
    // Corners of fracture plane
    {
        int ncorners = fFracplane.GetCornersX().Cols();
        for(int i = 0; i<ncorners; i++){
            int pznode = (int) fFracplane.CornerIndex(i) +shift;
            
            std::vector<size_t> node_identifiers;
            std::vector<double> coord;
            std::vector<double> parametricCoord;
            gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord,0,pznode,true);
            int gmshnode = (int) node_identifiers[0];
            mapGMshToPZ.insert({gmshnode,pznode});
        }
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
    
    // remember to use mapGMshToPZ to translate from GMsh node index to PZ nodeindex




    
    int64_t nels = gmesh->NElements();
    std::vector<std::pair<int, int> > dim_to_physical_groups;
    gmsh::model::getPhysicalGroups(dim_to_physical_groups);
   
    std::vector<std::pair<int, int> > entities_0d;
    std::vector<std::pair<int, int> > entities_1d;
    std::vector<std::pair<int, int> > entities_2d;
    std::vector<std::pair<int, int> > entities_3d;
    gmsh::model::getEntities(entities_0d,0);
    gmsh::model::getEntities(entities_1d,1);
    gmsh::model::getEntities(entities_2d,2);
    gmsh::model::getEntities(entities_3d,3);
   
    /// inserting the elements
    for (auto group: dim_to_physical_groups) {
       
        int dim = group.first;
        // only want elements of a given dimension
        if(dim != dimension) continue;
        int physical_identifier = group.second;
       
        std::vector< int > entities;
        gmsh::model::getEntitiesForPhysicalGroup(dim, physical_identifier, entities);
       
        for (auto tag: entities) {
           
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
                    for (int inode = 0; inode < n_nodes; inode++) {
                        // node_identifiers[inode] = group_node_identifiers[itype][iel*n_nodes+inode];
                        // Use mapGMshToPZ to translate from GMsh node index to PZ nodeindex
                        node_identifiers[inode] = mapGMshToPZ[group_node_identifiers[itype][iel*n_nodes+inode]];
                    }
                    TPZGeoMeshBuilder::InsertElement(gmesh, physical_identifier, el_type, el_identifier, node_identifiers);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
}