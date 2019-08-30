/*! 
 *  @brief     Compares a geomesh with fracture plane to find intersections.
 *  @details   Intersection search is performed after creation of skeleton
 *  elements with TRSFractureMesh::CreateSkeletonElements. Fracture plane should
 *  be a TRSFracPlane.
 *  @authors   Pedro Lima
 *  @authors   Jorge Ordoñez
 *  @date      2018-2019
 */

#include "TRSFractureMesh.h"
#include <math.h>
#include <cstdio>
#include <unordered_set>

// Empty Constructor
TRSFractureMesh::TRSFractureMesh(){
}

// Constructor with corner points, a geomesh and material ID
TRSFractureMesh::TRSFractureMesh(TRSFracPlane &FracPlane, TPZGeoMesh *gmesh, int matID){
    fFracplane = FracPlane;
    fGMesh = gmesh;
    fSurfaceMaterial = matID;

    // Maybe implement check for previously created skeleton
    // Create skeleton elements
    int materialSkeleton = 4;
    CreateSkeletonElements(2, materialSkeleton);
    CreateSkeletonElements(1, materialSkeleton);

    // Create FracPlane's geometric element into mesh
    // fFracplaneindex = fFracplane.CreateElement(fGMesh);
    fFracplane.SetPointsInGeomesh(fGMesh, fSurfaceMaterial);
    
}

// Copy constructor
TRSFractureMesh::TRSFractureMesh(const TRSFractureMesh &copy){
    fGMesh = copy.fGMesh;
    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
    fMidFaces = copy.fMidFaces;
	fEndFaces = copy.fEndFaces;
    fFracplane = copy.fFracplane;
	fFracplaneindex = copy.fFracplaneindex;
    fSurfaceMaterial = copy.fSurfaceMaterial;
    fTransitionMaterial = copy.fTransitionMaterial;
}

// Assignment operator
TRSFractureMesh &TRSFractureMesh::operator=(const TRSFractureMesh &copy){
    fGMesh = copy.fGMesh;
    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
    fMidFaces = copy.fMidFaces;
	fEndFaces = copy.fEndFaces;
    fFracplane = copy.fFracplane;
	fFracplaneindex = copy.fFracplaneindex;
    fSurfaceMaterial = copy.fSurfaceMaterial;
    fTransitionMaterial = copy.fTransitionMaterial;
    return *this;
}

/**
 * @brief Get the fracture plane
 * @return The plane corner coordinates
 */

TRSFracPlane TRSFractureMesh::GetPlane() const{
    return fFracplane;
}

/**
 * @brief Set the tolerance for the distance between a point-plane
 * @param Tolerance
 */

void TRSFractureMesh::SetTolerance(REAL tolerance){
    fTolerance = tolerance;
}

/**
 * @brief Get the tolerance
 * @return The tolerance
 */

REAL TRSFractureMesh::GetTolerance() const{
    return fTolerance;
}









/**
 * @brief Check if the neighbour has equal dimension
 * @param Geo element side
 */

bool TRSFractureMesh::HasEqualDimensionNeighbour(TPZGeoElSide &gelside){
    
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
  * @param Dimension
  * @param Material ID number
  */

void TRSFractureMesh::CreateSkeletonElements(int dimension, int matid)
{
    int nel = fGMesh->NElements();
    for (int iel = 0; iel < nel; iel++)
    {
        TPZGeoEl *gel = fGMesh->Element(iel);

        //40 is material id for fracture plane, and it should only have an 1D skeleton.
        if(gel->MaterialId() == fSurfaceMaterial && dimension == 2){continue;}

        int nsides = gel->NSides();
        for (int iside = 0; iside < nsides; iside++)
        {
            TPZGeoElSide gelside = gel->Neighbour(iside);

            if (gelside.Dimension() != dimension){continue;}
            bool haskel = HasEqualDimensionNeighbour(gelside);
            if (haskel == false)
            {
                TPZGeoElBC(gelside, matid);
            }
        }
    }
}











/**
 * @brief Sets the rib idexes in the map
 * @param Ribs to be set
 */

void TRSFractureMesh::AddRib(TRSRibs rib){
    int index= rib.ElementIndex();
    fRibs[index]=rib;
}

/**
 * @brief Add cut faces using indexes
 * @param Face to be set
 */

void TRSFractureMesh::AddMidFace(TRSFace &face){
    // iterate over ribs to find intersected ones
    int nribscut = 0;
    TPZManVector<int64_t,2> CutRibsIndex(2);
    TPZVec<int64_t> rib_index = face.GetRibs();
    int nribs = rib_index.size();
    for (int irib = 0; irib < nribs; irib++)
    {
        TRSRibs *ribtest = &fRibs[rib_index[irib]];
        if(ribtest->IsCut() == true){
            CutRibsIndex[nribscut]=rib_index[irib];
            nribscut++;
        }
    }
    // std::cout<<"first rib: "<<CutRibsIndex[0]<<std::endl;
    // std::cout<<"second rib: "<<CutRibsIndex[1]<<std::endl;

    // Connect intersection points
    TPZVec<int64_t> ipoints(2);
    int64_t gelpointIndex = Rib(CutRibsIndex[0])->IntersectionIndex();
    ipoints[0] = fGMesh->Element(gelpointIndex)->NodeIndex(0);
    gelpointIndex = Rib(CutRibsIndex[1])->IntersectionIndex();
    ipoints[1] = fGMesh->Element(gelpointIndex)->NodeIndex(0);
    int64_t nels = fGMesh->NElements();
    fGMesh->CreateGeoElement(EOned, ipoints, 46, nels);
    
    
    int index= face.ElementIndex();
    fMidFaces[index]=face;
}

/**
 * @brief Add faces that are cut at the edges of fracture (using indexes)
 * @param Face to be set
 */

void TRSFractureMesh::AddEndFace(TRSFace &face){
    // Create geometric element for intersection node beetween EndFace and fracture edge
        TPZVec<REAL> coords = FindEndFracturePoint(face);
        TPZVec<int64_t> nodeindex(1,0);
        nodeindex[0] = fGMesh->NodeVec().AllocateNewElement();
        fGMesh->NodeVec()[nodeindex[0]].Initialize(coords, *fGMesh);
        int64_t nels = fGMesh->NElements();
        fGMesh->CreateGeoElement(EPoint, nodeindex, 45, nels); // can't get 1D neighbours from nodes... so had to create GeoEl EPoint for all intersection points and track GeoEl index instead of nodeindex
        face.SetIntersectionIndex(nels);
        // std::cout<<"\n"<<face.ElementIndex
    //iterate over ribs to connect intersection points
    TPZVec<int64_t> rib_index = face.GetRibs();
    int nribs = rib_index.size();
    for (int irib = 0; irib < nribs; irib++)
    {
        TRSRibs *ribtest = &fRibs[rib_index[irib]];
        if(ribtest->IsCut() == true){
            // std::cout<<"single rib cut: "<<rib_index[irib]<<std::endl;
            // Connect intersection points
                nels++;
                TPZVec<int64_t> ipoints(2);
                int64_t gelpointIndex = Rib(rib_index[irib])->IntersectionIndex();
                ipoints[0] = fGMesh->Element(gelpointIndex)->NodeIndex(0);
                ipoints[1] = nodeindex[0];
                fGMesh->CreateGeoElement(EOned, ipoints, 46, nels);
        }
    }
    int index= face.ElementIndex();
    fEndFaces[index]=face;
}









TRSFace * TRSFractureMesh::Face(int64_t index){
    if(fMidFaces.count(index)){
        return &fMidFaces.at(index);
    }
    return &fEndFaces.find(index)->second;
}




/**
 * @brief Create cut surfaces
 * @param Material id
 * @comment This method was getting too long. Moved some of it into AddEndFace and AddMidFace.
 */

void TRSFractureMesh::SplitFaces(int matID){
    //CreateSkeletonElements(1, fSurfaceMaterial);
    // First, build a skeleton for fracplane
    fGMesh->BuildConnectivity();
    // TPZGeoEl *gel = fGMesh->Element(fFracplaneindex);
    // int nsides = gel->NSides();
    // for (int iside = 0; iside < nsides; iside++){
    //     TPZGeoElSide gelside = gel->Neighbour(iside);
    //     if (gelside.Dimension() != 1){continue;}
    //     bool haskel = HasEqualDimensionNeighbour(gelside);
    //     if (haskel == false){TPZGeoElBC(gelside, fSurfaceMaterial);}
    // }

    // iterate over all 2D elements and check their 1D neighbours for intersections
    int64_t nel = fGMesh->NElements();
    for(int iel = 0; iel<nel; iel++){
        TPZGeoEl *gel = fGMesh->Element(iel);
        if (gel->Dimension() != 2){continue;}
        //fSurfaceMaterial is MaterialID for fracture plane
        if(gel->MaterialId()==fSurfaceMaterial){continue;}
        
        int nribscut =0;
        int nedges = gel->NNodes();

        // vector with status for each node and rib of face
        TPZVec<bool> sidestatus(nedges*2,false);
        TPZManVector<int64_t,2> CutRibsIndex(2);
        // vector with indices of father ribs that outline the face
        TPZVec<int64_t> rib_index(nedges,-1);

        // iterate over ribs to check for intersection
        for(int iside = 0; iside < nedges; iside++){
            TPZGeoElSide gelside(gel,iside+nedges);
            TPZGeoElSide neig = gelside.Neighbour();
            while(neig.Element()->Dimension()!=1){
                neig=neig.Neighbour();
            }
            rib_index[iside] = neig.Element()->Index();
            TRSRibs *ribtest = &fRibs[rib_index[iside]];
            if(ribtest->IsCut()==true){
            //check if ribtest was divided into two ribs, or a rib and a point
                // get node where ribtest is cut
                int64_t cutnode = fGMesh->Element(ribtest->IntersectionIndex())->NodeIndex(0);
                TPZVec<int64_t> ribtestNodes(2);
                fGMesh->Element(ribtest->ElementIndex())->GetNodeIndices(ribtestNodes);
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
                
                // store cut rib (this might be discarded once I move "fracplane outline" out of this method)
                CutRibsIndex[nribscut]=rib_index[iside];
                
                nribscut++;
            }
        }

        // if there are ribs cut, create a face object
        if(nribscut == 0){continue;}
        // Create TRSFace
        TRSFace face(iel, true);
        face.SetRibs(rib_index); 
        gel->SetMaterialId(matID);
        face.SetStatus(sidestatus);
        face.SetFractureMesh(this);
        // if(nribscut == 1) {gel->SetMaterialId(matID+17);} // this is here for graphical debugging only... comment it on release

        // Add face to map
        switch (nribscut){
            case 2:AddMidFace(face);break;
            case 1:AddEndFace(face);break;
            default: std::cout<<"\nNo more than 2 ribs should've been cut\n";DebugStop();
        }
        face.DivideSurface(20);
        // // Give children a father and add them to map
        // for(int64_t j; j<children.size(); j++){
        //     TRSFace jchild(children[j], false);
        //     jchild.SetFather(gel);
        //     gel->SetFather ?
        // }
    }
}














TPZVec<REAL> TRSFractureMesh::FindEndFracturePoint(TRSFace &face){
    // Convert TPZGeoEl into TRSFracPlane
    TPZGeoEl *gelface = fGMesh->Element(face.ElementIndex());
    TPZFMatrix<REAL> corners;
    gelface->NodesCoordinates(corners);
    TRSFracPlane faceplane(corners);

    // Check fFracplane's ribs for intersection with faceplane
    int nribs = fFracplane.GetCornersX().Cols();
    for(int irib = 0; irib < nribs; irib++){
        TPZVec<REAL> p1(3);
        TPZVec<REAL> p2(3);
        for(int i = 0; i<3; i++){
            p1[i] = fFracplane.GetCornersX()(i, irib);
            p2[i] = fFracplane.GetCornersX()(i, (irib+1)%nribs);
        }
        if(faceplane.Check_rib(p1, p2)){
            return faceplane.CalculateIntersection(p1,p2);
        }
    }
    std::cout<<"\n TRSFractureMesh::FindEndFracturePoint\n";
    std::cout << "\nFailed to find intersection point in end-fracture face index: " << face.ElementIndex() << std::endl;
    DebugStop();
    return -1;
}







void TRSFractureMesh::SplitRibs(int matID){
    //search gmesh for cut ribs
    int64_t Nels = fGMesh->NElements();
    for (int iel = 0; iel < Nels; iel++){
        TPZGeoEl *gel = fGMesh->Element(iel);
        int nSides = gel->NSides();

        //skip all elements that aren't ribs
        if (gel->Dimension() != 1){continue;}

        // Get rib's vertices
        int64_t p1 = gel->SideNodeIndex(2, 0);
        int64_t p2 = gel->SideNodeIndex(2, 1);
        TPZVec<REAL> pp1(3);
        TPZVec<REAL> pp2(3);
        fGMesh->NodeVec()[p1].GetCoordinates(pp1);
        fGMesh->NodeVec()[p2].GetCoordinates(pp2);

        // Check rib
        bool resul = fFracplane.Check_rib(pp1, pp2);

        // Split rib
        if (resul == true){
            TRSRibs rib(iel, true);
            AddRib(rib);
            TPZVec<REAL> ipoint = fFracplane.CalculateIntersection(pp1, pp2);
            //5O is the material of children ribs
            Rib(iel)->DivideRib(fGMesh, ipoint, matID);
            gel->SetMaterialId(fTransitionMaterial);

            // std::cout<<"Element: "<<iel<<" Side: "<<side<<" Rib status: "<<resul<<" Fracture : 1"<<"\n";
        }
    }
}








void TRSFractureMesh::SplitFractureEdge(){  

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
        TRSFace *iface = &it->second;
        int64_t ipointindex = fGMesh->Element(iface->IntersectionIndex())->NodeIndex(0);
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
            // check if point is in edge by calculating if it's normal distance to the is zero
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
    
    //Once intersections on fracture-edges have been properly ordered and mapped by edge
	//iterate over edges to split them
	for (int iedge = 0; iedge < nedges; iedge++)
	{
        if(edgemap[iedge]->size() == 0){continue;}

		int64_t nels = fGMesh->NElements();
		TPZVec<int64_t> inodes(2);     //index of nodes to be connected
        //@FixIt
		// TPZGeoEl *geofrac = fGMesh->Element(fFracplaneindex); //pointer to fracture GeoEl

		// connect first end-face intersection to iedge's first node
		auto it = edgemap[iedge]->begin();
        inodes[0] = fGMesh->Element(fFracplane.PointElIndex(iedge))->NodeIndex(0);
		inodes[1] = it->second;
		fGMesh->CreateGeoElement(EOned, inodes, 46, nels);

        
        // iterate over iedge's map
        while(++it != edgemap[iedge]->end()){
            nels++;
            inodes[0] = inodes[1];
            inodes[1] = it->second;
            fGMesh->CreateGeoElement(EOned, inodes, 46, nels);
        }

		// connect last end-intersection to edge last node
		nels++;
		inodes[0] = inodes[1];
        inodes[1] = fGMesh->Element(fFracplane.PointElIndex((iedge+1)%nedges))->NodeIndex(0);
		fGMesh->CreateGeoElement(EOned, inodes, 46, nels);
	}
	
// fGMesh->BuildConnectivity();  ?
}












void TRSFractureMesh::SplitFracturePlane(){
    //first, list all ipoints (points of intersection of mesh and fracplane)
    std::unordered_set<int64_t> ipoints;
	
    // iterate over ribs and get their ipoints
    for(auto itr = fRibs.begin(); itr != fRibs.end(); itr++){
        TRSRibs *irib = &itr->second;
		if (irib->IsCut() == false){continue;}
        int64_t index = irib->IntersectionIndex();
        ipoints.insert(index);
    }
    // iterate over endFaces and get their ipoints
    for(auto itr = fEndFaces.begin(); itr != fEndFaces.end(); itr++){
        TRSFace *iface = &itr->second;
		if (iface->IsCut() == false){continue;}
        int64_t index = iface->IntersectionIndex();
        ipoints.insert(index);
    }
    // points of fracture plane (corners) require repetition of somewhat expensive computations and triangulation will always terminate without using them as center point anyway because mesh is composed of convex polyhedra 
	// // List points from fracture plane
    // int ncorners = fFracplane.GetCornersX().Cols();
    // for(int i = 0; i<ncorners; i++){
    //     ipoints.insert(fFracplane.PointElIndex(i));
    // }
	

    //iterate over ipoints
	fGMesh->BuildConnectivity();
    for(auto itr = ipoints.begin(); itr != ipoints.end(); itr++){
		// map of opposite nodes sorted by (cosine of) angle measured from a reference segment
        std::map<REAL, int64_t> oppositenodes;
		TPZGeoEl *ipoint = fGMesh->Element(*itr);
		int64_t centerindex = ipoint->NodeIndex(0);
        TPZGeoElSide neig0 = ipoint->Neighbour(0);
        // Materials from 40 through 50 are being used for fracture mesh during debbuging
        // while (neig0.Element()->Dimension() != 1 || neig0.Element()->MaterialId() != fSurfaceMaterial){
        while (neig0.Element()->Dimension() != 1 || neig0.Element()->MaterialId() < fSurfaceMaterial){
            neig0 = neig0.Neighbour();
        }

        // Reference segment vector
        TPZManVector<REAL,3> ipointcoord(3);
        // Sort edges using angle to the reference vector (scope to keep some variables local)
        {
            // Cosine of angle of segment from inode to opposite node (oppnode) relative to reference segment
            REAL gamma = 1.0; 
            // Opposite node index
            int64_t oppnode;
            if(neig0.Element()->NodeIndex(0) == centerindex){
                oppnode = neig0.Element()->NodeIndex(1);
            }
            else{ oppnode = neig0.Element()->NodeIndex(0);}
            // Insert reference opposite node in map
            int64_t reference = oppnode;
            oppositenodes.insert({gamma,reference});
            // Reference segment vector
            // TPZManVector<REAL,3> ipointcoord(3);
            fGMesh->NodeVec()[centerindex].GetCoordinates(ipointcoord);
            TPZManVector<REAL,3> vref(3);
            fGMesh->NodeVec()[oppnode].GetCoordinates(vref);
            vref[0] -= ipointcoord[0];
            vref[1] -= ipointcoord[1];
            vref[2] -= ipointcoord[2];
            // normalize reference segment
            REAL norm = sqrtl(vref[0]*vref[0]+vref[1]*vref[1]+vref[2]*vref[2]);
            vref[0] = vref[0]/norm;
            vref[1] = vref[1]/norm;
            vref[2] = vref[2]/norm;

            // Iterate over next 1D neighbours to sort segments from it according to angle to reference
            TPZGeoElSide neig_i = neig0.Neighbour();
            while(neig_i != neig0){
                // if(neig_i.Element()->Dimension() == 1 && neig_i.Element()->MaterialId() == fSurfaceMaterial){
                // materials over 40 are being used during development to ifentify elements at fracture surface
                if(neig_i.Element()->Dimension() == 1 && neig_i.Element()->MaterialId() >= fSurfaceMaterial){
                    // get opposite node 
                    if(neig_i.Element()->NodeIndex(0) == centerindex){
                        oppnode = neig_i.Element()->NodeIndex(1);
                    }
                    else{ oppnode = neig_i.Element()->NodeIndex(0);}
                    // Compute segment vector
                    TPZManVector<REAL,3> vi(3);
                    fGMesh->NodeVec()[oppnode].GetCoordinates(vi);
                    vi[0] -= ipointcoord[0];
                    vi[1] -= ipointcoord[1];
                    vi[2] -= ipointcoord[2];
                    // normalize segment vector
                    norm = sqrtl(vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
                    vi[0] = vi[0]/norm;
                    vi[1] = vi[1]/norm;
                    vi[2] = vi[2]/norm;
                    // compute cosine of angle
                    gamma = vref[0]*vi[0]
                            +vref[1]*vi[1]
                            +vref[2]*vi[2];
                    // cross product to determine quadrant of angle
                    TPZManVector<REAL,3> cross(3);
                    cross[0] = vref[1]*vi[2] - vref[2]*vi[1];
                    cross[1] = vref[2]*vi[0] - vref[0]*vi[2];
                    cross[2] = vref[0]*vi[1] - vref[1]*vi[0];
                    // Dot product between cross and fracplane's normal vector
                    REAL dot=cross[0]*fFracplane.axis(0,2)
                            +cross[1]*fFracplane.axis(1,2)
                            +cross[2]*fFracplane.axis(2,2);
                    if(dot < 0){
                        gamma = -(2+gamma);
                    }
                    // Insert opposite node in the map
                    oppositenodes.insert({gamma, oppnode});
                }
                neig_i = neig_i.Neighbour();
            }
        }

        // iterate over sorted opposite nodes 
        for(auto iedge = oppositenodes.rbegin(); iedge != oppositenodes.rend(); iedge++){
			auto nextedge = iedge;
			nextedge++;
			// nextedge = (iedge+1)%nedges
			if(nextedge == oppositenodes.rend()){nextedge = oppositenodes.rbegin();}
            
			// check if element to be created already exists
			bool test = false;
			TPZGeoElSide neighbour = neig0.Neighbour();
			while(neighbour != neig0){
				if(neighbour.Element()->Dimension() == 2 && neighbour.Element()->MaterialId() >=fSurfaceMaterial){
					TPZManVector<int64_t,3> testnodes(3);
					neighbour.Element()->GetNodeIndices(testnodes);
					for(int jnode = 0; jnode < 3; jnode++){
						if(testnodes[jnode] == centerindex ||
							testnodes[jnode] == iedge->second ||
							testnodes[jnode] == nextedge->second)
						{	
							test = true;
							continue;
						}
						else{
							test = false; 
							break;
						}
					}
				}
				if(test){break;}
				neighbour = neighbour.Neighbour();
			}
			if(test){continue;}
			
			// if angle between consectutive segments is 180, then it's the edge of the mesh and no face should be created
			TPZManVector<REAL,3> v1(3), v2(3);
			fGMesh->NodeVec()[iedge->second].GetCoordinates(v1);
				v1[0] -= ipointcoord[0];
				v1[1] -= ipointcoord[1];
				v1[2] -= ipointcoord[2];
				REAL norm = sqrtl(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
				v1[0] = v1[0]/norm;
				v1[1] = v1[1]/norm;
				v1[2] = v1[2]/norm;
			fGMesh->NodeVec()[nextedge->second].GetCoordinates(v2);
				v2[0] -= ipointcoord[0];
				v2[1] -= ipointcoord[1];
				v2[2] -= ipointcoord[2];
				norm = sqrtl(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
				v2[0] = v2[0]/norm;
				v2[1] = v2[1]/norm;
				v2[2] = v2[2]/norm;
				// compute cosine of angle
				REAL gamma = v1[0]*v2[0]
					    	+v1[1]*v2[1]
					    	+v1[2]*v2[2];
			// cos(180) = -1
			if(fabs(gamma+1.)<fTolerance){continue;}
			// passed verifications, then
            // create GeoEl triangle
			TPZManVector<int64_t,3> nodeindices(3);
			nodeindices[0] = centerindex;
			nodeindices[1] = iedge->second;
			nodeindices[2] = nextedge->second;
			int64_t nelements = fGMesh->NElements();
			fGMesh->CreateGeoElement(ETriangle, nodeindices, 47, nelements);
            
			// map it as part of the fracture surface maybe?
                // int nedges = map.size
                // { ipoint, opposite-node[i], opposite-node[(i+1)%nedges] }
			
			
			// create 1D GeoEl for opposite nodes
			// check if segment to be created already exists first
			fGMesh->BuildConnectivity();
			TPZGeoEl *newface = fGMesh->Element(nelements);
			for (int iside = 3; iside < 6; iside++){
				TPZGeoElSide gelside = newface->Neighbour(iside);
				if (gelside.Dimension() != 1){continue;}
				bool haskel = HasEqualDimensionNeighbour(gelside);
				if (haskel == false)
				{
					TPZGeoElBC(gelside, 48);
				}
			}
            
		}
	oppositenodes.clear();
    }
}


void TRSFractureMesh::AddVolume(TRSVolume volume){
    int index = volume.ElementIndex();
    fVolumes[index] = volume;
}



/// Sets material for elements at surface of fracture (40 is default)
void TRSFractureMesh::SetSurfaceMaterial(int matID = 40){
    int64_t nels = fGMesh->NElements();
    for (int64_t iel = 0; iel < nels; iel++)
    {
        if(fGMesh->Element(iel)->MaterialId() == fSurfaceMaterial){
            fGMesh->Element(iel)->SetMaterialId(matID);
        }
    }
    fSurfaceMaterial = matID;
}




void TRSFractureMesh::CreateVolumes(){
    // map all volumes that have neighbour faces cut
    int64_t nels = fGMesh->NElements();
    for (int64_t iel = 0; iel < nels; iel++){
        TPZGeoEl *gel = fGMesh->Element(iel);
        if(gel->Dimension() != 3){continue;}
        int nsides = gel->NSides();
        // 2D sides wont start before index 9 for any 3D element
        for (int iside = 9; iside < nsides; iside++){
            TPZGeoElSide gelside = gel->Neighbour(iside);
            if (gelside.Dimension() != 2){continue;}
            if(Face(gelside.Element()->Index())->IsCut()){
                TRSVolume volume(iel,true);
                AddVolume(volume);
                break;
            }            
        }
    }

    // search through each element of the triangulated fracture surface to find their enclosing volume
    // iterate over fracplane's elements created at SplitFracturePlane
    for(int64_t iel = 0; iel < nels; iel++){
        TPZGeoEl *gel = fGMesh->Element(iel);
        // if(gel->MaterialID() != fSurfaceMaterial){continue;}
        // During development, elements at fracture surface have material id over 40
        if(gel->MaterialId() <= fSurfaceMaterial){continue;}
        if(gel->Dimension() != 2){continue;}

        // Find volume that encloses that element
        FindEnclosingVolume(gel);
    }

    // Use GMSH to tetrahedralize volumes
    std::ofstream outfile("fracture1.geo");
    WriteGMSH(outfile);

}









bool TRSFractureMesh::FindEnclosingVolume(TPZGeoEl *ifracface){
    // get coordinates of geometric center of face
    TPZVec<REAL> faceCenter(3);
    {
        TPZGeoElSide geliside(ifracface, ifracface->NSides()-1);
        geliside.CenterX(faceCenter);
    }

    // map of indices for volumes that could contain the face
    std::map<REAL, int64_t> candidates;
    
    // iterate over ifracface 1D sides 
    int nsides = ifracface->NSides();
    for(int iside = 0; iside < nsides; iside++){
        if(ifracface->SideDimension(iside) != 1){continue;}
        TPZGeoElSide geliside(ifracface, iside);

        // iterate over neighbours through iside
        TPZGeoElSide ineig = geliside.Neighbour();
        for( ; ineig != geliside; ineig = ineig.Neighbour()){
            // ignore elements at fracture surface
            // if(ineig.Element()->MaterialId() == fSurfaceMaterial){continue;}
            // During development, elements at fracture surface have material id over 40
            if(ineig.Element()->MaterialId() >= fSurfaceMaterial){continue;}
            
            // find 2-dimensional neighbour that has a father
            if(ineig.Element()->Dimension() != 2){continue;}
            TPZGeoEl *father = ineig.Element()->Father();
            if(!father){continue;}

            // get father's center coordinates
            TPZManVector<REAL,3> fatherCenter(3);
            TPZGeoElSide fatherfaceside(father, father->NSides()-1);
            fatherfaceside.CenterX(fatherCenter);
            // construct vector from center of father to center of ifracface
            TPZManVector<REAL,3> v1(3,0);
                v1[0] = faceCenter[0] - fatherCenter[0];
                v1[1] = faceCenter[1] - fatherCenter[1];
                v1[2] = faceCenter[2] - fatherCenter[2];
            // Normalize v1
            REAL norm = sqrtl(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
                v1[0] = v1[0]/norm;
                v1[1] = v1[1]/norm;
                v1[2] = v1[2]/norm;

            // iterate over volumetric neighbours through father's face
            TPZGeoElSide ivolume = fatherfaceside.Neighbour();
            for( ; ivolume != fatherfaceside; ivolume = ivolume.Neighbour()){
                if(ivolume.Element()->Dimension() != 3){continue;}
                // get coordinates for center of volume
                TPZManVector<REAL,3> volumeCenter(3);
                {
                    TPZGeoElSide gelsidevolume (ivolume.Element(),ivolume.Element()->NSides()-1);
                    gelsidevolume.CenterX(volumeCenter);
                }

                // construct vector from center of ifracface to center of volume
                TPZManVector<REAL,3> v2(3,0);
                    v2[0] = volumeCenter[0] - fatherCenter[0];
                    v2[1] = volumeCenter[1] - fatherCenter[1];
                    v2[2] = volumeCenter[2] - fatherCenter[2];
                // Normalize v2
                norm = sqrtl(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
                    v2[0] = v2[0]/norm;
                    v2[1] = v2[1]/norm;
                    v2[2] = v2[2]/norm;
                
                // if dot product between the vectors constructed for centers is
                // positive, that volume is a candidate
                REAL dot = 0;
                for(int ico = 0; ico < 3; ico++){dot += v1[ico]*v2[ico];}
                if(dot>0){
                    candidates[dot] = ivolume.Element()->Index();
                }
            }
        }
    }
    
    // return best candidate 
    if(candidates.size() > 0){
        // reverse iterator (rbegin) gives biggest key in map
        int64_t volumeindex = candidates.rbegin()->second;
        fVolumes[volumeindex].SetFaceInVolume(ifracface->Index());
        return true;
    }

    // degeneracy: ifracface's edges are completely enclosed by volume
    std::set<TPZGeoElSide *> verified;
    int nnodes = ifracface->NNodes();
    for(int iside = 0; iside<nnodes; iside++){
        if(ifracface->SideDimension(iside) != 0){continue;}
        TPZGeoElSide gelsidenode(ifracface, iside);

        // iterate over neighbours through gelsidenode
        TPZGeoElSide ineig = gelsidenode.Neighbour();
        for( ; ineig != gelsidenode; ineig = ineig.Neighbour()){
            // ignore elements at fracture surface
            // if(ineig.Element()->MaterialId() == fSurfaceMaterial){continue;}
            // During development, elements at fracture surface have material id over 40
            if(ineig.Element()->MaterialId() >= fSurfaceMaterial){continue;}

            // find 2-dimensional neighbour that has a father
            if(ineig.Element()->Dimension() != 2){continue;}
            TPZGeoEl *father = ineig.Element()->Father();
            if(!father){continue;}
            TPZGeoElSide fatherfaceside = father->Neighbour(father->NSides()-1);

           // iterate over volumetric neighbours through father's plane
            TPZGeoElSide ivolume = fatherfaceside.Neighbour();
            for( ; ivolume != fatherfaceside; ivolume = ivolume.Neighbour()){
                if(ivolume.Element()->Dimension() != 3){continue;}
                if(verified.find(&ivolume) != verified.end()){continue;}
                TPZVec<REAL> ksi(3,2.0);
                bool test = ivolume.Element()->ComputeXInverse(faceCenter, ksi, fTolerance);
                if(test == true){
                    int64_t volumeindex = ivolume.Element()->Index();
                    fVolumes[volumeindex].SetFaceInVolume(ifracface->Index());
                    return true;
                }
            }
        }
    }

    std::cout<<"\n TRSFractureMesh::FindEnclosingVolume found no enclosing volume for element #"<<ifracface->Index()<<"\n";
    DebugStop();
    return false;
}










void TRSFractureMesh::WriteGMSH(std::ofstream &outfile){
    
    // Giving fGMesh another name for readability's sake
    TPZGeoMesh *pzgmesh = fGMesh;
    pzgmesh->BuildConnectivity();
    CreateSkeletonElements(1,19);
    // Title
    outfile<<"//  Geo file generated by Discrete Fracture Network methods \n"
            <<"// Fracture #1 \n\n";
    
    // write nodes
    outfile<< "// POINTS DEFINITION \n\n";
    int64_t nnodes = pzgmesh->NNodes();
    for (int64_t inode = 0; inode < nnodes; inode++){
        TPZManVector<REAL, 3> co(3,0.);
        pzgmesh->NodeVec()[inode].GetCoordinates(co);
        outfile << "Point(" << inode << ") = {" << co[0] << ',' << co[1] << ',' << co[2] << "};\n";
    }
    
    // write edges
    int64_t nels = pzgmesh->NElements();
    outfile << "\n\n// LINES DEFINITION \n\n";
    {
        std::list<int64_t> groupSurface;
        std::list<int64_t> groupTransition;
        std::list<int64_t> groupIntact;
        for (int64_t iel = 0; iel < nels; iel++){
            TPZGeoEl *gel = pzgmesh->Element(iel);
            if(gel->Dimension() != 1
                || gel->MaterialId() == fTransitionMaterial) continue;
            outfile << "Line(" << iel << ") = {" << gel->NodeIndex(0) << ',' << gel->NodeIndex(1) << "};\n";
            // this is a mess, but only for debugging
            if(gel->MaterialId() >= 19 && gel->MaterialId() < 40){groupTransition.push_back(iel);}
            else if(gel->MaterialId() >= 40){groupSurface.push_back(iel);}
            else groupIntact.push_back(iel);
        }
        // write physical groups
        outfile<<"\nPhysical Curve("<<19<<") = {";
        for(auto itr = groupTransition.begin(); itr != groupTransition.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupTransition.end()? "," : "};\n");
        }
        outfile<<"\nPhysical Curve("<<40<<") = {";
        for(auto itr = groupSurface.begin(); itr != groupSurface.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupSurface.end()? "," : "};\n");
        }
        outfile<<"\nPhysical Curve("<<1<<") = {";
        for(auto itr = groupIntact.begin(); itr != groupIntact.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupIntact.end()? "," : "};\n");
        }
    }
    // write faces
    outfile << "\n\n// FACES DEFINITION \n\n";
    {
        std::list<int64_t> groupSurface;
        std::list<int64_t> groupTransition;
        std::list<int64_t> groupIntact;
        for (int64_t iel = 0; iel < nels; iel++){
            TPZGeoEl *gel = pzgmesh->Element(iel);
            if(gel->Dimension() != 2
                || gel->MaterialId() == fTransitionMaterial) continue;
            
            int nnodes = gel->NCornerNodes();
            int nedges = nnodes; //for readability 
            TPZManVector<int64_t,4> facenodevec(nnodes);
            gel->GetNodeIndices(facenodevec);
            // line loop
            outfile << "Line Loop(" << iel << ") = {";
            // line loops require a proper orientation of lines
            for(int iside = nedges; iside<2*nedges; iside++){
                TPZGeoElSide gelside(gel,iside);
                TPZGeoElSide side = gelside.Neighbour();
                // find line element
                while(side.Element()->Dimension() != 1) {side = side.Neighbour();}
                // find first node of line
                int inode = 0;
                while(facenodevec[inode] != side.SideNodeIndex(0)) ++inode;
                // check orientation by comparing second node of line with next node of face
                int64_t index=0;
                    if(side.SideNodeIndex(1)==facenodevec[(inode+1)%nedges]){
                        index = side.Element()->Index();
                    }
                    else{
                        index = -side.Element()->Index();
                    }
                outfile << index <<(iside < 2*nedges-1? "," : "};\n");
            }
            // surface
            outfile << "Surface("<<iel<<") = {"<<iel<<"};\n";
            // this is a mess, but only for debugging
            if(gel->MaterialId() >= 19 && gel->MaterialId() < 40){groupTransition.push_back(iel);}
            else if(gel->MaterialId() >= 40){groupSurface.push_back(iel);}
            else groupIntact.push_back(iel);
        }
        // write physical groups
        outfile<<"\nPhysical Surface("<<19<<") = {";
        for(auto itr = groupTransition.begin(); itr != groupTransition.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupTransition.end()? "," : "};\n");
        }
        outfile<<"\nPhysical Surface("<<40<<") = {";
        for(auto itr = groupSurface.begin(); itr != groupSurface.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupSurface.end()? "," : "};\n");
        }
        outfile<<"\nPhysical Surface("<<1<<") = {";
        for(auto itr = groupIntact.begin(); itr != groupIntact.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupIntact.end()? "," : "};\n");
        }
    }

    // write volumes
    outfile << "\n\n// VOLUMES DEFINITION \n\n";
    {

    }
    
    outfile<<"\nTransfinite Surface \"*\";\n";
    outfile<<"Recombine Surface \"*\";\n";
    outfile<<"\nTransfinite Volume \"*\";\n";
}


















