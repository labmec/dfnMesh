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

// Constructor with corner points and a geomesh
TRSFractureMesh::TRSFractureMesh(TRSFracPlane &FracPlane, TPZGeoMesh *gmesh){
    fFracplane = FracPlane;
    fGMesh = gmesh;

    // Maybe implement check for previously created skeleton
    // Create skeleton elements
    CreateSkeletonElements(2, 4);
    CreateSkeletonElements(1, 4);

    // Create FracPlane's geometric element into mesh
    fFracplaneindex = fFracplane.CreateElement(fGMesh);
    
    
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
 * @return True if has a lower dimension
 * @return False if has a upper dimension
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
        if(gel->MaterialId() == 40 && dimension == 2){continue;}

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

void TRSFractureMesh::AddMidFace(TRSFace face){
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

void TRSFractureMesh::AddEndFace(TRSFace face){
    // Create geometric element for intersection node
        TPZVec<REAL> coords = FindEndFracturePoint(face);
        TPZVec<int64_t> nodeindex(1,0);
        nodeindex[0] = fGMesh->NodeVec().AllocateNewElement();
        fGMesh->NodeVec()[nodeindex[0]].Initialize(coords, *fGMesh);
        int64_t nels = fGMesh->NElements();
        fGMesh->CreateGeoElement(EPoint, nodeindex, 45, nels);
        face.SetIntersectionIndex(nels);
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

/**
 * @brief Gives the ribs map
 * @return A map that contains the ribs information
 */

std::map<int64_t ,TRSRibs> TRSFractureMesh::GetRibs(){
    return fRibs;    
}

/**
 * @brief Create cut surfaces
 * @param Material id
 * @comment This method is getting too long. Might move some of it into AddEndFace and AddMidFace.
 */

void TRSFractureMesh::SplitFaces(int matID){
    //CreateSkeletonElements(1, 40);
    // First, build a skeleton for fracplane
    fGMesh->BuildConnectivity();
    TPZGeoEl *gel = fGMesh->Element(fFracplaneindex);
    int nsides = gel->NSides();
    for (int iside = 0; iside < nsides; iside++){
        TPZGeoElSide gelside = gel->Neighbour(iside);
        if (gelside.Dimension() != 1){continue;}
        bool haskel = HasEqualDimensionNeighbour(gelside);
        if (haskel == false){TPZGeoElBC(gelside, 40);}
    }

    // iterate over all 2D elements and check their 1D neighbours for intersections
    int64_t nel = fGMesh->NElements();
    for(int iel = 0; iel<nel; iel++){
        TPZGeoEl *gel = fGMesh->Element(iel);
        if (gel->Dimension() != 2){continue;}
        //40 is MaterialID for fracture plane
        if(gel->MaterialId()==40){continue;}
        
        int nribscut =0;
        int nedges = gel->NNodes();

        // vector with status for each rib and node of face
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
                
                // store cut rib (this might be discarded once I fracplane outline out of this method)
                CutRibsIndex[nribscut]=rib_index[iside];
                nribscut++;
            }
        }
// debugging_______________________________________________________
        // if(nribscut == 1){std::cout<<"\nEndFace";}
        // if(nribscut == 2){std::cout<<"\nMidFace";}
        // std::cout<<"\n sidestatus = {";
        // for(int i=0; i<8; i++){std::cout<<sidestatus[i]<<" ";}
        // std::cout<<"}";
// debugging_______________________________________________________

        if(nribscut == 0){continue;}
        // Create TRSFace
        TRSFace face(iel, true);
        face.SetRibs(rib_index); 
        gel->SetMaterialId(matID);
        face.SetStatus(sidestatus);
        if(nribscut == 1) {gel->SetMaterialId(matID+15);} // this is here for graphical debugging only... comment it on release

        // Add face to map
        switch (nribscut){
            case 2:AddMidFace(face);break;
            case 1:AddEndFace(face);break;
            default: std::cout<<"\nNo more than 2 ribs should've been cut\n";DebugStop();
        }
        face.DivideSurface(this, 2);
    }
}






TRSRibs *TRSFractureMesh::Rib(int index){
    return &fRibs[index];
}







TPZVec<REAL> TRSFractureMesh::FindEndFracturePoint(TRSFace &face){
    // Convert TPZGeoEl into TRSFracPlane
    TPZGeoEl *gelface = fGMesh->Element(face.ElementIndex());
    TPZFMatrix<REAL> corners;
    gelface->NodesCoordinates(corners);
    TRSFracPlane faceplane(corners);

    // Check fFracplane's ribs for intersection with faceplane
    int nribs = fFracplane.GetCorners().Cols();
    for(int irib = 0; irib < nribs; irib++){
        TPZVec<REAL> p1(3);
        TPZVec<REAL> p2(3);
        for(int i = 0; i<3; i++){
            p1[i] = fFracplane.GetCorners()(i, irib);
            p2[i] = fFracplane.GetCorners()(i, (irib+1)%nribs);
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
		//unnecessary "for" loop... delete it later
        for (int side = 0; side < nSides; side++){
			// this isn't doing anything... delete it later
            if (gel->NSideNodes(side) == 2){

                // Get rib's vertexes
                int64_t p1 = gel->SideNodeIndex(side, 0);
                int64_t p2 = gel->SideNodeIndex(side, 1);
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
                    gel->SetMaterialId(12);

                    // std::cout<<"Element: "<<iel<<" Side: "<<side<<" Rib status: "<<resul<<" Fracture : 1"<<"\n";
                }
            }
        }
    }
}








void TRSFractureMesh::SplitFractureEdge(){  

    // Compute fracture edges' lenghts
    int nedges = fFracplane.GetCorners().Cols();
    TPZVec<REAL> edgelength(nedges,0);
    Matrix fraccorners(3,nedges);
    fraccorners = fFracplane.GetCorners();
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
		TPZGeoEl *geofrac = fGMesh->Element(fFracplaneindex); //pointer to fracture GeoEl

		// connect first end-face intersection to iedge's first node
		auto it = edgemap[iedge]->begin();
		inodes[0] = geofrac->SideNodeIndex(iedge+nedges,0);;
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
		inodes[1] = geofrac->SideNodeIndex(iedge+nedges,1);
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
    // can't get 1D neighbours from nodes... so had to create GeoEl EPoint for all ipoints and track GeoEl index instead of nodeindex
	
	
	// GeoPoints for fracture vertices aren't needed if fracture plane is all contained by gmesh
    // (since all polygons generated are convex)
		// int ncorners = fGMesh->Element(fFracplaneindex)->NCornerNodes();
		// for(int i = 0; i<ncorners; i++){
		// 	TPZVec<int64_t> cornerindex(1);
		// 	cornerindex[0] = fGMesh->Element(fFracplaneindex)->NodeIndex(i);
		// 	int64_t nels = fGMesh->NElements();
		// 	fGMesh->CreateGeoElement(EPoint, cornerindex, 45, nels);
		// 	ipoints.insert(nels);
		// }

    //iterate over ipoints
	fGMesh->BuildConnectivity();
    for(auto itr = ipoints.begin(); itr != ipoints.end(); itr++){
		// map of opposite nodes sorted by (cosine of) angle measured from a reference segment
        std::map<REAL, int64_t> oppositenodes;
		TPZGeoEl *ipoint = fGMesh->Element(*itr);
		int64_t centerindex = ipoint->NodeIndex(0);
        TPZGeoElSide neig0 = ipoint->Neighbour(0);
		// Materials from 40 to 50 are being used for fracture mesh during debbuging
        while (neig0.Element()->Dimension() != 1 || neig0.Element()->MaterialId() < 40){
            neig0 = neig0.Neighbour();
        }

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
        TPZVec<REAL> ipointcoord(3);
        fGMesh->NodeVec()[centerindex].GetCoordinates(ipointcoord);
        TPZVec<REAL> vref(3);
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
			if(neig_i.Element()->Dimension() == 1 && neig_i.Element()->MaterialId() >= 40){
				// get opposite node 
				if(neig_i.Element()->NodeIndex(0) == centerindex){
					oppnode = neig_i.Element()->NodeIndex(1);
				}
				else{ oppnode = neig_i.Element()->NodeIndex(0);}
				// Compute segment vector
				TPZVec<REAL> vi(3);
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
				TPZVec<REAL> cross(3);
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

        // iterate over sorted opposite nodes 
        for(auto iedge = oppositenodes.begin(); iedge != oppositenodes.end(); iedge++){
			auto nextedge = iedge;
			nextedge++;
			// nextedge = (iedge+1)%nedges
			if(nextedge == oppositenodes.end()){nextedge = oppositenodes.begin();}
            
			// check if element to be created already exists
			bool test = false;
			TPZGeoElSide neighbour = neig0.Neighbour();
			while(neighbour != neig0){
				if(neighbour.Element()->Dimension() == 2 && neighbour.Element()->MaterialId() >=40){
					TPZVec<int64_t> testnodes(3);
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
			TPZVec<REAL> v1(3), v2(3);
			fGMesh->NodeVec()[iedge->second].GetCoordinates(v1);
				v1[0] -= ipointcoord[0];
				v1[1] -= ipointcoord[1];
				v1[2] -= ipointcoord[2];
				norm = sqrtl(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
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
				gamma =  v1[0]*v2[0]
						+v1[1]*v2[1]
						+v1[2]*v2[2];
			// cos(180) = -1
			if(fabs(gamma+1.)<fTolerance){continue;}

			// passed verifications, then
            // create GeoEl triangle
			TPZVec<int64_t> nodeindices(3);
			nodeindices[0] = centerindex;
			nodeindices[1] = iedge->second;
			nodeindices[2] = nextedge->second;
			int64_t nelements = fGMesh->NElements();
			fGMesh->CreateGeoElement(ETriangle, nodeindices, 47, nelements);
            
			// map it as part of the fracture surface maybe?
                // int nedges = map.size
                // { ipoint, opposite-node[i], opposite-node[(i+1)%nedges] }
			
			
			// create 1D GeoEl for opposite nodes
			// check if segment to be created already exists
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





















//void TRSFractureMesh::CreateTransitionVolumes(){

/*
    iterate over fMidFaces
        map 3D neighbours
    iterate over fEndFaces
        map 3D neighbours
 
    this map is called fTranVolumes (Transition volumes)
    maybe create a class for TranVolumes
    implement its splitting method in there
 
    build connectivity
    iterate over fTranVolumes (ivolume)
        iterate over ivolume 2D neighbours (iplane)
            if iplane has subfaces {continue;}
            create set
            iterate over iplane's 2D neighbours (first_neighbour)
                if set.size == 3 {breake;}
                if neighbourhood happens at a node {continue;}
                if first_neighbour has subface {continue}
                if first_neighbour is in fracplane {continue;}
                iterate over iplane's 2D neighbours from (first+1)_neighbour (second_neighbour)
                    if neighbourhood happens at a node {continue;}
                    if second_neighbour is also neighbour of first_neighbour 
                        if the neighbourhood first-to-second happens on a different edge than that of first-to-iplane
                            insert iplane in set
                            insert first_neighbour in set
                            insert second_neighbour in set
                                // only 3 faces are needed to get a tetrahedron's corner nodes to create GeoElement
                            breake
            if (iplane.nnodes == 4 || first_neighbour.nnodes == 4 || second_neighbour.nnodes == 4)
                if set.size == 4 {breake;}
                iterate over iplane's 2D neighbours (third_neighbour)
                    if neighbourhood happens at a node {continue;}
                    if third_neighbour is already in set {continue;}
                    if third_neighbour is neighbour of either first or second neighbours
                        if the neighbourhood third_to_(first or second) happens on a different edge than that of (first or second)_to_iplane
                            insert third_neighbour in set
                            breake;
        
            // get corner nodes
            TPZVec<int64_t> nodeindices(0)
            iterate over set
                index = get corner nodes for each face in set
                nodeindices.push_back(index)
            remove duplicates in nodeindices

            // determine element type
            switch nodeindices.size
                case 4: tetraedro
                case 6: pyramid
                case 8: hexaedro
            
            // check if volume to be created already exists
            bool test = false;
            iterate over iplane's 3D neighbours (volume_neig)
                if volume_neig.elementType != newvolume.elementType {continue;}
                get volume_neig nodes indices
                compare all volume_neig nodes with each nodeindices[i]
                    test = true;
                    breake loop if they match
            if(test == true){continue;}


            // create volume element
            int64_t nelements = fGMesh->NElements();
            fGMesh->CreateGeoElement(ElementType, nodeindices, 49, nelements);
            

            // then check if new volume has skeleton
            fGMesh->BuildConnectivity();
            TPZGeoEl *newvolume = fGMesh->Element(nelements);
            int nsides = newvolume->NSides();
            for (int iside = 0; iside < nsides; iside++){
                TPZGeoElSide gelside = newvolume->Neighbour(iside);
                if (gelside.Dimension() != 2){continue;}
                bool haskel = HasEqualDimensionNeighbour(gelside);
                if (haskel == false)
                {
                    TPZGeoElBC(gelside, 49);
                }
            }
 */

// }