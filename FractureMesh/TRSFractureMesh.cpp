/*! 
 *  @brief     Compares a geomesh with fracture plane to find intersections.
 *  @details   Intersection search is performed after creation of skeleton
 *  elements with TRSFractureMesh::CreateSkeletonElements. Fracture plane should
 *  be a TRSFracPlane.
 *  @authors   Jorge Ordoñez
 *  @authors   Pedro Lima
 *  @date      2018-2019
 */

#include "TRSFractureMesh.h"
#include <math.h>
#include <cstdio>

// Empty Constructor
TRSFractureMesh::TRSFractureMesh(){
}

// Constructor with corner points and a geomesh
TRSFractureMesh::TRSFractureMesh(TRSFracPlane &FracPlane, TPZGeoMesh *gmesh){
    fFracplane = FracPlane;
    fGMesh = gmesh;

    // Implement check for previously created skeleton
    // Create skeleton elements
    CreateSkeletonElements(2, 4);
    CreateSkeletonElements(1, 4);

    // Create FracPlane's geometric element into mesh
    int nnodes =  fGMesh->NNodes();
    int ncorners = fFracplane.GetCorners().Cols();
    fGMesh->NodeVec().Resize(nnodes + ncorners);
    TPZVec<int64_t> CornerIndexes(ncorners);
    for (int i = 0; i < ncorners; i++)
    {
        TPZVec<REAL> nodeX(3, 0);
        nodeX[0] = fFracplane.GetCorners()(0,i);
        nodeX[1] = fFracplane.GetCorners()(1,i);
        nodeX[2] = fFracplane.GetCorners()(2,i);

        fGMesh->NodeVec()[nnodes + i].Initialize(nodeX, *fGMesh);
        CornerIndexes[i] = nnodes + i;
    }
    MElementType elemtype;
    switch (ncorners)
    {
        case 3: elemtype = ETriangle; break;
        case 4: elemtype = EQuadrilateral; break;
        default: DebugStop();
    }
    fFracplaneindex = fGMesh->NElements();
    fGMesh->CreateGeoElement(elemtype, CornerIndexes, 40, fFracplaneindex);
    

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
 * @return The tolerance
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



// /**
//  * @brief Checks if a surface needs to be divided
//  * @param Surface index (integer)
//  * @param Cut ribs vector
//  * @return True if the surface needs to be divided
//  * @return False if the surface does not need to be divided
//  */

// bool TRSFractureMesh::NeedsSurface_Divide(int64_t suface_index, TPZVec<int64_t> interribs) {
//     bool state= false;   //By definition does not need to be divided
//     int nribs = interribs.size();
//     int nribscut =0;
//     for(int i=0; i< nribs; i++){
//         int index_anal = interribs[i];
//         TRSRibs rib_an = fRibs[index_anal];
//         if(rib_an.CutsPlane()==true){
//             nribscut++;
//         }
//     }
//     if(nribscut > 0){     //Checks if a surface has ribs cut
//         return true;   //The surface needs to be divided
//     }
//     return false;
// }

/**
 * @brief Check if the neighbour has a equal dimension
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
                int nel_mesh = fGMesh->NElements();
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
    int index= face.ElementIndex();
    fMidFaces[index]=face;
}

/**
 * @brief Add faces that are cut at the edges of fracture (using indexes)
 * @param Face to be set
 */

void TRSFractureMesh::AddEndFace(TRSFace face){
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
    //CreateSkeletonElements(1, 4);
    //CreateSkeletonElements(2, 4);
    int64_t nel = fGMesh->NElements();
    
    for(int iel = 0; iel<nel; iel++){
        TPZGeoEl *gel = fGMesh->Element(iel);
        if (gel->Dimension() != 2){continue;}
        //40 is MaterialID for fracture plane
        if(gel->MaterialId()==40){continue;}
        
        int nribscut =0;
        int nsides = gel->NNodes();
        TPZVec<bool> ribstatus(nsides,false);
        TPZManVector<int64_t,2> CutRibsIndex(2);
        TPZVec<int64_t> rib_index(nsides,-1);

        for(int iside = 0; iside < nsides; iside++){
            TPZGeoElSide gelside(gel,iside+nsides);
            TPZGeoElSide neig = gelside.Neighbour();
            while(neig.Element()->Dimension()!=1){
                neig=neig.Neighbour();
            }
            rib_index[iside] = neig.Element()->Index();
            TRSRibs *ribtest = &fRibs[rib_index[iside]];
            if(ribtest->IsCut()==true){
                ribstatus[iside] = true;
                // store the rest of the ribs
                CutRibsIndex[nribscut]=rib_index[iside];
                nribscut++;
            }
        }
        
        switch (nribscut)
        {
            case 0: {break;}
            case 2:{ //mid-fracture element
            // Create TRSFace
                TRSFace face(iel, true);
                AddMidFace(face);
				std::cout<<"first rib: "<<CutRibsIndex[0]<<std::endl;
				std::cout<<"second rib: "<<CutRibsIndex[1]<<std::endl;
				face.SetRibs(rib_index); 
				gel->SetMaterialId(matID);
				
			// Connect intersection points
				TPZVec<int64_t> ipoints(2);
				ipoints[0] = Rib(CutRibsIndex[0])->IntersectionIndex();
				ipoints[1] = Rib(CutRibsIndex[1])->IntersectionIndex();
				int64_t nels = fGMesh->NElements();
				fGMesh->CreateGeoElement(EOned, ipoints, 46, nels);

                break;
            }
            case 1:{ //end-fracture element
                TRSFace face(iel, true);
                std::cout<<"single rib cut: "<<CutRibsIndex[0]<<std::endl;
                face.SetRibs(rib_index); 
                gel->SetMaterialId(matID+15);
                //Is the fracture skeleton built? Code won't work otherwise.
                TPZVec<REAL> coords = FindEndFracturePoint(face);
                
                // Create geometric element for intersection node
                TPZVec<int64_t> nodeindex(1,0);
                nodeindex[0] = fGMesh->NodeVec().AllocateNewElement();
                fGMesh->NodeVec()[nodeindex[0]].Initialize(coords, *fGMesh);
                int64_t nels = fGMesh->NElements();
                fGMesh->CreateGeoElement(EPoint, nodeindex, 45, nels);
				face.SetIntersectionIndex(nodeindex[0]);
                AddEndFace(face);

                // Connect intersection points
                nels++;
                TPZVec<int64_t> ipoints(2);
                ipoints[0] = Rib(CutRibsIndex[0])->IntersectionIndex();
                ipoints[1] = nodeindex[0];
                fGMesh->CreateGeoElement(EOned, ipoints, 46, nels);
                break;
            }
            default: {std::cout<<"\nNo more than 2 ribs should've been cut\n";DebugStop();}
        }
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
		//unnecessary for loop... delete it later
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
            // check if point is in edge by calculating if it's normal distance to the is zero
            REAL dist = temp/edgelength[iedge];
			if(dist<fTolerance){
				// compute local 1D coordinate (alpha)
				REAL norm = sqrtl(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
				REAL alpha = norm/edgelength[iedge];
				if(alpha > 1+fTolerance){std::cout<<"\nProblem with alpha\n";DebugStop();}
				// map intersection indexes from smallest alpha to biggest
                // map <alpha, index>
				edgemap[iedge]->insert(std::make_pair(alpha, ipointindex));
			}
        }
    }
    
    //Once intersections on fracture-edges have been properly ordered and mapped by edge
	//iterate over edges to split them
	for (int iedge = 0; iedge < nedges; iedge++)
	{
        if(edgemap[iedge]->size() == 0){continue;}

		int64_t nels = fGMesh->NElements();
		TPZVec<int64_t> ipoints(2);     //index of nodes to be connected
		TPZGeoEl *geofrac = fGMesh->Element(fFracplaneindex); //pointer to fracture GeoEl

		// connect first end-face intersection to iedge's first node
		auto it = edgemap[iedge]->begin();
		ipoints[0] = geofrac->SideNodeIndex(iedge+nedges,0);;
		ipoints[1] = it->second;
		fGMesh->CreateGeoElement(EOned, ipoints, 46, nels);

        
        // iterate over iedge's map
        while(++it != edgemap[iedge]->end()){
            nels++;
            ipoints[0] = ipoints[1];
            ipoints[1] = it->second;
            fGMesh->CreateGeoElement(EOned, ipoints, 46, nels);
        }

		// connect last end-intersection to edge last node
		nels++;
		ipoints[0] = ipoints[1];
		ipoints[1] = geofrac->SideNodeIndex(iedge+nedges,1);
		fGMesh->CreateGeoElement(EOned, ipoints, 46, nels);
	}
	
// fGMesh->BuildConnectivity();  ?
}












void TRSFractureMesh::SplitFracturePlane(){
    //need to map ipoints first

    //iterate over ipoints
        // build connectivity
        // std::map<REAL, int64_t> Sorted1DNeighbours
        // store neig_0
        // get opposite node to neig_0
        // map opposite node of neig_0 with angle = 0 (use it as reference)
        // iterate over next 1D neighbours
            // neighbour neig_i
            // get opposite node 

            // using angle
            // REAL angle
                // angle = arccos(dotProduct(neig_0.normalize(), neig_i.normalize()));
                // cross = crossProduct(neig_0, neig_i);
                // if (dotProduct(fAxis2, cross) < 0) {
                //   angle = 2π - angle;
                // }
            // arccos of dot product to get angle from neig_0 to neig_i
                // (is arccos too expensive? is there another solution?)


            //using cossine
            //REAL cos = dotProduct(neig_0.normalize(), neig_i.normalize());
                // cross = crossProduct(neig_0, neig_i);
                // if (dotProduct(fAxis2, cross) < 0) {
                //   cos = -2 - cos;
                // }


            // insert angle and opposite node in map
        // iterate over Sorted1DNeighbours 
            // check if 2D element already exists {continue;}
            // create GeoEl triangles
            // map it as part of the fracture surface
                // { ipoint, opposite-node[i], opposite-node[(i+1)%map.size] }
            // check if 1D element exists
            // create 1D GeoEl from opposite-node[i] to opposite-node[(i+1)%map.size]

}








//void TRSFractureMesh::CreateTransitionVolumes(){
    // iterate over
// }