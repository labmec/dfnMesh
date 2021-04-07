/*! 
 *	DFNFace.cpp
 *  @authors   Pedro Lima
 *  @authors   Philippe Devloo
 *  @date      2018-2020
 */

#include "DFNFace.h"
#include "DFNFracture.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

#ifdef LOG4CXX
    static LoggerPtr logger(Logger::getLogger("dfn.mesh"));
#endif

//Constructor
DFNFace::DFNFace(TPZGeoEl *gel, DFNFracture *fracture, TPZVec<DFNRib *> &ribvec) :
	fGeoEl(gel),
	fFracture(fracture),
	fCoord(0)
{
	int nsides = gel->NSides();
	fStatus.Resize(nsides,0);
//	int nedges = gel->NCornerNodes();
//	fRibs.Resize(nedges,nullptr);
    SetRibs(ribvec);
	UpdateStatusVec();
	// UpdateRefMesh();
}

//Constructor
DFNFace::DFNFace(TPZGeoEl *gel, DFNFracture *fracture) :
    fGeoEl(gel),
    fFracture(fracture),
    fCoord(0)
{
	// DebugStop();
    int nsides = gel->NSides();
    fStatus.Resize(nsides,0);
    int nedges = gel->NCornerNodes();
    fRibs.Resize(nedges,nullptr);
}

/// Copy constructor
DFNFace::DFNFace(const DFNFace &copy){
    this->operator=(copy);
}

/// Assignment operator
DFNFace &DFNFace::operator=(const DFNFace &copy){
    fGeoEl = copy.fGeoEl;
	fStatus = copy.fStatus;
	fCoord = copy.fCoord;
	fIntersectionIndex = copy.fIntersectionIndex;
	fFracture = copy.fFracture;
	fRibs = copy.fRibs;
    fRefMesh = copy.fRefMesh;
    return *this;
}



DFNRib* DFNFace::Rib(int i) const {return fRibs[i];}



//@todo I'm not sure how robust this is. I haven't considered some 
// possibilities that may arise from fractures smaller than elements... gotta run some tests
bool DFNFace::IsOnBoundary(){
	int n_intersections = 0;
	int n_lowdim_sides = fGeoEl->NSides()-1;
	for(int i=0; i<n_lowdim_sides; i++){
		n_intersections += fStatus[i];
		if(n_intersections > 1) return false;
	}
	n_intersections = 0;
	for(int irib=0; irib < fRibs.size(); irib++){
		DFNRib* rib = fRibs[irib];
		if(!rib) continue;
		for(int i=0; i<3; i++){
			n_intersections++;
			if(n_intersections > 1) return false;
		}
	}
	return n_intersections == 1;
}

bool DFNFace::NeedsRefinement() const{
	int nsides = fGeoEl->NSides();
	if(fStatus[nsides-1]) return true;
	int nnodes = fGeoEl->NCornerNodes();
	int nedges = nnodes;
	int cutedges = 0;
	int cutnodes = 0;
	bool consecutive_nodes = false;
	// Basically return False if only one node or only two consecutive nodes have been intersected
	for(int i=0; i<nedges; i++){
		// Count number of nodes toward the intersection was snapped
		cutnodes += fStatus[i];
		// Count number of non-snapped intersections
		cutedges += fStatus[i+nnodes];
		// Check if snaps were made down to 2 consecutive nodes
		if(fStatus[i]&&fStatus[(i+1)%nnodes]){
			// check anterior and posterior edges for intersection
			consecutive_nodes = true;
			// consecutive_nodes = (fStatus[i-1+nnodes]==0)&&(fStatus[(i+1)%nnodes+nnodes]==0);
			// //                             anterior edge                   posterior edge
		}
	}
	if(cutnodes > 2){ consecutive_nodes = false; DebugStop();}

	if(consecutive_nodes) return false;
	if(cutnodes == 1 && cutedges == 0) return false;
	return true;
}

int DFNFace::CheckAssimilatedSide() const{
	int nnodes = fGeoEl->NSides(0);
	for(int i=0; i<nnodes; i++){
		if(fStatus[i] && fStatus[ (i+1)%nnodes ]) return i + nnodes;
	}
	return -1;
}

void DFNFace::UpdateRefMesh(){
	if(!this->NeedsRefinement()) return;
	TPZManVector<TPZManVector<int64_t,4>,6> child(0); 
	TPZManVector<TPZManVector<REAL,3>> newnode(0);
	FillChildrenAndNewNodes(child,newnode);
    
#ifdef PZDEBUG
	for(int in=0; in<newnode.size(); in++){
		if(newnode[in].size() != 3){
			DebugStop();
		}
	}
#endif
    
    
	int ncornersfather = fGeoEl->NCornerNodes();
	int nchildren = child.size();
	int matid = fGeoEl->MaterialId();
	TPZGeoMesh *gmesh = fGeoEl->Mesh();

	
	// Number of nodes for refinement pattern
	int refNNodes = ncornersfather + newnode.size();

	// insert nodes
	fRefMesh.CleanUp();
	fRefMesh.NodeVec().Resize(refNNodes);
	int i;
	for(i=0; i<ncornersfather; i++){
		TPZManVector<REAL,3> coord(3);
		fGeoEl->NodePtr(i)->GetCoordinates(coord);
		fRefMesh.NodeVec()[i].Initialize(coord,fRefMesh);
	}
	i = ncornersfather;
	for(auto nodeco : newnode){
		fRefMesh.NodeVec()[i].Initialize(nodeco,fRefMesh);
		i++;
	}

	// insert father
	{
		MElementType elemtype = fGeoEl->Type();
		TPZManVector<int64_t,4> cornerindices(ncornersfather);
		for(int i = 0; i<ncornersfather; i++) cornerindices[i] = i;
		int64_t index = 0;
		fRefMesh.CreateGeoElement(elemtype, cornerindices, matid, index);
	}
	// insert children
	for (int i = 0; i < nchildren; i++)
	{
		int ncorners = child[i].size();
		MElementType elemtype;
		switch (ncorners){
			case 3: elemtype = ETriangle; break;
			case 4: elemtype = EQuadrilateral; break;
		}
		int64_t index = i+1;
		// int64_t index = i;
		TPZManVector<int64_t,4> refchild(ncorners);
		for(int k = 0; k<ncorners; k++) refchild[k] = child[i][k];
		fRefMesh.CreateGeoElement(elemtype, refchild, matid, index);
	}
	fRefMesh.BuildConnectivity(); 
}

void DFNFace::FillChildrenAndNewNodes(
	TPZManVector<TPZManVector<int64_t,4>,6> &child, 
	TPZManVector<TPZManVector<REAL,3>> &newnode)
{
	if(!this->NeedsRefinement()) return;
    // this method will determine how a face needs to be refined
    // most complex operation of fracture insertion!!
	int splitcase = GetSplitPattern(fStatus);
	int nodeA = fGeoEl->NCornerNodes();
	int nodeB = fGeoEl->NCornerNodes()+1;
	// Vector of face's node indices
	TPZManVector<int64_t,4> node;
	fGeoEl->GetNodeIndices(node);
	
	switch(splitcase){
		case 1:{ // A quadrilateral with 2 opposite refined edges
			child.resize(2);
			newnode.resize(2);
			int i = 0;
			while(fStatus[i+4]==false){i++;}
			// Intersection node at rib i
			DFNRib *ribA = fRibs[i];
			newnode[0] = ribA->AntCoord();
			// Intersection node at rib opposite to rib i
			DFNRib *ribB = fRibs[i+2];
			newnode[1] = ribB->AntCoord();
			
			child[0].Resize(4);
				child[0][0] = nodeB;
				child[0][1] = nodeA;
				child[0][2] = i+1;
				child[0][3] = i+2;
			child[1].resize(4);
				child[1][0] = nodeA;
				child[1][1] = nodeB;
				child[1][2] = (i+3)%4;
				child[1][3] = i;
			break;}
		case 2:{ // A quadrilateral with 2 adjancent refined edges
			child.resize(4);
			newnode.resize(2);
			int i = 0;
			while(fStatus[i+4]==false || fStatus[(i+3)%4+4]==false){i++;}
			// Intersection node at rib clockwise adjacent to rib i
			DFNRib *ribA = fRibs[(i+3)%4];
			newnode[0] = ribA->AntCoord();
			// Intersection node at rib i
			DFNRib *ribB = fRibs[i];
			newnode[1] = ribB->AntCoord();

			child[0].Resize(3);
				child[0][0] = nodeB;
				child[0][1] = nodeA;
				child[0][2] = i;
			child[1].resize(3);
				child[1][0] = nodeB;
				child[1][1] = (i+1)%4;
				child[1][2] = (i+2)%4;
			child[2].resize(3);
				child[2][0] = nodeA;
				child[2][1] = nodeB;
				child[2][2] = (i+2)%4;
			child[3].resize(3);
				child[3][0] = nodeA;
				child[3][1] = (i+2)%4;
				child[3][2] = (i+3)%4;
			break;}
		case 3:{ // A quadrilateral with 2 opposite intersected corners
			child.resize(2);
			int i = fStatus[1];
			child[0].Resize(3);
				child[0][0] = i;
				child[0][1] = i+1;
				child[0][2] = i+2;
			child[1].resize(3);
				child[1][0] = i;
				child[1][1] = i+2;
				child[1][2] = (i+3)%4;
			break;}
		case 4: // A quadrilateral with a single refined edge
		case 7:{ // A quadrilateral with an intersected corner and a refined edge
			child.resize(3);
			newnode.resize(1);
			int i = 0;
			while(fStatus[i+4]==false) i++;
			// Intersection node at rib i
			DFNRib *ribA = fRibs[i];
			newnode[0] = ribA->AntCoord();

			child[0].resize(3);
				child[0][0] = nodeA;
				child[0][1] = (i+1)%4;
				child[0][2] = (i+2)%4;
			child[1].resize(3);
				child[1][0] = nodeA;
				child[1][1] = (i+2)%4;
				child[1][2] = (i+3)%4;
			child[2].resize(3);
				child[2][0] = nodeA;
				child[2][1] = (i+3)%4;
				child[2][2] = i;
			break;}
		case 5:{ // A quadrilateral with a refined edge and a mid-face intersection node
			DebugStop();
			// child.resize(5);
			// newnode.resize(2);
			// int i = 0;
			// while(fStatus[i+4]==false){i++;}
			// // Intersection node at rib i
			// DFNRib *ribA = fRibs[i];
			// newnode[0] = ribA->AntCoord();
			// // In-plane itersection node
			// newnode[1] = fCoord;

			// child[0].Resize(3);
			// // @warning: if child[0][0] != internal-node for cases 5, 6 and 8, DFNFace::Refine() will break
			// 	child[0][0] = nodeB;
			// 	child[0][1] = nodeA;
			// 	child[0][2] = (i+1)%4;
			// child[1].resize(3);
			// 	child[1][0] = nodeB;
			// 	child[1][1] = (i+1)%4;
			// 	child[1][2] = (i+2)%4;
			// child[2].resize(3);
			// 	child[2][0] = nodeB;
			// 	child[2][1] = (i+2)%4;
			// 	child[2][2] = (i+3)%4;
			// child[3].resize(3);
			// 	child[3][0] = nodeB;
			// 	child[3][1] = (i+3)%4;
			// 	child[3][2] = i;
			// child[4].resize(3);
			// 	child[4][0] = nodeB;
			// 	child[4][1] = i;
			// 	child[4][2] = nodeA;
			break;}
		case 6: // A quadrilateral with an intersected corner and a mid-face intersection node
		// algorithm for case 8 == case 6
		// algorithm for case 7 == case 4
		case 8:{ // A quadrilateral with a single mid-face intersection node 
			DebugStop();
			/// @deprecated
			// @warning: if child[0][0] != internal-node for cases 5, 6 and 8, DFNFace::Refine() will break
			// In-plane itersection node
			// child.resize(4);
			// newnode.resize(1);
			// newnode[0] = fCoord;
			// for(int i=0; i<4; i++) {
			// 	child[i].resize(3);
			// 	child[i][0] = nodeA;
			// 	child[i][1] = i;
			// 	child[i][2] = (i+1)%4;
			// }
			break;}
		case 9:{ // A quadrilateral with 2 mid-face intersection nodes
			/// @deprecated
			DebugStop();
			break;}
		case 10:{// A triangle with 2 adjancent refined edges
			child.resize(2);
			newnode.resize(2);
			int i = 0;
			while(fStatus[i+3]==false || fStatus[(i+2)%3+3]==false){i++;}
			// Intersection node at rib right before (counter-clockwise) to rib i
			DFNRib *ribA = fRibs[(i+2)%3];
			newnode[0] = ribA->AntCoord();
			// Intersection node at rib opposite to rib i
			DFNRib *ribB = fRibs[i];
			newnode[1] = ribB->AntCoord();
			
			child[0].Resize(3);
				child[0][0] = nodeB;
				child[0][1] = nodeA;
				child[0][2] = i;
			child[1].resize(4);
				child[1][0] = nodeA;
				child[1][1] = nodeB;
				child[1][2] = (i+1)%3;
				child[1][3] = (i+2)%3;
			break;}
		case 11: // A triangle with an intersected corner and a refined edge
		// algorithm for case 11 == case 14
		case 14:{ // A triangle with a single refined edge
			child.resize(2);
			newnode.resize(1);
			int i=0;
			while(fStatus[i+3] == false) i++;
			// Intersection node at rib i
			DFNRib *ribA = fRibs[i];
			newnode[0] = ribA->AntCoord(); 

			child[0].resize(3);
				child[0][0] = nodeA;
				child[0][1] = (i+1)%3;
				child[0][2] = (i+2)%3;
			child[1].resize(3);
				child[1][0] = nodeA;
				child[1][1] = (i+2)%3;
				child[1][2] = i;
			break;}
		case 12:{ // A triangle with a refined edge and a mid-face intersection node
			DebugStop();
			// child.resize(4);
			// newnode.resize(2);
			// int i = 0;
			// while(fStatus[i+3]==false){i++;}
			// // Intersection node at rib i
			// DFNRib *ribA = fRibs[i];
			// newnode[0] = ribA->AntCoord();
			// // In-plane itersection node
			// newnode[1] = fCoord;

			// child[0].Resize(3);
			// 	child[0][0] = nodeB;
			// 	child[0][1] = nodeA;
			// 	child[0][2] = (i+1)%3;
			// child[1].resize(3);
			// 	child[1][0] = nodeB;
			// 	child[1][1] = (i+1)%3;
			// 	child[1][2] = (i+2)%3;
			// child[2].resize(3);
			// 	child[2][0] = nodeB;
			// 	child[2][1] = (i+2)%3;
			// 	child[2][2] = i;
			// child[3].resize(3);
			// 	child[3][0] = nodeB;
			// 	child[3][1] = i;
			// 	child[3][2] = nodeA;
			break;}
		case 13: // A triangle with an intersected corner and a mid-face intersection node
		// algorithm for case 13 == case 15
		// algorithm for case 14 == case 11
		case 15:{ // A triangle with a single mid-face intersection node
			DebugStop();
			// In-plane itersection node
			// child.resize(3);
			// newnode.resize(1);
			// for(int i=0; i<3; i++) {
			// 	child[i].resize(3);
			// 	child[i][0] = nodeA;
			// 	child[i][1] = i;
			// 	child[i][2] = (i+1)%3;
			// }
			break;}
		case 16:{ // A triangle with 2 mid-face intersection nodes
			DebugStop();
			// @ToDo
			break;}
		default: DebugStop();
	}
}







/// Check if should be refined and generate the subelements of material id matID
void DFNFace::Refine(){
	if(!this->NeedsRefinement()) return;
	// if(fGeoEl->MaterialId()!=DFNMaterial::Efracture) fGeoEl->SetMaterialId(DFNMaterial::Erefined);
	// define refPattern
	TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(fRefMesh);
	fGeoEl->SetRefPattern(refpat);
	TPZManVector<TPZGeoEl*,6> children(fRefMesh.NElements()-1);
	// refine
	fGeoEl->Divide(children);

	// let children inherit polyhedral indices
	/** @note I thought this would be done more efficiently if we updated all polyhedra at the same moment, so I wrote the method DFNMesh::InheritPolyhedra() (no parameters). 
	I'm trying to avoid excessive dynamic memmory allocation in the vector DFNMesh::fPolyh_per_face, which would be increased every time a face gets refined 
	and its subelements inherit the polyhedral index. 
	This could also be approached with a chunk vector, but I've also tried and failed to employ PZChunkVector in these classes. Maybe in the future.*/
	// DFNMesh* dfn = fFracture->dfnMesh();
	// dfn->InheritPolyhedra(fGeoEl);

	// update internal node index
	int splitcase = GetSplitPattern(fStatus);
	switch(splitcase){
		case 5:
		case 6:
		case 8:	fIntersectionIndex = children[0]->NodeIndex(0);
		// @warning: if child[0][0] != internal-node for cases 5, 6 and 8, DFNFace::Refine() will break
		default: break;
	}
}





// void DFNFace::InheritPolyhedra(){
// 	if(!fGeoEl->HasSubElement()) return;
// 	DFNMesh* dfn = fFracture->dfnMesh();
// 	if(dfn->Dimension() < 3) return;
// 	dfn->

// 	int polyh_front = dfn->GetPolyhedralIndex({fGeoEl->Index(),+1});
// 	int polyh_back  = dfn->GetPolyhedralIndex({fGeoEl->Index(),-1});

// 	TPZGeoEl* child = nullptr;
// 	int nchildren = fGeoEl->NSubElements();
// 	for(int i=0; i<nchildren; i++){
// 		child = fGeoEl->SubElement(i);
// 		dfn->SetPolyhedralIndex({child->Index(),+1},polyh_front);
// 		dfn->SetPolyhedralIndex({child->Index(),-1},polyh_back);
// 	}
// }










/**
 * @brief Returns the split pattern that should be used to split this face
 * @param Status vector (boolean) that indicates which ribs and/or nodes are cut
 * @return Integer that indicates which split pattern to use. (check documentation)
 */
int DFNFace::GetSplitPattern(const TPZManVector<int> &status) const{
	if(!this->NeedsRefinement()) return 0;
    // Count number of ribs and nodes cut
    int ribscut = 0;
	int nodescut = 0;
	int n = fGeoEl->NCornerNodes();
	for (int i = 0; i < n; i++){
		ribscut += status[i+n];
		nodescut += status[i];
	}
	bool inplane_point = fStatus[fGeoEl->NSides()-1];

	// Get split case
    // Check documentation for consult on what each splitcase means
    
    // @TODO phil : where is the documentation
	int splitcase;
	if (n == 4){ //quadrilateral
		switch(ribscut){
			case 0:
				switch(nodescut){
					case 0: splitcase = 8;break; //or 9
					case 1: splitcase = 6;break;
					case 2: splitcase = 3;break; 
				}
				break;
			case 1:
				switch(nodescut){
					case 0: if(inplane_point){splitcase = 5;break;}
							else{splitcase = 4;break;}
					case 1: splitcase = 7;break;
				}
				break;
			case 2:{
				if((status[n+0]&&status[n+2])||(status[n+1]&&status[n+3])){splitcase =1;}
				else {splitcase = 2;}
				break;}
			default: DebugStop();break;
		}
	}
	else{ //triangle
		switch(ribscut){
			case 0:
				switch(nodescut){
					case 0: splitcase = 15;break; //or 16
					case 1: splitcase = 13;break;
				}
				break;
			case 1:
				switch(nodescut){
					case 0: if(inplane_point){splitcase = 12;break;}
							else{splitcase = 14;break;}
					case 1: splitcase = 11;break;
				}
				break;
			case 2:
				splitcase = 10;
				break;
			default: DebugStop();break;
		}
	}
	return splitcase;
}

// Naive vector comparison
template <class T, int NumExtAlloc1, int NumExtAlloc2>
bool operator!=(TPZManVector<T,NumExtAlloc1>& v1,TPZManVector<T,NumExtAlloc2>& v2){
	int64_t size1 = v1.size();
	int64_t size2 = v2.size();
	if(size1 != size2) return true;
	for(int64_t i = 0; i<size1; i++){
		if(v1[i]!=v2[i]) return true;
	}
	return false;
}
// // fill vector with zeroes
// template< class T, int NumExtAlloc >
// void Zero(TPZManVector<T,NumExtAlloc> &vec){
// 	int64_t size = vec.size();
// 	for(int i=0; i<size; i++){vec[i]=0;}
// }


bool DFNFace::UpdateStatusVec(){
	if(fRibs.size()<1) DebugStop();
	
	int nnodes = fGeoEl->NCornerNodes();
	int nribs = fGeoEl->NSides(1);
	DFNRib *rib;
	
	TPZManVector<int> old_fStatus = fStatus;
    fStatus.Fill(0);

	// loop through ribs and match their intersection side index into DFNFace::fStatus
	for(int i=0; i<nribs; i++){
		int irib = i;
		rib = fRibs[irib];
		if(!rib) continue;
		switch(rib->IntersectionSide()){													// Get index of side through which fracture passes (before or after snap)
			// If through side 2, simply mark that side in the status vector
			case 2: fStatus[i+nnodes] = 1; break;											
			// If it was snapped to a node, use rib orientation to permute accordingly
			case 1:																			
			case 0:{ /* fStatus[(i+(!orientation+rib->IntersectionSide())%2)%nnodes] = 1; break; */
					// Check if orientation of rib matches orientation of the faceside it occupies. If node0 of the i-th rib matches the i-th node of the face, then their orientation match.
					int orientation = rib->GeoEl()->NodeIndex(0) == fGeoEl->NodeIndex(i);
					int permuted_rib_node = (rib->IntersectionSide() + !orientation) % 2;	// Note the ! (not operator) that swaps orientation before sum
					int face_node_index = (permuted_rib_node + irib) % nnodes;				// This line is basically the opposite of TPZGeoEl::SideNodeLocIndex()
					int old = (i+(!orientation+rib->IntersectionSide())%2)%nnodes;
					fStatus[face_node_index] = 1;
					break;
			}
			default: DebugStop();
		}
	}
	bool changed = old_fStatus != fStatus;
#ifdef LOG4CXX
	// @todo This feels too verbose. Maybe make it logger->IsTraceEnabled()
	if(changed && logger->isDebugEnabled()){
		std::stringstream sstream;
		sstream << "Updated StatusVec in Face #" << fGeoEl->Index();
		sstream << "\nOld [" << old_fStatus << "]";
		DFN::SketchStatusVec(old_fStatus,sstream);
		sstream << "\nNew [" << fStatus << "]";
		DFN::SketchStatusVec(fStatus,sstream);
		LOG4CXX_DEBUG(logger,sstream.str());
	}
#endif // LOG4CXX
	return changed;
}




bool DFNFace::FindInPlanePoint(){
    // Convert TPZGeoEl into DFNPolygon
    TPZGeoEl *gelface = fGeoEl;
    TPZFMatrix<REAL> corners(3,4);
    int n;
    // check if face is quadrilateral
    if(gelface->Type() == MElementType::ETriangle){
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
        DFNPolygon faceplane(subcorners, fGeoEl->Mesh());
        // Check fPolygon's ribs for intersection with faceplane
        int nribs = fFracture->Polygon().GetCornersX().Cols();
        for(int irib = 0; irib < nribs; irib++){
            TPZManVector<REAL,3> p1(3);
            TPZManVector<REAL,3> p2(3);
            for(int i = 0; i<3; i++){
                p1[i] = fFracture->Polygon().GetCornersX().g(i, irib);
                p2[i] = fFracture->Polygon().GetCornersX().g(i, (irib+1)%nribs);
            }
            if(faceplane.Check_pair(p1, p2, fCoord)){
				fStatus[fGeoEl->NSides()-1] = 1;
                return true;
            }
        }
    }
    return false;
}

void DFNFace::UpdateMaterial(){
	DebugStop(); // deprecated
	if(fGeoEl->MaterialId() == DFNMaterial::Efracture) return;

	bool two_in_a_row = false;
	int cutnodes = 0;
	int cutedges = 0;
	int nnodes = fGeoEl->NCornerNodes();
	for(int i=0; i<nnodes; i++){
		cutnodes += fStatus[i];
		if(fStatus[i]+fStatus[(i+1)%nnodes]==2){two_in_a_row = true;}
	}

	if(cutnodes==nnodes){fGeoEl->SetMaterialId(DFNMaterial::Efracture);}
	else if(cutnodes==2 && two_in_a_row){fGeoEl->SetMaterialId(DFNMaterial::Eintact);}
	else {fGeoEl->SetMaterialId(DFNMaterial::Erefined);}
}

int64_t DFNFace::LineInFace() const{
	// if(!fGeoEl->HasSubElement()) return -1;
	int nsides = fGeoEl->NSides();
	int nintersections=0;
	TPZStack<int> sides;
	for(int i=0; i<nsides; i++){
		if(fStatus[i]) sides.push_back(i);
		nintersections += fStatus[i];
	}
	if(nintersections < 2) return -1;
	if(nintersections > 2) DebugStop(); // @todo.. as of 16/july/2020 we're still decided on this not happening. So I'll have to come back here if we change our minds

	int nnodes = fGeoEl->NCornerNodes();
	TPZStack<int64_t> framenodes;
	for(int i=0; i<nsides; i++){
		if(fStatus[i]){
			if(i<nnodes){
				framenodes.push_back(fGeoEl->NodeIndex(i));
			}else if(i<nsides-1){
				framenodes.push_back(this->Rib(i-nnodes)->IntersectionIndex());
			}else if(i == nsides-1){
				framenodes.push_back(this->IntersectionIndex());
			}
		}
	}
	int nchildren = fGeoEl->NSubElements();
	if(nchildren == 0){nchildren++;}
	// queue all possible lines by checking 1D neighbours of children
	std::set<TPZGeoEl *> candidate_ribs;
	for(int ichild=0; ichild<nchildren; ichild++){
		TPZGeoEl* child;
		if(fGeoEl->HasSubElement()){
			child = fGeoEl->SubElement(ichild);
		}else{
			child = fGeoEl;
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
	if(candidate_ribs.size() == 0) return -1;
	if(candidate_ribs.size() == 1) return (*candidate_ribs.begin())->Index();
	int nframenodes = framenodes.size();
	for(auto candidate : candidate_ribs){
		for(int inode=0; inode<nframenodes; inode++){
			if(framenodes[inode]==candidate->NodeIndex(0)){
				for(int nextnode=0; nextnode<nframenodes; nextnode++){
					if(framenodes[nextnode]==candidate->NodeIndex(1)){
						return candidate->Index();
					}
				}
			}
		}
	}
	DebugStop();
	return -1;
}

// method returns true is a side has been divided in two
bool DFNFace::CanBeSnapped() const{
	const int ncorners = fGeoEl->NCornerNodes();
	const int nsides = fGeoEl->NSides();
	// If fStatus indicates an intersection that wasn't snapped, then it could potentially be snapped, so it returns true.
	// If all intersections were already snapped, there's nothing else to snap and method returns false.
	for(int iside = ncorners; iside<nsides; iside++){
		if(fStatus[iside]) break;			// if this condition evaluates to true, there's an intersection that was not snapped yet.
		if(iside==nsides-1) return false;
	}
	return true;
}





bool DFNFace::AllSnapsWentToSameNode() const{
	const int nsides = fGeoEl->NSides();
	int sumAll = 0;
	for(int i=0; i<nsides; i++){
		sumAll += fStatus[i];
	}
	int sumNodes = fStatus[0]+fStatus[1]+fStatus[2]+(fGeoEl->NCornerNodes()==4)*fStatus[3];
	return (sumNodes==1 && sumAll==1);
}




bool DFNFace::NeedsSnap(REAL tolDist, REAL tolAngle_cos){
    // Consistency checks
	if(!this->NeedsRefinement()) {return false;}
	int nels = fRefMesh.NElements();
	if(nels < 1) {PZError<<"\n\n"<<__PRETTY_FUNCTION__<<"\nSnap face refinements requires a proper refinement mesh to describe the refinement\n";DebugStop();}
	if(tolAngle_cos > 1 + gDFN_SmallNumber || tolAngle_cos < -1 - gDFN_SmallNumber) {
		PZError << "\n\n Tolerable angle cosine doesn't seem to be a cosine value.\n"; DebugStop();
	}

	// If intersections have already been snapped to corners, there's nothing else to snap
	if(!this->CanBeSnapped()) return false;

    // this stretch of code determines if an element has a bad angle or aspect ratio
	for(int iel=1; iel <nels; iel++){
		TPZGeoEl* gel = fRefMesh.Element(iel);
		// Check if any angle violates the tolerable angle
		int ncorners = gel->NCornerNodes();
		for(int icorner=0; icorner < ncorners; icorner++){
			REAL corner_cosine = DFN::CornerAngle_cos(gel,icorner);
			if(corner_cosine > tolAngle_cos) return true;
		}
		// Check if any edge violates the tolerable distance
		int nedges = gel->NSides(1);
		int nsides = gel->NSides();
		for(int iedge=ncorners; iedge<nsides-1; iedge++){
			if(gel->SideArea(iedge) < tolDist) return true;
		}
	}
	return false;
}


bool DFNFace::SnapIntersection_try(REAL tolDist, REAL tolAngle_cos){
    // this method determines if the mesh has elements with bad aspect ratio
	// while(this->NeedsSnap(tolDist, tolAngle_cos))
	if(this->NeedsSnap(tolDist, tolAngle_cos))
    {
        // @pedro : how are you sure that "all" sides need to be snapped?
        // can it not be that a single rib needs to be snapped?
		// @reply: I agree, this was a mistake of mine, I'll rewrite the method.
        // thanks!!
		this->SnapIntersection_force();
        return true;
    }
	else
    {
		return false;
    }
}

// this method will snap the nodes of all ribs
bool DFNFace::SnapIntersection_force(){
	int nribs = fGeoEl->NSides(1);
    int ncorner = fGeoEl->NCornerNodes();
	for(int irib=0; irib<nribs; irib++){
		DFNRib* rib = fRibs[irib];
		if(!rib) continue;
        // will snap the node to a closest point
		rib->SnapIntersection_force();
        // adjust the neighbouring faces to acount for the snapped node
		UpdateNeighbours(irib+ncorner);
	}
	UpdateStatusVec();
	UpdateRefMesh();
	return true;
}

void DFNFace::UpdateNeighbours(int iside){
	TPZGeoElSide gelside(fGeoEl, iside);
	TPZGeoElSide neig = gelside.Neighbour();
	DFNFace* neig_face = nullptr;
	for(/*void*/;neig!=gelside; neig = neig.Neighbour()){
		if(neig.Element()->Dimension() != 2) continue;
		//@todo skeletonMesh material would enter here
		neig_face = fFracture->Face(neig.Element()->Index());
		if(!neig_face) continue;
        // this probably changes the values of the undocumented fStatus vector
		if(!neig_face->UpdateStatusVec()) continue;
		neig_face->UpdateRefMesh();
		REAL tolDist = fFracture->dfnMesh()->TolDist();
		REAL tolAngle_cos = fFracture->dfnMesh()->TolAngle_cos();
		neig_face->SnapIntersection_try(tolDist,tolAngle_cos);
		// neig_face->SnapIntersection_force();
	}
}

void DFNFace::AddRib(DFNRib* ribptr, int side){
	if(fGeoEl->SideDimension(side) != 1) DebugStop();
	TPZGeoElSide ribside = {ribptr->GeoEl(),2};
	TPZGeoElSide faceside = {fGeoEl,side};
	if(!faceside.NeighbourExists(ribside)) DebugStop();

	int ribindex = side - fGeoEl->NSides(0);
	fRibs[ribindex] = ribptr;
}

/// Print the data structure
void DFNFace::Print(std::ostream &out, bool print_refmesh) const
{
    out<<"\nFace GeoEl index # "<<fGeoEl->Index();
	out<<"\nIntersected ribs:\n";
	for(auto rib : fRibs){
		out<<"\t";
		if(!rib) out<<"---";
		else out<<rib->Index();
		out<<"\n";
	}

	out<<"Status Vector : {"<< fStatus<<"}";
	out<<"\nSplit Case: "; 
	int splitcase = GetSplitPattern(fStatus);
	if(!splitcase){	
		out<<"\'unrefined\'";
	}else{  		
		out<<splitcase;
	}

	out<<"\nCorner Nodes:\n";
	for(int i=0; i<fGeoEl->NCornerNodes(); i++){
		out<<"\t"<<i<<": index: "<<fGeoEl->NodeIndex(i)<<" "; 
		fGeoEl->Node(i).Print(out);
	}

	SketchStatusVec(out);

	if(print_refmesh){
		// int nels = fGeoEl->Mesh()->NElements();
		// int w = (int) std::log10(nels)+1;
		out << "[Start][Refinement Mesh][Face Index " 
			// << std::setw(w) 
			// << std::right 
			<< fGeoEl->Index()
			<< "]\n\n";
		fRefMesh.Print(out);
		out << "[End][Refinement Mesh][Face Index " 
			// << std::setw(w) 
			// << std::right 
			<< fGeoEl->Index()
			<< "]\n\n";
	}
	out<<"\n";
}



/// @brief Returns the other 1D side with a DFNRib
int DFNFace::OtherRibSide(int inletside) const{
	if(fGeoEl->SideDimension(inletside) != 1) DebugStop();
	int nnodes = fGeoEl->NSides(0);
	int inletedge = inletside-nnodes;
	if(!fRibs[inletedge]) DebugStop();

	int nedges = fGeoEl->NSides(1);
	for(int outedge=0; outedge < nedges; outedge++){
		if(outedge == inletedge) continue;
		if(fRibs[outedge]) return outedge + nnodes;
	}
	#ifdef PZDEBUG
	fFracture->dfnMesh()->DumpVTK();
	#endif // PZDEBUG
	DebugStop();
	return -1;
}

int DFNFace::NIntersectedRibs() const{
	int n=0;
	for(int i=0,nribs=fGeoEl->NSides(1); i<nribs; i++) {
		n += fRibs[i]!=nullptr;
	}
	return n;
}
int DFNFace::NInboundRibs() const{
	int n=0;
	for(int i=0,nribs=fGeoEl->NSides(1); i<nribs; i++) {
		n += fRibs[i]!=nullptr && !fRibs[i]->IsOffbound();
	}
	return n;
}



void DFNFace::SketchStatusVec(std::ostream& out) const{
	std::stringstream sout;
	sout << "\nStatusVec Sketch:";
	switch(fGeoEl->Type()){
		case MElementType::ETriangle:{
			sout << "\n  " << (fStatus[2]?"x":" ");
			sout << "\n | \\";
			sout << "\n" << (fStatus[5]?"x":" ") <<"|  \\" << (fStatus[4]?"x":" ");
			sout << "\n |___\\";
			sout << "\n" << (fStatus[0]?"x":" ") << "  " << (fStatus[3]?"x":" ") << "  " << (fStatus[1]?"x":" ");
			break;
		}
		case MElementType::EQuadrilateral:{
			sout << "\n" << (fStatus[3]?"x":" ") << " __" << (fStatus[6]?"x":"_") << "__ " << (fStatus[2]?"x":" ");
			sout << "\n |     | ";
			sout << "\n" << (fStatus[7]?"x":" ") << "|     |" << (fStatus[5]?"x":" ");
			sout << "\n |_____| ";
			sout << "\n" << (fStatus[0]?"x":" ") << "   " << (fStatus[4]?"x":" ") << "   " << (fStatus[1]?"x":" ");
			break;
		}
		default: DebugStop();
	}

	out << sout.str();
}