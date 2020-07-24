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


//Constructor
DFNFace::DFNFace(TPZGeoEl *gel, DFNFracture *fracture) :
	fGeoEl(gel),
	fFracture(fracture),
	fCoord(0)
{
	int nsides = gel->NSides();
	fStatus.Resize(nsides,0);
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



DFNRib* DFNFace::Rib(int i){return fFracture->Rib(fRibs[i]);}



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
	for(int64_t irib : fRibs){
		DFNRib* rib = fFracture->Rib(irib);
		if(!rib) continue;
		for(int i=0; i<3; i++){
			n_intersections += rib->StatusVec()[i];
			if(n_intersections > 1) return false;
		}
	}
	return n_intersections == 1;
}





bool DFNFace::UpdateRefMesh(){
	if(!this->NeedsRefinement()) return false;
	TPZManVector<TPZManVector<int64_t,4>,6> child(0); 
	TPZManVector<TPZManVector<REAL,3>> newnode(0);
	FillChildrenAndNewNodes(child,newnode);
	int ncornersfather = fGeoEl->NCornerNodes();
	int nchildren = child.size();
	TPZGeoMesh *gmesh = fGeoEl->Mesh();

	
	// Number of nodes for refinement pattern
	int refNNodes = ncornersfather + newnode.size();

	// insert nodes
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
		fRefMesh.CreateGeoElement(elemtype, cornerindices, DFNMaterial::Erefined, index);
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
		fRefMesh.CreateGeoElement(elemtype, refchild, DFNMaterial::Erefined, index);
	}
	fRefMesh.BuildConnectivity(); 
}

void DFNFace::FillChildrenAndNewNodes(
	TPZManVector<TPZManVector<int64_t,4>,6> &child, 
	TPZManVector<TPZManVector<REAL,3>> &newnode)
{
	if(!this->NeedsRefinement()) return;
	int splitcase = GetSplitPattern(fStatus);
	int nodeA = fGeoEl->NCornerNodes();
	int nodeB = fGeoEl->NCornerNodes()+1;
	// Vector of face's node indices
	TPZManVector<int64_t,4> node;
	fGeoEl->GetNodeIndices(node);
	
	switch(splitcase){
		case 1:{
			child.resize(2);
			newnode.resize(2);
			int i = 0;
			while(fStatus[i+4]==false){i++;}
			// Intersection node at rib i
			DFNRib *ribA = fFracture->Rib(fRibs[i]);
			newnode[0] = ribA->AntCoord();
			// Intersection node at rib opposite to rib i
			DFNRib *ribB = fFracture->Rib(fRibs[i+2]);
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
		case 2:{
			child.resize(4);
			newnode.resize(2);
			int i = 0;
			while(fStatus[i+4]==false || fStatus[(i+3)%4+4]==false){i++;}
			// Intersection node at rib clockwise adjacent to rib i
			DFNRib *ribA = fFracture->Rib(fRibs[(i+3)%4]);
			newnode[0] = ribA->AntCoord();
			// Intersection node at rib i
			DFNRib *ribB = fFracture->Rib(fRibs[i]);
			newnode[1] = ribB->AntCoord();

			child[0].Resize(3);
				child[0][0] = nodeA;
				child[0][1] = nodeB;
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
		case 3:{
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
		case 4:
		case 7:{
			child.resize(3);
			newnode.resize(1);
			int i = 0;
			while(fStatus[i+4]==false) i++;
			// Intersection node at rib i
			DFNRib *ribA = fFracture->Rib(fRibs[i]);
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
		case 5:{
			child.resize(5);
			newnode.resize(2);
			int i = 0;
			while(fStatus[i+4]==false){i++;}
			// Intersection node at rib i
			DFNRib *ribA = fFracture->Rib(fRibs[i]);
			newnode[0] = ribA->AntCoord();
			// In-plane itersection node
			newnode[1] = fCoord;

			child[0].Resize(3);
			// @warning: if child[0][0] != internal-node for cases 5, 6 and 8, DFNFace::Refine() will break
				child[0][0] = nodeB;
				child[0][1] = nodeA;
				child[0][2] = (i+1)%4;
			child[1].resize(3);
				child[1][0] = nodeB;
				child[1][1] = (i+1)%4;
				child[1][2] = (i+2)%4;
			child[2].resize(3);
				child[2][0] = nodeB;
				child[2][1] = (i+2)%4;
				child[2][2] = (i+3)%4;
			child[3].resize(3);
				child[3][0] = nodeB;
				child[3][1] = (i+3)%4;
				child[3][2] = i;
			child[4].resize(3);
				child[4][0] = nodeB;
				child[4][1] = i;
				child[4][2] = nodeA;
			break;}
		case 6:
		// case 7 == case 4
		case 8:{ // case 8 == case 6
			// @warning: if child[0][0] != internal-node for cases 5, 6 and 8, DFNFace::Refine() will break
			// In-plane itersection node
			child.resize(4);
			newnode.resize(1);
			newnode[0] = fCoord;
			for(int i=0; i<4; i++) {
				child[i].resize(3);
				child[i][0] = nodeA;
				child[i][1] = i;
				child[i][2] = (i+1)%4;
			}
			break;}
		case 9:{
			// @ToDo
			DebugStop();
			break;}
		case 10:{
			child.resize(2);
			newnode.resize(2);
			int i = 0;
			while(fStatus[i+3]==false || fStatus[(i+2)%3+3]==false){i++;}
			// Intersection node at rib right before (counter-clockwise) to rib i
			DFNRib *ribA = fFracture->Rib(fRibs[(i+2)%3]);
			newnode[0] = ribA->AntCoord();
			// Intersection node at rib opposite to rib i
			DFNRib *ribB = fFracture->Rib(fRibs[i]);
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
		case 11:
		case 14:{ //case 11 == case 14
			child.resize(2);
			newnode.resize(1);
			int i=0;
			while(fStatus[i+3] == false) i++;
			// Intersection node at rib i
			DFNRib *ribA = fFracture->Rib(fRibs[i]);
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
		case 12:{
			child.resize(4);
			newnode.resize(2);
			int i = 0;
			while(fStatus[i+3]==false){i++;}
			// Intersection node at rib i
			DFNRib *ribA = fFracture->Rib(fRibs[i]);
			newnode[0] = ribA->AntCoord();
			// In-plane itersection node
			newnode[1] = fCoord;

			child[0].Resize(3);
				child[0][0] = nodeB;
				child[0][1] = nodeA;
				child[0][2] = (i+1)%3;
			child[1].resize(3);
				child[1][0] = nodeB;
				child[1][1] = (i+1)%3;
				child[1][2] = (i+2)%3;
			child[2].resize(3);
				child[2][0] = nodeB;
				child[2][1] = (i+2)%3;
				child[2][2] = i;
			child[3].resize(3);
				child[3][0] = nodeB;
				child[3][1] = i;
				child[3][2] = nodeA;
			break;}
		case 13: // case 13 == case 15
		// case 14 == case 11
		case 15:{
			// In-plane itersection node
			child.resize(3);
			newnode.resize(1);
			for(int i=0; i<3; i++) {
				child[i].resize(3);
				child[i][0] = nodeA;
				child[i][1] = i;
				child[i][2] = (i+1)%3;
			}
			break;}
		case 16:{
			// @ToDo
			break;}
		default: DebugStop();
	}
}







/// Check if should be refined and generate the subelements of material id matID
void DFNFace::Refine(){
	if(!this->NeedsRefinement()) return;
	fGeoEl->SetMaterialId(DFNMaterial::Erefined);
	// define refPattern
	TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(fRefMesh);
	fGeoEl->SetRefPattern(refpat);
	TPZManVector<TPZGeoEl*,6> children(fRefMesh.NElements()-1);
	// refine
	fGeoEl->Divide(children);
	// set father for children?
	// create sskeleton?

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
















/**
 * @brief Returns the split pattern that should be used to split this face
 * @param Status vector (boolean) that indicates which ribs and/or nodes are cut
 * @return Integer that indicates which split pattern to use. (check documentation)
 */
int DFNFace::GetSplitPattern(TPZManVector<int> &status){
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
// fill vector with zeroes
template< class T, int NumExtAlloc >
void Zero(TPZManVector<T,NumExtAlloc> &vec){
	int64_t size = vec.size();
	for(int i=0; i<size; i++){vec[i]=0;}
}


bool DFNFace::UpdateStatusVec(){
	if(fRibs.size()<1 || fRibs[0]<0) DebugStop();
	
	TPZGeoMesh *gmesh = fGeoEl->Mesh();
	int nnodes = fGeoEl->NCornerNodes();
	int nribs = nnodes;
	DFNRib *rib;
	TPZGeoEl *rib_gel;
	
	TPZManVector<int> old_fStatus = fStatus;
	Zero(fStatus);

	int orientation = 0;
	for(int i=0; i<nribs; i++){
		rib = fFracture->Rib(fRibs[i]);
		if(rib){
			rib_gel = gmesh->Element(fRibs[i]);
			orientation = (rib_gel->NodeIndex(0) == fGeoEl->NodeIndex(i)) ? 0 : 1;
			fStatus[i] = MAX(fStatus[i],rib->StatusVec()[orientation]);
			fStatus[(i+1)%nnodes] = MAX(fStatus[(i+1)%nnodes],rib->StatusVec()[(orientation+1)%2]);
			fStatus[i+nnodes] = rib->StatusVec()[2];
		}
	}

	return old_fStatus != fStatus;
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
        DFNPolygon faceplane(subcorners);
        // Check fPolygon's ribs for intersection with faceplane
        int nribs = fFracture->Polygon()->GetCornersX().Cols();
        for(int irib = 0; irib < nribs; irib++){
            TPZManVector<REAL,3> p1(3);
            TPZManVector<REAL,3> p2(3);
            for(int i = 0; i<3; i++){
                p1[i] = fFracture->Polygon()->GetCornersX()(i, irib);
                p2[i] = fFracture->Polygon()->GetCornersX()(i, (irib+1)%nribs);
            }
            if(faceplane.Check_rib(p1, p2, &fCoord)){
				fStatus[fGeoEl->NSides()-1] = 1;
                return true;
            }
        }
    }
    return false;
}

void DFNFace::UpdateMaterial(){
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
	// @todo... this is incomplete. 2_nodes_in_a_row may mean intact
}

int64_t DFNFace::LineInFace(){
	if(!fGeoEl->HasSubElement()) return -1;
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