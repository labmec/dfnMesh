/*! 
 *	TRSFace.cpp
 *  @authors   Pedro Lima
 *  @authors   Jorge Ordoñez
 *  @date      2018-2019
 */

#include "TRSFace.h"

// Empty constructor
TRSFace::TRSFace(){
    fFaceIndex = -1;
    fIsCut = false;
}

//Constructor
TRSFace::TRSFace(int64_t index, bool iscut){
    fFaceIndex = index;
    fIsCut = iscut;
}

/// Copy constructor
TRSFace::TRSFace(const TRSFace &copy){
    fFaceIndex = copy.fFaceIndex;
    fIsCut = copy.fIsCut;
    fStatus = copy.fStatus;
    fRibs = copy.fRibs;
    fSubFaces = copy.fSubFaces;
    fIntersection = copy.fIntersection;
}

/// Assignment operator
TRSFace &TRSFace::operator=(const TRSFace &copy){
    fFaceIndex = copy.fFaceIndex;
    fIsCut = copy.fIsCut;
    fStatus = copy.fStatus;
    fRibs = copy.fRibs;
    fSubFaces = copy.fSubFaces;
    fIntersection = copy.fIntersection;
    return *this;
}

/// Define the element index and whether it cuts the plane
void TRSFace::SetElementIndex(int64_t elindex, bool iscut){
    fFaceIndex = elindex;
    fIsCut = iscut;
}

/// Element index
int64_t TRSFace::ElementIndex() const{
    return fFaceIndex;
}

/// Intersects the plane or not
bool TRSFace::IsCut() const{
    return fIsCut;
}

/// Return the subelement indices
TPZVec<int64_t> TRSFace::SubElements() const{
    return fSubFaces;
}

/// Set the subelement indices
void TRSFace::SetChildren(const TPZVec<int64_t> &subels){
    fSubFaces = subels;
}











/// Divide the given this surface and generate the subelements
void TRSFace::DivideSurface(TRSFractureMesh *fracmesh, int matid){
	TPZGeoMesh *gmesh = fracmesh->GetgeoMesh();
   	TPZGeoEl *face = gmesh->Element(fFaceIndex);

	// check if face exist in mesh
   	if(!face){DebugStop();}
	
	TPZVec<int64_t> node;
	face->GetNodeIndices(node);
	
	// Determine pattern of refinement
	int splitcase = GetSplitPattern(fStatus);

	// Vector of children elements
	TPZVec<TPZVec<int64_t>> child;

	// Divide surface according to split pattern (these algorithms make no sense without documentation)
	switch(splitcase){
		case 1:{
			child.resize(2);
			int i = 0;
			while(fStatus[i+4]==false){i++;}
			// Intersection node at rib i
			int64_t nodeA = gmesh->Element(fracmesh->Rib(fRibs[i])->IntersectionIndex())->NodeIndex(0);
			// Intersection node at rib opposite to rib i
			int64_t nodeB = gmesh->Element(fracmesh->Rib(fRibs[i+2])->IntersectionIndex())->NodeIndex(0);
			
			TPZVec<int64_t> child1(4);
				child1[0] = nodeB;
				child1[1] = nodeA;
				child1[2] = node[i+1];
				child1[3] = node[i+2];
			TPZVec<int64_t> child2(4);
				child2[0] = nodeA;
				child2[1] = nodeB;
				child2[2] = node[(i+3)%4];
				child2[3] = node[i];
			
			child[0] = child1;
			child[1] = child2;
			break;}
		case 2:{
			child.resize(4);
			int i = 0;
			while(fStatus[i+4]==false || fStatus[(i+3)%4+4]==false){i++;}
			// Intersection node at rib i
			int64_t nodeA = gmesh->Element(fracmesh->Rib(fRibs[(i+3)%4])->IntersectionIndex())->NodeIndex(0);
			// Intersection node at rib opposite to rib i
			int64_t nodeB = gmesh->Element(fracmesh->Rib(fRibs[i])->IntersectionIndex())->NodeIndex(0);

			TPZVec<int64_t> child1(3);
				child1[0] = nodeA;
				child1[1] = nodeB;
				child1[2] = node[i];
			TPZVec<int64_t> child2(3);
				child2[0] = nodeB;
				child2[1] = node[(i+1)%4];
				child2[2] = node[(i+2)%4];
			TPZVec<int64_t> child3(3);
				child3[0] = nodeA;
				child3[1] = nodeB;
				child3[2] = node[(i+2)%4];
			TPZVec<int64_t> child4(3);
				child4[0] = nodeA;
				child4[1] = node[(i+2)%4];
				child4[2] = node[(i+3)%4];

			child[0] = child1;
			child[1] = child2;
			child[2] = child3;
			child[3] = child4;
			break;}
		case 3:{
			break;}
		case 4:{
			break;}
		case 5:{
			child.resize(5);
			int i = 0;
			while(fStatus[i+4]==false){i++;}
			// Intersection node at rib i
			int64_t nodeA = gmesh->Element(fracmesh->Rib(fRibs[i])->IntersectionIndex())->NodeIndex(0);
			// In-plane itersection node
			int64_t nodeB = gmesh->Element(fIntersection)->NodeIndex(0);

			TPZVec<int64_t> child1(3);
				child1[0] = nodeB;
				child1[1] = nodeA;
				child1[2] = node[(i+1)%4];
			TPZVec<int64_t> child2(3);
				child2[0] = nodeB;
				child2[1] = node[(i+1)%4];
				child2[2] = node[(i+2)%4];
			TPZVec<int64_t> child3(3);
				child3[0] = nodeB;
				child3[1] = node[(i+2)%4];
				child3[2] = node[(i+3)%4];
			TPZVec<int64_t> child4(3);
				child4[0] = nodeB;
				child4[1] = node[(i+3)%4];
				child4[2] = node[i];
			TPZVec<int64_t> child5(3);
				child5[0] = nodeB;
				child5[1] = node[i];
				child5[2] = nodeA;

			child[0] = child1;
			child[1] = child2;
			child[2] = child3;
			child[3] = child4;
			child[4] = child5;
			break;}
		case 6:{
			break;}
		case 7:{
			break;}
		case 8:{
			break;}
		case 9:{
			break;}
		case 10:{
			break;}
		case 11:{
			break;}
		case 12:{
			break;}
		case 13:{
			break;}
		case 14:{
			break;}
		case 15:{
			break;}
		case 16:{
			break;}
		default: DebugStop();
	}
}
 















/**
 * @brief Returns the split pattern that should be used to split this face
 * @param Status vector (boolean) that indicates which ribs and/or nodes are cut
 * @return Integer that indicates which split pattern to use. (check documentation)
 */
int TRSFace::GetSplitPattern(TPZVec<bool> &status){
    // Count number of ribs and nodes cut
    int ribscut = 0;
	int nodescut = 0;
	int n = (int) status.size()/2;
	for (int i = 0; i < n; i++){
		ribscut += status[i+n];
		nodescut += status[i];
	}

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
					case 0: if(fIntersection == -1){splitcase = 4;break;}
							else{splitcase = 5;break;}
					case 1: splitcase = 7;break;
				}
				break;
			case 2:{
				int i=0;
				while(status[i+4]==false){i++;}
				if(status[i+5]==true||status[(i+3)%4+4]==true){splitcase = 2;}
				else splitcase = 1;

				// std::vector<bool> test1 = {0,0,0,0,1,0,1,0};
				// if(status == test1){splitcase = 1;break;}
				// std::vector<bool> test2 = {0,0,0,0,0,1,0,1};
				// if(status == test2){splitcase = 1;break;}
				// splitcase = 2;
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
					case 0: if(fIntersection == -1){splitcase = 14;break;}
							else{splitcase = 12;break;}
					case 1: splitcase = 11;break;
				}
				break;
			case 2:
				splitcase = 10;
				break;
		}
	}
	return splitcase;
}
  