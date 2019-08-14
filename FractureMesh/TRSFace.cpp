/*! 
 *	TRSFace.cpp
 *  @authors   Jorge Ordoñez
 *  @authors   Pedro Lima
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

///
///Divide the given this surface and generate the subelements
void TRSFace::DivideSurface(TPZGeoMesh *gmesh, int matid){}
    //@TODO Check consistency
// 

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
				break;
			}
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

//    if(!(gmesh->Element(fFaceIndex))){DebugStop();}
//    TPZGeoEl *gel = gmesh->Element(fFaceIndex);
//    int no_intr_P1 = fRibInter[0].second;
//    int no_intr_P2 = fRibInter[1].second;
//    int P1 = fRibInter[0].first;
//    int P2 = fRibInter[1].first;
//    TPZGeoEl *gel1 = gmesh->Element(P1);
//    TPZGeoEl *gel2 = gmesh->Element(P2);
//    TPZVec<int> status(4,0);
//    for(int iside=4; iside<8; iside++){
//        TPZGeoElSide gelside(gel,iside);
//        TPZGeoElSide neig = gelside.Neighbour();
//        int index_neig = neig.Element()->Index();
//        if(index_neig==P1){
//            status[iside]=1;
//        }
//        if(index_neig==P2){
//            status[iside]=2;
//        }
//    }
//    int cases=-1;
//    if(status[0]!=0 && status[2]!=0){
//        cases=0;
//    }
//    if(status[1]!=0 && status[3]!=0){
//         cases=1;
//    }
//    if (cases < 0){
//        
//    }
    
  
    
    
