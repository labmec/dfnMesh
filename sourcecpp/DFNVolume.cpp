/*! 
 *	TRSVolume.cpp
 *  @author     Pedro Lima
 *  @date       2019
 */

#include "TRSVolume.h"

// Empty constructor
TRSVolume::TRSVolume(){
    fVolumeIndex = -1;
    fIsCut = false;
}

//Constructor
TRSVolume::TRSVolume(int64_t index, bool iscut){
    fVolumeIndex = index;
    fIsCut = iscut;
}

/// Copy constructor
TRSVolume::TRSVolume(const TRSVolume &copy){
    fVolumeIndex = copy.fVolumeIndex;
    fIsCut = copy.fIsCut;
    fSubEls = copy.fSubEls;
    fIntersection = copy.fIntersection;
}

/// Assignment operator
TRSVolume &TRSVolume::operator=(const TRSVolume &copy){
    fVolumeIndex = copy.fVolumeIndex;
    fIsCut = copy.fIsCut;
    fSubEls = copy.fSubEls;
    fIntersection = copy.fIntersection;
    return *this;
}

/// Define the element index and whether it cuts the plane
void TRSVolume::SetElementIndex(int64_t elindex, bool iscut){
    fVolumeIndex = elindex;
    fIsCut = iscut;
}

/// Element index
int64_t TRSVolume::ElementIndex() const{
    return fVolumeIndex;
}

/// Intersects the plane or not
bool TRSVolume::IsCut() const{
    return fIsCut;
}

/// Return the subelement indices
TPZVec<int64_t> TRSVolume::SubElements() const{
    return fSubEls;
}

/// Set the subelement indices
void TRSVolume::SetChildren(const TPZVec<int64_t> &subels){
    fSubEls = subels;
}

/// Set a 2D element enclosed by the volume
void TRSVolume::SetFaceInVolume(int64_t Elindex){
    int N = fEnclosedFaces.size();
    fEnclosedFaces.resize(N+1);
    fEnclosedFaces[N] = Elindex;
}