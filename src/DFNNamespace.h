/*! 
 *  @authors   Pedro Lima
 *  @date      2020.09
 */

#ifndef DFNNamespace_h
#define DFNNamespace_h


#include "pzgmesh.h"



/**
 * @brief Contains methods and variables of the namespace DFN
*/
namespace DFN{
// 2*3.1415...
const float _2PI = 6.2831853071795865;
// A small number for geometric tolerances
static const double gSmallNumber = 1.e-3;



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
 * @param sideorientation decides if method returns angle or 2*pi-angle, as a right-handed axis 
 * orientation, with the right thumb place over the shared 1D side, and considering the first element node distribution.
 * If thumb orientation matches the orientation of gelside, use 1, else, use -1.
 */
static float DihedralAngle(TPZGeoElSide &gelside, TPZGeoElSide &neighbour, int sideorientation = 1){

    // Consistency checks
    if(gelside.Element()->Dimension() != 2)     DebugStop();
    if(gelside.Dimension() != 1)                DebugStop();
    if(neighbour.Element()->Dimension() !=2)    DebugStop();
    if(neighbour.Dimension() != 1)              DebugStop();

    if(gelside == neighbour) return 0.;
    // {
    //     switch(sideorientation){
    //         case  1: return 0.;
    //         case -1: return DFN::_2PI;
    //     }
    // }

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
    float angle = atan2(y,x);
    
    return (angle >= 0? angle : angle + _2PI);
}

/**
 * @brief Get a vector from node 0 to node 1 of a 1D side
 */
static void GetSideVector(TPZGeoElSide &gelside, TPZManVector<REAL,3>& vector){
    if(gelside.Dimension() != 1) DebugStop();
    int node0 = gelside.SideNodeLocIndex(0);
    int node1 = gelside.SideNodeLocIndex(1);
    
    TPZManVector<REAL,3> coord0(3,0);
    TPZManVector<REAL,3> coord1(3,0);
    gelside.Element()->Node(node0).GetCoordinates(coord0);
    gelside.Element()->Node(node1).GetCoordinates(coord1);
    
    vector = coord1 - coord0;
}

/**
 * @brief Check if the side that connects 2 neighbours has the same orientation in each element
 * @note currently exclusive to 1D sides
 */
static bool OrientationMatch(TPZGeoElSide &neig1, TPZGeoElSide &neig2){
    if(neig1.Dimension() != 1) DebugStop();
    if(!neig1.NeighbourExists(neig2)) DebugStop();
    return (neig1.SideNodeIndex(0) == neig2.SideNodeIndex(0));
}

/**
 * @brief Computes the cossine of the angle at a corner of a 2D element
*/
static REAL CornerAngle_cos(TPZGeoEl *gel, int corner){
	int ncorners = gel->NCornerNodes();
	if(corner >= ncorners) DebugStop();

	int nsides = gel->NSides();

	TPZManVector<REAL,3> point_corner(3);
	TPZManVector<REAL,3> point_anterior(3);
	TPZManVector<REAL,3> point_posterior(3);

	gel->Node(corner).GetCoordinates(point_corner);
	gel->Node((corner+1)%ncorners).GetCoordinates(point_posterior);
	gel->Node((corner-1+ncorners)%ncorners).GetCoordinates(point_anterior);

	TPZManVector<REAL,3> vec1 = point_posterior - point_corner;
	TPZManVector<REAL,3> vec2 = point_anterior  - point_corner;

	REAL cosine = DotProduct_f(vec1,vec2)/(Norm_f(vec1)*Norm_f(vec2));
	return cosine;
}

    /**
     * @brief Get a pointer to an element that is superposed in a lower dimensional side of a geometric element
     * @return nullptr if there is no element in that side
    */
    static TPZGeoEl* GetSkeletonNeighbour(TPZGeoEl* gel, int side){
        if(gel->SideDimension(side) == gel->Dimension()) return nullptr;
        TPZGeoElSide gelside(gel,side);
        TPZGeoElSide neig;
        for(neig = gelside.Neighbour(); neig!=gelside; neig=neig.Neighbour()){
            if(neig.Element()->Dimension() == gel->SideDimension(side)){
                return neig.Element();
            }
        }
        return nullptr;
    }
}






    // A 2D element sorted around an edge (a rolodex card)
    struct TRolodexCard{
        /// index of the geometric 2d element/side represented by this card
        int64_t fgelindex;
        int fSide;
        /// angle measured around the one dimensional rib
        REAL fangle_to_reference;
        /// integer indicating whether the card with next higher angle is in the direction of the
        // normal to the face or contrary
        int forientation = 0; //1 or -1

        /// empty constructor and destructor
        TRolodexCard(int64_t id=-1, REAL angle=0.0, int side=0):fgelindex(id),fangle_to_reference(angle), fSide(side){};
        ~TRolodexCard(){};
        /// copy constructor and attribution operator
        TRolodexCard& operator=(const TRolodexCard& copy){
            fangle_to_reference = copy.fangle_to_reference;
            fgelindex = copy.fgelindex;
            fSide = copy.fSide;
            forientation = copy.forientation;
            return *this;
        }
        TRolodexCard(const TRolodexCard& copy){
            this->operator=(copy);
        }
        /// less than operator
        bool operator<(const TRolodexCard& other){
            return fgelindex < other.fgelindex;
        }
        bool operator==(const TRolodexCard& other){
            return fgelindex == other.fgelindex;
        }
        /** @brief Print method for logging */
        void Print(std::ostream& out = std::cout) const{
            out <<" "<< fgelindex<< " | "<< fSide << " | " << fangle_to_reference << "\n";
        }
    };
    /// A set of 2D elements around an edge, sorted by an angle to a reference
    // The reference is the first card = fcards.begin()
    struct TRolodex{
        std::vector<TRolodexCard> fcards;
        /// index of the one dimensional geometric index
        int64_t fedgeindex;

        /// default constructor and destructor
        TRolodex(): fedgeindex(-1){fcards.resize(0);};
        ~TRolodex(){};
        /// copy constructor and attribution operator
        TRolodex& operator=(const TRolodex& copy){
            fcards = copy.fcards;
            fedgeindex = copy.fedgeindex;
            return *this;
        }
        TRolodex(const TRolodex& copy){
            this->operator=(copy);
        }
        /**
         * @brief Get the next face forward or backwards according to direction and fill the angle between faces
         * @param current_index = index of current element
         * @param angle to be filled with angle between current and next face
         * @param direction 1 or -1
        */
        TRolodexCard& NextFace(int64_t current_index, REAL& angle, int direction=1){
            int jcard; //position where current card appears in the card vector
            TRolodexCard& current_card = Card(current_index,jcard);
            int ncards = fcards.size();
            
            int64_t next_id = -1;
            switch(direction){
                case  1: next_id = (jcard+1)%ncards;        break;
                case -1: next_id = (jcard-1+ncards)%ncards; break;
                default: DebugStop();
            }
            TRolodexCard& next_card = fcards[next_id];
            angle = next_card.fangle_to_reference - current_card.fangle_to_reference;
            if(angle < 0.) {angle = angle + DFN::_2PI;}
            return next_card;
        }

        /**
         * @brief Get a reference to a card using an index
         * @param index of the geometric element of that card
         * @param position to fill with the position where this card appears in the card vector
         * @note This methods involve a search (std::find) through a vector
        */
        TRolodexCard& Card(int64_t index, int & position){
            TRolodexCard dummycard(index,0.0); // just to use std::find
            auto it = std::find(fcards.begin(),
                                fcards.end(),
                                dummycard);
            if(it == fcards.end()) DebugStop();
            position = it - fcards.begin();
            return *it;
        }

        /** @brief Print method for logging */
        void Print(std::ostream& out = std::cout) const{
            out << "\n\nRolodex around edge # "<< fedgeindex<<"\n";
            out << "gel index | side index | angle\n";
            for(auto& card : fcards){
                card.Print(out);
            }
        }
    };



#endif /* DFNNamespace_h */
