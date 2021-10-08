/*! 
 *	DFNRolodex.hpp
 *  @authors    Pedro Lima
 *  @authors    Philippe Devloo
 *  @date       2020
 *  @brief      Contains the implementation of structs DFNRolodex and DFNRolodexCard
 */

#ifndef DFNRolodex_h
#define DFNRolodex_h

#include "DFNNamespace.h"






// A 2D element sorted around an edge (a rolodex card)
struct TRolodexCard{
  public: // Data structure
    /// index of the geometric 2d element/side represented by this card
    int64_t fgelindex = -1;
    int fSide = -1;
    /// angle measured around the one dimensional rib
    REAL fangle_to_reference = 0.0;
    /// integer indicating whether the card with next higher angle is in the direction of the
    // normal to the face or contrary
    int forientation = 0; //1 or -1
    /// position of this card in its rolodex
    int fposition = -1;

  public: // Member functions
    /// empty constructor and destructor
    // TRolodexCard(){};
    TRolodexCard(int64_t id=-1, REAL angle=0.0, int side=0)
            :fgelindex(id),fangle_to_reference(angle), fSide(side){};
    ~TRolodexCard(){};
    /// copy constructor and attribution operator
    TRolodexCard& operator=(const TRolodexCard& copy){
        fangle_to_reference = copy.fangle_to_reference;
        fgelindex = copy.fgelindex;
        fSide = copy.fSide;
        forientation = copy.forientation;
        fposition = copy.fposition;
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
        switch(forientation){
            case +1: out<<"  (+)  ";break;
            case -1: out<<"  (-)  ";break;
            default: out<<"   "<<forientation<<"   ";
        }
        out<<" | ";
        // out << std::setw(9) << std::right << forientation*fgelindex << " | ";
        out << std::setw(9) << std::right << fgelindex           << " | ";
        out << std::setw(4) << std::right << fSide               << " | ";
        out                 << std::left  << fangle_to_reference << "\n";
    }
};










/// A set of 2D elements around an edge, sorted by an angle to a reference
// The reference is the first card = fcards.begin() and angles are right-hand oriented according to the orientation of the edge.
struct TRolodex{

public:
    // a vector of cards in this rolodex
    std::vector<TRolodexCard> fcards;
    /// index of the one dimensional geometric element (the axle of the rolodex)
    int64_t fedgeindex;

public:
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

#ifdef PZDEBUG
    /// Check if angles used to initialize a rolodex are coherent
    static void CoherentAngleTest(const int64_t AxleIndex, const std::map<REAL, TPZGeoElSide>& SortedCards);
#endif // PZDEBUG

    /// @brief Constructor from an std::map of sorted cards {angle, face_edgeside}
    void Initialize(const int64_t AxleIndex, const std::map<REAL,TPZGeoElSide>& SortedCards){
#ifdef PZDEBUG
        TRolodex::CoherentAngleTest(AxleIndex, SortedCards);
#endif // PZDEBUG

        fedgeindex = AxleIndex;
		fcards.clear();
		fcards.resize(SortedCards.size());

        TPZGeoEl* AxleGel = SortedCards.begin()->second.Element()->Mesh()->Element(AxleIndex);
        TPZGeoElSide edgeside(AxleGel,2);

        int j = 0;
		for(auto& iterator : SortedCards){
			TPZGeoElSide faceside = iterator.second;
			TRolodexCard& facecard = fcards[j];
			facecard.fgelindex = faceside.Element()->Index();
			facecard.fSide = faceside.Side();
			facecard.fangle_to_reference = iterator.first;
			facecard.forientation = (DFN::OrientationMatch(edgeside,faceside)?1:-1);
			facecard.fposition = j;
			j++;
		}
    }

    /**
     * @brief Get the next face forward or backwards according to direction and fill the angle between faces
     * @param card_orientation a pair of a card in the rolodex and the direction of the desired next card (backwards -1 or forwards +1)
     * @param angle (output) to be filled with angle between current and next card
    */
    std::pair<TRolodexCard,int> FacingCard(std::pair<TRolodexCard,int> card_orientation, REAL& angle){
        int jcard = card_orientation.first.fposition;
        int ncards = fcards.size();
        int facing_id = -1;
        int direction = card_orientation.second*card_orientation.first.forientation;
        switch(direction){
            case  1: facing_id = (jcard+1)%ncards;        break;
            case -1: facing_id = (jcard-1+ncards)%ncards; break;
            default: DebugStop();
        }
            
        std::pair<TRolodexCard,int> output;
        output.first = fcards[facing_id];
        output.second = -direction*output.first.forientation;
        // compute angle between faces
        angle = direction*(output.first.fangle_to_reference - card_orientation.first.fangle_to_reference);
        if(angle < 0.) {angle = angle + DFN::_2PI;}
        
        return output;
    }

    /** @brief Print method for logging */
    void Print(std::ostream& out = std::cout) const{
        out << "\nRolodex around edge # "<< fedgeindex<<"\n";
        out << "orient. | gel index | side | angle\n";
        for(auto& card : fcards){
            card.Print(out);
        }
    }

    /**
     * @brief Get a reference to a card using an index
     * @param index of the geometric element of that card
     * @param position to fill with the position where this card appears in the card vector
     * @note This method involves a search (std::find) through a vector
    */
    TRolodexCard& Card(int64_t index){
        TRolodexCard dummycard(index,0.0); // just to use std::find
        auto it = std::find(fcards.begin(),
                            fcards.end(),
                            dummycard);
        if(it == fcards.end()){
            PZError << "\nYou're looking for a card that doesn't exist in this rolodex.\nCard index : " << index; 
            this->Print(PZError); 
            DebugStop();
        }
        return *it;
    }


    /** @return Number of cards*/
    int NCards() const{return fcards.size();}

    /// Plot VTK for unitialized Rolodex
    static void PlotVTK(const std::string filepath, const TPZGeoElSide Axle, const bool UnrefinedOnly = true, const bool orientationMatch = false);
    
    
    void PlotVTK(const std::string filepath, TPZGeoMesh* gmesh){
    
        std::cout << "\n =====> Plotting rolodex to " << filepath << std::endl;
        
        TPZGeoMesh auxMesh;
        auxMesh.ElementVec().Resize(gmesh->NElements());
        for (int i = 0; i < auxMesh.NElements(); i++) {
            auxMesh.ElementVec()[i] = nullptr;
        }
        auxMesh.NodeVec() = gmesh->NodeVec();
        

        // Axle
        TPZGeoEl* axle = gmesh->Element(fedgeindex);
        if (!axle) {
            DebugStop();
        }
        TPZGeoEl* newaxle = axle->Clone(auxMesh);
        newaxle->SetMaterialId(0);
        TPZVec<REAL> elData(gmesh->NElements(),-9.0);
        elData[newaxle->Index()] = 0.0;
        
        // All the preexisting cards in rolodex
        for (auto& card : fcards) {
            const int index = card.fgelindex;
            if (index < 0) {
                DebugStop();
            }
            TPZGeoEl *gelInRolodex = gmesh->Element(index);
            if (!gelInRolodex) {
                DebugStop();
            }
            TPZGeoEl* newel = gelInRolodex->Clone(auxMesh);
            newel->SetMaterialId(card.forientation);
            elData[newel->Index()] = card.fangle_to_reference;
        }
        
        std::ofstream out(filepath);
        TPZVTKGeoMesh::PrintGMeshVTK(&auxMesh, out, elData);
        // TPZVTKGeoMesh::PrintGMeshVTK(&auxMesh, out, true, true);

    }
};



inline std::ostream& operator<<(std::ostream &out, const TRolodex& rolodex){
    rolodex.Print(out);
    return out;
}

#endif /* DFNRolodex_h */
