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

    /// empty constructor and destructor
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
    // /**
    //  * @brief Get the next face forward or backwards according to direction and fill the angle between faces
    //  * @param current_index = index of current element
    //  * @param angle to be filled with angle between current and next face
    //  * @param direction 1 or -1
    // */
    // TRolodexCard& NextCard(int64_t current_index, REAL& angle, int direction=1){
    //     TRolodexCard& current_card = Card(current_index);
    //     int ncards = fcards.size();
    //     int jcard = current_card.fposition; //position where current card appears in the card vector
        
    //     int64_t next_id = -1;
    //     switch(direction){
    //         case  1: next_id = (jcard+1)%ncards;        break;
    //         case -1: next_id = (jcard-1+ncards)%ncards; break;
    //         default: DebugStop();
    //     }
    //     TRolodexCard& next_card = fcards[next_id];
    //     angle = next_card.fangle_to_reference - current_card.fangle_to_reference;
    //     if(angle < 0.) {angle = angle + DFN::_2PI;}
    //     return next_card;
    // }
    /**
     * @brief Get the next face forward or backwards according to direction and fill the angle between faces
     * @param current_index = index of current element
     * @param angle to be filled with angle between current and next face
     * @param direction 1 or -1
    */
    std::pair<TRolodexCard,int> FacingCard(std::pair<TRolodexCard,int> input, REAL& angle){
        int jcard = input.first.fposition;
        int ncards = fcards.size();
        int facing_id = -1;
        int direction = input.second*input.first.forientation;
        switch(direction){
            case  1: facing_id = (jcard+1)%ncards;        break;
            case -1: facing_id = (jcard-1+ncards)%ncards; break;
            default: DebugStop();
        }
            
        std::pair<TRolodexCard,int> output;
        output.first = fcards[facing_id];
        output.second = -direction*output.first.forientation;
        // compute angle between faces
        angle = direction*(output.first.fangle_to_reference - input.first.fangle_to_reference);
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
};



inline std::ostream& operator<<(std::ostream &out, const TRolodex& rolodex){
    rolodex.Print(out);
    return out;
}





#endif /* DFNRolodex_h */