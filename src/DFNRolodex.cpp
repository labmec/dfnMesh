#include "DFNRolodex.h"


#ifdef PZDEBUG
    void TRolodex::CoherentAngleTest(const int64_t AxleIndex, const std::map<REAL, TPZGeoElSide>& SortedCards){
        constexpr REAL _359degrees = 359.*M_PI/180.;
        const REAL MaxAngle = SortedCards.rbegin()->first;
        if(MaxAngle > _359degrees){
            PZError << "[FATAL] Something went wrong when computing angles to sort faces around edge #" << AxleIndex
                    << "\nAngles should go from zero to 2pi. And these were the angles you asked me to initialize:\n"
                    << "gel index | side | angle";
            for(auto& card_it : SortedCards){
                PZError << std::setw(9) << std::right << card_it.second.Element()->Index() << " | "
                        << std::setw(4) << std::right << card_it.second.Side() << " | "
                        << std::left  << card_it.first << "\n";
            }
            PZError << "\nMaybe an orientation match problem."
                    << "\nI'll try to plot the rolodex to a VTK file.";
            constexpr char plotpath[] = "./LOG/twoDNeighbours.vtk";
            TPZGeoMesh* gmesh = SortedCards.begin()->second.Element()->Mesh();
            TPZGeoElSide edgeside(gmesh->Element(AxleIndex),2);
            DFN::PlotNeighbours(plotpath,edgeside,2,true,true);
            PZError << "\n2D neighbours of edge #" << AxleIndex << " plotted to: " << plotpath << std::endl;
            DebugStop();
        }
    }
#endif // PZDEBUG

