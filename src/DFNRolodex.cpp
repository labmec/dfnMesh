#include "DFNRolodex.h"

#if PZ_LOG
    static TPZLogger logger("dfn.mesh");
#endif

#ifdef PZDEBUG
    void TRolodex::CoherentAngleTest(const int64_t AxleIndex, const std::map<REAL, TPZGeoElSide>& SortedCards){
        constexpr REAL _359degrees = 359.*M_PI/180.;
        const REAL MaxAngle = SortedCards.rbegin()->first;
        if(MaxAngle > _359degrees){
            std::stringstream msg;
            msg << "[FATAL] Something went wrong when computing angles to sort faces around edge #" << AxleIndex
                << "\nAngles should go from zero to 2pi. And these were the angles you asked me to initialize:\n"
                << "gel index | side | angle\n";
            for(auto& card_it : SortedCards){
                msg << std::setw(9) << std::right << card_it.second.Element()->Index() << " | "
                    << std::setw(4) << std::right << card_it.second.Side() << " | "
                    << std::left  << card_it.first << "\n";
            }
            msg << "\nMaybe an orientation match problem."
                << "\nMaybe Gmsh created a sliver."
                << "\nI'll try to plot the rolodex to a VTK file.";
            LOGPZ_FATAL(logger,msg.str());
            
            constexpr char plotpath[] = "./LOG/twoDNeighbours.vtk";
            TPZGeoMesh* gmesh = SortedCards.begin()->second.Element()->Mesh();
            const TPZGeoElSide edgeside(gmesh->Element(AxleIndex),2);
            TRolodex::PlotVTK(plotpath,edgeside,true,true);
            LOGPZ_FATAL(logger,"\n2D neighbours of edge #" << AxleIndex << " plotted to: " << plotpath << '\n');
            
            constexpr char _3Dplotpath[] = "./LOG/threeDNeighbours.vtk";
            DFN::PlotNeighbours(_3Dplotpath,edgeside,3,true,false);
            LOGPZ_FATAL(logger,"\n3D neighbours of edge #" << AxleIndex << " plotted to: " << _3Dplotpath << '\n');
            DebugStop();
        }
    }
#endif // PZDEBUG




void TRolodex::PlotVTK(const std::string filepath, const TPZGeoElSide Axle, const bool UnrefinedOnly, const bool orientationMatch){
    // Consistency
    if(!Axle.Element()) DebugStop();
    if(Axle.Dimension() != 1) DebugStop();

    // Initialize graphics mesh
    TPZGeoMesh* gmesh = Axle.Element()->Mesh();
    TPZGeoMesh graphicMesh;
    graphicMesh.ElementVec().Resize(gmesh->NElements());
    for (int i = 0; i < graphicMesh.NElements(); i++) {
        graphicMesh.ElementVec()[i] = nullptr;
    }
    graphicMesh.NodeVec() = gmesh->NodeVec();

    // Copy elements to graphics mesh
    Axle.Element()->Clone(graphicMesh)->SetMaterialId(1);
    TPZGeoElSide neig = Axle.Neighbour();
    for(/*void*/; neig != Axle; neig++){
        if(neig.Element()->Dimension() > 2) continue;
        if(neig.Element()->Dimension() < 1) continue;
        if(UnrefinedOnly && neig.Element()->HasSubElement()) continue;
        
        TPZGeoEl* newel = neig.Element()->Clone(graphicMesh);
        if(newel->HasSubElement()) newel->SetSubElement(0,nullptr);
        if(orientationMatch){
            int orientation = DFN::OrientationMatch(Axle,neig);
            newel->SetMaterialId(orientation);
        }
    }

    // choose a reference
    TPZGeoElSide reference{nullptr,-1};
    int reference_orientation = 0;
    if(Axle.Element()->Dimension() == 2 ){
        reference = Axle;
        reference_orientation = 1;
    }else{
        neig = Axle.Neighbour();
        for(/*void*/; neig != Axle; neig++){
            if(neig.Element()->Dimension() != 2) continue;
            if(UnrefinedOnly && neig.Element()->HasSubElement()) continue;
            reference = neig;
            reference_orientation = (DFN::OrientationMatch(Axle,reference)?1:-1);
            break;
        }
        if(!reference.Element()){
            if(UnrefinedOnly){LOGPZ_FATAL(logger,"\nPlot failed. May you meant to plot a rolodex with faces that were refined. Change your 'UnrefinedOnly' argument.")}
            LOGPZ_FATAL(logger,"\nI don't know how to create and plot a Rolodex for an Axle that does not have any 2D neighbour.\n")
            DebugStop();
        }
    }
    // Compute angle to reference
    TPZVec<REAL> elData(gmesh->NElements(),-9.0);
    elData[Axle.Element()->Index()] = 0.0;
    neig = Axle.Neighbour();
    for(/*void*/; neig != Axle; neig++){
        if(neig.Element()->Dimension() != 2) continue;
        if(UnrefinedOnly && neig.Element()->HasSubElement()) continue;
        
        const REAL angle = DFN::DihedralAngle(reference,neig,reference_orientation);
        elData[neig.Element()->Index()] = angle;
    }
    
    std::ofstream out(filepath);
    TPZVTKGeoMesh::PrintGMeshVTK(&graphicMesh, out, elData);
}