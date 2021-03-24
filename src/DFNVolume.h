/*! 
 *  @brief     Describes a volume within DFN scope of application.
 *  @details   
 *  @author    Pedro Lima
 *  @date      2019
 *  @deprecated WE ARE NOT USING THIS CLASS. IT WAS CONCEIVED AS A TRACKER FOR COARSE ELEMENT AND SUB-MESH
 */
// #ifndef DFNVolume_h
// #define DFNVolume_h

// #include <stdio.h>
// #include "pzmanvector.h"
// #include "pzmatrix.h"

// /*! 
//  * @brief      Describes a volume within the DFN scope of application.
//  *  @details    Has information on the planes it encloses and the children volumes
//  *  @author     Pedro Lima
//  *  @date       2019
//  */
// class DFNVolume
// {
//   private:
// 	/// Volume GeoElement index within it's gmesh
//     int64_t fVolumeIndex;

//     /// Indicates whether it is cut by fracture plane
//     bool fIsCut;

//     /// Vector of sub elements
//     TPZManVector<int64_t> fSubEls;
    
//     /// Index of in-volume intersection EPoint
//     int64_t fIntersection = -1;
    
//     // @pedro there is no such thing as a "plane" in the datastructure, do you mean 2d geometric
//     // elements that form the volume boundary
//     /// Indices of planes enclosed by volume
//     TPZManVector<int64_t> fEnclosedFaces;

//     // /// Submesh
//     // TPZGeoMesh *submesh;

//   public:
//     /// Empty constructor
//     DFNVolume(){DebugStop();/*[DEPRECATED]*/}
//     /// Empty destructor
//     ~DFNVolume(){};
    
//     /// Constructor
//     DFNVolume(int64_t index, bool iscut);
    
//     /// Copy constructor
//     DFNVolume(const DFNVolume &copy);
    
//     /// Assignment operator
//     DFNVolume &operator=(const DFNVolume &copy);
    
//     /// Define the element index and whether it cuts the plane
//     void SetElementIndex(int64_t elindex, bool iscut);
    
//     /// Element index
//     int64_t ElementIndex() const;
    
//     /// Intersects the plane or not
//     bool IsCut() const;
    
//     /// Return the subelement indices
//     TPZVec<int64_t> SubElements() const;
    
//     /// Set the subelement indices
//     void SetChildren(const TPZVec<int64_t> &subels);

//     /// Set a 2D element enclosed by the volume
//     void SetFaceInVolume(int64_t Elindex);

//     /// Get 2D elements enclosed by volume
//     TPZVec<int64_t> GetFacesInVolume(){return fEnclosedFaces;}

//     /** @brief Print method for logging */
//     void Print(std::ostream& out = std::cout) const;
// };
// #endif /* DFNVolume_h */