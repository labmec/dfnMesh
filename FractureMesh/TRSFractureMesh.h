/*! 
 *  @authors   Pedro Lima
 *  @authors   Jorge Ordoñez
 *  @date      2018-2019
 */

#ifndef TRSFractureMesh_h
#define TRSFractureMesh_h

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include <stdio.h>
#include <tuple>
#include "pzmatrix.h"
#include "pzcompel.h"
#include "pzgeoelbc.h"

#include "TRSRibs.h"
#include "TRSFace.h"
#include "TRSVolume.h"
#include "TRSFracPlane.h"


typedef TPZFMatrix<REAL> Matrix;

/*! 
 *  @brief     Compares a geomesh with fracture plane to find intersections.
 *  @details   Intersection search is performed after creation of skeleton
 * elements with TRSFractureMesh::CreateSkeletonElements. Fracture plane should
 *  be a TRSFracPlane.
 *  @authors   Pedro Lima
 *  @authors   Jorge Ordoñez
 *  @date      2018-2019
 */
class TRSFractureMesh
{
private:
    /// Define a default tolerance
    REAL fTolerance = 1.e-4;
    
    /// Map of intersected ribs
    std::map<int64_t, TRSRibs> fRibs;
    
    /// Map of intersected faces
    std::map<int64_t, TRSFace> fMidFaces;
    
    /// Map of end-fracture faces
    std::map<int64_t, TRSFace> fEndFaces;

    /// Map of intersected volumes
    std::map<int64_t, TRSVolume> fVolumes;

    /// Pointer for the geometric mesh
    TPZGeoMesh *fGMesh;

    /// Bounded plane from a fracture
    TRSFracPlane fFracplane;

    /// Fracplane's geometric element index
    int64_t fFracplaneindex;

    /// Material of elements at fracture surface
    int fSurfaceMaterial;

public:
    
    /// Empty constructor
    TRSFractureMesh();
    
    /**Define the fracture plane from 3 to 4 points
     * Points should be coplanar
     * The matrix should be dimension 3xN, each column defining the coordinates
     * of a point
     *  
     */
    TRSFractureMesh(TRSFracPlane &FracPlane, TPZGeoMesh *gmesh, int matID);
    
    /// Copy constructor
    TRSFractureMesh(const TRSFractureMesh &copy);
    
    /// Assignment operator
    TRSFractureMesh &operator=(const TRSFractureMesh &copy);
    
    /// Associate the geometric mesh
    void SetgeoMesh(TPZGeoMesh *gmesh){
        fGMesh=gmesh;
    }
    
    /// Access the geomesh
    TPZGeoMesh* GetGeoMesh(){
        return fGMesh;
    }
    
    /// Return the corner nodes of the fracture
    TRSFracPlane GetPlane() const;
    
    /// Modify the default tolerance
    void SetTolerance(REAL tolerance);
    REAL GetTolerance() const;
    
private:
        
    /// Checks neighbour's dimension and returns true if it is equal
    bool HasEqualDimensionNeighbour(TPZGeoElSide &gelside);
    
    /// Finds intersection point of fracture boundaries and geometric mesh faces
    TPZVec<REAL> FindEndFracturePoint(TRSFace &face);
    
public:
    
    /// Insert intersection elements of lower dimension in the geometric mesh.
    void CreateSkeletonElements(int dimension, int matid);
    
    /// Access the ribs data structure
    void AddRib(TRSRibs rib);
    
    /// Access mid-fracture faces' data structure
    void AddMidFace(TRSFace &face);
    
    /// Access end-fracture faces' data structure
    void AddEndFace(TRSFace &face);

    /// Insert new volume in data structure
    void AddVolume(TRSVolume volume);

    /// Pointer to rib of index 'index'
    TRSRibs *Rib(int64_t index){return &fRibs[index];}

    /// Pointer to face of index 'index'
    TRSFace *Face(int64_t index);
    
    TRSVolume *Volume(int64_t index){return &fVolumes[index];}
    /// Find and split intersected faces
    void SplitFaces(int matID);
    
    /// Find and split intersected ribs
    void SplitRibs(int matID);

    /// Connects fracture-edge intersections (temporary name for lack of better one)
    void SplitFractureEdge();

    /// Triangulates fracture plane
    void SplitFracturePlane();

    /// Write mesh elements to .geo file
    void WriteGMSH(TPZGeoMesh *pzgmesh);

    /// Uses Gmsh to mesh volumes cut by fracture plane
    void CreateVolumes();

    /// Sets material for elements at surface of fracture
    void SetSurfaceMaterial(int matID);

    /// Get material ID for elements at surface of fracture
    int MaterialID(){return fSurfaceMaterial;}

    /// Find the volumetrical element that encloses a 2D element
    bool FindEnclosingVolume(TPZGeoEl *ifracface);
};

#endif /* TRSFractureMesh_h */



/*
bool TRSFractureMesh::FindEnclosingVolume(TPZGeoEl ifracface)
{
    TPZVec<REAL> faceCenter = ifracface.center;
    std::map<REAL, int64_t> candidates;
    iterate over ifracface 1D sides (iside)
    {
        iterate over neighbours through iside (ineig)
        {
            if(ineig.Dimension != 2){continue;}
            if(ineig.material == FracSurfaceMaterial){continue;}
            if(ineig.HasFather == false){continue;}

            TPZGeoEl *father = ineig.GetFather()
            TPZVec<REAL> fatherCenter = father->center;
            TPZGeoElSide *fatherSide = father->side(father->NSides - 1);

            iterate over neighbours through fatherside (ivolume)
            {
                if(ivolume.Dimension != 3){continue;}
                TPZVec<REAL> volumecenter = ivolume->center;
                TPZVec<REAL> v1 = faceCenter - fatherCenter;
                TPZVec<REAL> v2 = volumecenter - fatherCenter;

                REAL dot = DotProduct(v1,v2);
                if(dot >=0)
                {
                    dot.Normalize();
                    candidates[dot] = ivolume->Index
                }
            }
        }
    }

    if(candidates.size() > 0){
        // 'reverse iterator begin' gives biggest key in map
        int64_t volumeindex = candidates.rbegin()->second;
        volumes[volumeindex]->SetFaceInVolume(ifracface.index);
        return true;
    }

    // degeneracy: ifracface's edges are completely enclosed by volume
    
    std::set<TPZGeoElSide *> verified;
    iterate over nodes (iside)
    {
        iterate over neighbours through inode (ineig)
        {
            if(ineig.Dimension != 2){continue;}
            if(ineig.material == FracSurfaceMaterial){continue;}
            if(ineig.HasFather == false){continue;}

            TPZGeoEl *father = ineig.GetFather()
            TPZVec<REAL> fatherCenter = father->center;
            TPZGeoElSide *fatherSide = father->side(father->NSides - 1);

            iterate over neighbours through fatherside (ivolume)
            {
                if(ivolume.Dimension != 3){continue;}
                if(verified.find(ivolume)){continue;}
                TPZVec<REAL> ksi(3,2);
                bool test = ivolume->ComputeXinverse(faceCenter,qsi,fTolerance)
                if(test == true){
                    volumes[ivolume->Index]->SetFaceInVolume(ifracface.index);
                    return true;
                }
            }
        }
    }

    std::cout<<"\n TRSFractureMesh::FindEnclosingVolume found no enclosing volume for element #"<<ifracface.index<<"\n";
    return false;
}
*/