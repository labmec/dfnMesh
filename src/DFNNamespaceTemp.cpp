/*! 
 *  @authors   Pedro Lima
 *  @date      2020.10
 *  @brief     Namespace DFN template definitions
 */

#ifndef DFNNamespaceTemp_cpp
#define DFNNamespaceTemp_cpp


#include "pzgmesh.h"
#include "DFNNamespace.h"
#include <array>



/**
 * @brief Contains methods and variables of the namespace DFN
*/
namespace DFN{



    template<typename TReturnType,typename Ttype1,typename Ttype2>
    TReturnType DotProduct(TPZManVector<Ttype1,3> &vec1, TPZManVector<Ttype2,3> &vec2){
        int size1 = vec1.size();
        int size2 = vec2.size();
        if(size1 != size2){DebugStop();}
        TReturnType dot = 0.;
        for(int j=0; j<size1; j++){
            dot += vec1[j]*vec2[j];
        }
        return dot;
    }

    template<typename Ttype1, typename Ttype2>
    float DotProduct_f(TPZManVector<Ttype1,3> &vec1, TPZManVector<Ttype2,3> &vec2){
        int size1 = vec1.size();
        int size2 = vec2.size();
        if(size1 != size2){throw std::bad_exception();}
        float dot = 0.;
        for(int j=0; j<size1; j++){
            dot += vec1[j]*vec2[j];
        }
        return dot;
    }

    /** @brief Returns the norm of a vector with template precision*/
    template<typename Ttype>
    Ttype Norm(TPZManVector<Ttype, 3> &vec){
        Ttype norm = 0.;
        for(int j=0, size=vec.size(); j<size; j++){
            norm += vec[j]*vec[j];
        }
        return std::sqrt(norm);
    }

    /** @brief Returns the norm of a vector with float precision*/
    template<typename Ttype>
    float Norm_f(TPZManVector<Ttype, 3> &vec){
        float norm = 0.;
        for(int j=0, size=vec.size(); j<size; j++){
            norm += vec[j]*vec[j];
        }
        return std::sqrt(norm);
    }

    /** 
     * @brief Vector cross product with template return type
     * @param ReturnType CrossProduct<ReturnType>(vec1,vec2)
    */
    template<typename T1, typename T2>
    TPZVec<T1> CrossProduct(TPZManVector<T2,3> &vec1, TPZManVector<T2,3> &vec2){
        if(vec1.size() != 3){throw std::invalid_argument("CrossProduct requires vector of dimension 3\n");}
        if(vec2.size() != 3){throw std::invalid_argument("CrossProduct requires vector of dimension 3\n");}
        
        TPZManVector<T1,3> result(3,0.);
        result[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
        result[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
        result[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
        return result;
    }

    template<typename Ttype>
    TPZVec<float> CrossProduct_f(TPZManVector<Ttype,3> &vec1, TPZManVector<Ttype,3> &vec2){
        if(vec1.size() != 3){throw std::invalid_argument("CrossProduct requires vector of dimension 3\n");}
        if(vec2.size() != 3){throw std::invalid_argument("CrossProduct requires vector of dimension 3\n");}
        
        TPZManVector<float,3> result(3,0.);
        result[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
        result[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
        result[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
        return result;
    }
    template <typename T>
    TPZVec<T> operator-(TPZVec<T>& v1,TPZVec<T>& v2){
        int64_t size1 = v1.size();
        int64_t size2 = v2.size();
        if(size1 != size2) DebugStop();
        TPZVec<T> result(size1);
        for(int64_t i = 0; i<size1; i++){
            result[i] = v1[i] - v2[i];
        }
        return result;
    }
    template <typename T>
    TPZVec<T> operator+(TPZVec<T>& v1,TPZVec<T>& v2){
        int64_t size1 = v1.size();
        int64_t size2 = v2.size();
        if(size1 != size2) DebugStop();
        TPZVec<T> result(size1);
        for(int64_t i = 0; i<size1; i++){
            result[i] = v1[i] + v2[i];
        }
        return result;
    }

    /// @brief Get an orthogonal projection of x onto an arbitrary plane represented by a normal vector and a reference point in the plane
    /// @warning It assumes the given normal vector has norm == 1 and won't verify (to keep it efficient)
    template<typename Ttype>
    TPZManVector<Ttype, 3> GetProjectedX(TPZManVector<Ttype, 3> &x,TPZManVector<Ttype,3>& inplane_point,TPZManVector<Ttype,3>& normal){
        TPZManVector<Ttype, 3> projection(3,0.);
        TPZManVector<Ttype, 3> dist_to_inplane = x - inplane_point;
        Ttype orth_dist = DotProduct<Ttype>(dist_to_inplane,normal);
        projection[0] = x[0] - orth_dist*normal[0];
        projection[1] = x[1] - orth_dist*normal[1];
        projection[2] = x[2] - orth_dist*normal[2];
        return projection;
    }



    /** @brief Set material ids for a set of element indices and skip negative entries*/
    template<class TContainer>
    void BulkSetMaterialId(TPZGeoMesh* gmesh, TContainer elementindices, int matid){
        if(elementindices.size() < 1) DebugStop();
        for(int64_t iel : elementindices){
            if(iel < 0) continue;
            gmesh->Element(iel)->SetMaterialId(matid);
        }
    }

    /** @brief returns the intersection of 2 sets*/
    template<typename Ttype>
    std::set<Ttype> set_intersection(const std::set<Ttype>& set1,const std::set<Ttype>& set2){
        std::set<Ttype> intersection;
        auto& smaller_set =  set1.size() > set2.size() ? set2 : set1;
        auto& bigger_set = !(set1.size() > set2.size()) ? set2 : set1;

        if(smaller_set.size() < 1) return intersection; //empty set

        auto end = bigger_set.end();
        for(auto& iel : smaller_set){
            auto itr = bigger_set.find(iel);
            if(itr == end) continue;
            intersection.insert(*itr);
        }
        return intersection;
    }



    /// builds a loop of oriented 1D elements occupying the 1D sides of a 2D el
    /// @param shift: indices will get shifted by a constant 
    template<class Tcontainer>
    void GetLineLoop(TPZGeoEl* face_el, Tcontainer& lineloop, const int shift){
        if(face_el->Dimension() != 2) DebugStop();
        int nsides = face_el->NSides();
        TPZManVector<int,4> lineloop_debug(4,gDFN_NoIndex);
        lineloop.resize(face_el->NSides(1));
        for(int iside = face_el->NSides(0); iside<nsides-1; iside++){
            TPZGeoElSide gelside(face_el,iside);
            TPZGeoElSide neig;
            for(neig = gelside.Neighbour(); neig != gelside; neig = neig.Neighbour()){
                if(neig.Element()->Dimension()!=1) continue;
                int orientation = OrientationMatch(gelside,neig)?1:-1;
                lineloop[iside-face_el->NSides(1)] = orientation*(neig.Element()->Index()+shift);
                lineloop_debug[iside-face_el->NSides(1)] = orientation*neig.Element()->Index();
                break;
            }
        }
        return;
    }

    /// Removes negative integers from a container
    template<class Tcontainer>
    void ClearNegativeEntries(Tcontainer& list){
        Tcontainer copy(list);
        #ifdef PZDEBUG
            list.Fill(gDFN_NoIndex);
        #endif // PZDEBUG
        list.clear();
        for(auto& index : copy){
            if( !(index < 0) ) list.push_back(index);
        }
    }

    template<class Tcontainer>
    TPZGeoEl* GetLoopedFace(const Tcontainer& edges, TPZGeoMesh* gmesh){
        // Consistency checks
        int nedges = edges.size();
        if(nedges < 3 || nedges > 4) return nullptr;

        std::array<TPZGeoElSide,4> gelside;
        int i = 0;
        for(const int64_t index : edges){
            TPZGeoEl* gel = gmesh->Element(index);
            if(gel->Dimension() != 1) DebugStop();
            gelside[i] = {gel,2};
            i++;
        }

        TPZGeoEl* CandidateFace = FindCommonNeighbour(gelside[0],gelside[1],gelside[2],2);

        if(!CandidateFace) return nullptr;
        if(nedges == 3) return CandidateFace;
        // return CandidateFace;

        TPZGeoElSide neig = gelside[3].Neighbour();
        for(/*void*/; neig != gelside[3]; ++neig){
            if(neig.Element()->Index() == CandidateFace->Index())
                {return CandidateFace;}
        }

        return nullptr;
    }

} /*namespace DFN*/



#endif /* DFNNamespaceTemp_cpp */
