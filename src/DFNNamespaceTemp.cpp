/*! 
 *  @authors   Pedro Lima
 *  @date      2020.10
 *  @brief     Namespace DFN template definitions
 */

#ifndef DFNNamespaceTemp_cpp
#define DFNNamespaceTemp_cpp


#include "pzgmesh.h"



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

    /** @brief Takes a simple oriented lineloop with 3 or 4 edges and create a geometric element
     * @param lineloop an oriented loop of 3 or 4 edges
    */  
    template<class Tcontainer>
    static TPZGeoEl* MeshSimplePolygon(TPZGeoMesh* gmesh, Tcontainer lineloop, int matid){
        int nnodes = lineloop.size();
        if(nnodes < 3) DebugStop();
        TPZManVector<int64_t,4> cornerindices(nnodes,-1);
        int i=0;
        for(int64_t edge : lineloop){
            int64_t index = std::abs(edge);
            if(edge < 0){
                cornerindices[i] = gmesh->Element(index)->NodeIndex(1);
            }else{
                cornerindices[i] = gmesh->Element(index)->NodeIndex(0);
            }
            i++;
        }
        int64_t index = -1;
        MElementType eltype;
        switch (nnodes){
            case  2: eltype = MElementType::EInterfaceLinear; break;
            case  3: eltype = MElementType::ETriangle; break;
            case  4: eltype = MElementType::EQuadrilateral; break;
            default: DebugStop();
        }
        TPZGeoEl* new_el = gmesh->CreateGeoElement(eltype,cornerindices,matid,index);
        return new_el;
    }

    /** @brief Check if a set of mesh nodes are coplanar
    */
    template<class TContainer>
    bool AreCoPlanar(TPZGeoMesh* gmesh, TContainer nodeindices, REAL tolerance){
        int nnodes = nodeindices.size();
        if(nnodes < 4) return true;

        TPZManVector<TPZManVector<REAL,3>,3> p(3,{0.,0.,0.});
        int j=0;
        for(int64_t inode : nodeindices){
            gmesh->NodeVec()[inode].GetCoordinates(p[j]);
            j++;
            if(j>2) break;
        }
        
        TPZManVector<REAL,3> vec0 = p[0] - p[1];
        TPZManVector<REAL,3> vec2 = p[2] - p[1];
        TPZManVector<REAL,3> vecj = {0., 0., 0.};

        TPZManVector<REAL,3> normal = CrossProduct<REAL>(vec0,vec2);
        REAL norm = DFN::Norm<REAL>(normal);
        if(fabs(norm) < tolerance){
            PZError << "\n" << __PRETTY_FUNCTION__;
            PZError << "\n\t Failed due to co-linear points";
            DebugStop();
        }
        normal[0] /= norm;
        normal[1] /= norm;
        normal[2] /= norm;

        TPZManVector<REAL,3> pj = {0., 0., 0.};
        for(int64_t inode : nodeindices){
            if(j<3) {j++;continue;}
            gmesh->NodeVec()[inode].GetCoordinates(pj);
            vecj = pj - p[1];
            REAL orth_dist = DotProduct<REAL>(normal,vecj);
            if(fabs(orth_dist) > tolerance) return false; 
        }
        return true;
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
    std::set<Ttype> set_intersection(std::set<Ttype>& set1, std::set<Ttype>& set2){
        std::set<Ttype> intersection;
        std::set<Ttype>& smaller_set =  set1.size() > set2.size() ? set2 : set1;
        std::set<Ttype>& bigger_set = !(set1.size() > set2.size()) ? set2 : set1;

        if(smaller_set.size() < 1) return intersection; //empty set

        auto end = bigger_set.end();
        for(auto& iel : smaller_set){
            auto itr = bigger_set.find(iel);
            if(itr == end) continue;
            intersection.insert(*itr);
        }
        return intersection;
    }
} /*namespace DFN*/



#endif /* DFNNamespaceTemp_cpp */
