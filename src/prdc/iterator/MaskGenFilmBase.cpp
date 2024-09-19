/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaskGenFilmBase.tpp"

namespace Pscf {
namespace Prdc
{

   using namespace Util;

   // Explicit Specializations for checkLatticeVectors
   
   /*
   * Check that user-defined lattice basis vectors are compatible with 
   * the thin film constraint.
   * 
   * In 1D, there is nothing to do; the lattice basis vector is correct in 
   * all cases.
   */
   template <>
   void MaskGenFilmBase<1>::checkLatticeVectors() const 
   {} // do nothing

   
   // In 2D, we require that gamma = 90°.
   template <>
   void MaskGenFilmBase<2>::checkLatticeVectors() const 
   {
      RealVec<2> a, b;
      a = systemLatticeVector(0);
      b = systemLatticeVector(1);

      double gamma = dot(a,b);
      if (gamma > 1e-8) { // Dot product between a and b should be 0
         UTIL_THROW("ERROR: Lattice basis vectors must be orthogonal");
      }
   } 

   /*
   * In 3D, we require that there be one lattice basis vector that is 
   * orthogonal to the walls (parameter with index normalVecId), and two
   * that are parallel to the walls.
   */
   template <>
   void MaskGenFilmBase<3>::checkLatticeVectors() const 
   {
      RealVec<3> a, b, c;
      a = systemLatticeVector(0);
      b = systemLatticeVector(1);
      c = systemLatticeVector(2);

      double alpha, beta, gamma;
      gamma = dot(a,b);
      beta = dot(a,c);
      alpha = dot(b,c);

      if (normalVecId() == 0) {
         if (beta > 1e-8 || gamma > 1e-8) {
            UTIL_THROW("ERROR: beta and gamma must be 90 degrees");
         }
      } else if (normalVecId() == 1) {
         if (alpha > 1e-8 || gamma > 1e-8) {
            UTIL_THROW("ERROR: alpha and gamma must be 90 degrees");
         }
      } else { // normalVecId == 2
         if (alpha > 1e-8 || beta > 1e-8) {
            UTIL_THROW("ERROR: alpha and beta must be 90 degrees");
         }
      }
   }

   template <>
   int MaskGenFilmBase<1>::convertNormalVecIdToParamId(int normalVecIndex) 
   const
   {
      UTIL_CHECK(normalVecIndex == 0);
      return 0;
   }

   template <>
   int MaskGenFilmBase<2>::convertNormalVecIdToParamId(int normalVecIndex) 
   const
   {
      UTIL_CHECK((normalVecIndex == 0) || (normalVecIndex == 1));
      if (normalVecIndex == 0) {
         return 0;
      } else { // normalVecIndex == 1
         UnitCell<2>::LatticeSystem lattice = systemLatticeSystem();
         UTIL_CHECK(lattice != UnitCell<2>::Null);
         if ((lattice == UnitCell<2>::Rectangular) ||
             (lattice == UnitCell<2>::Oblique)) {
            return 1;
         } else {
            return 0;
         }
      }
   }

   template <>
   int MaskGenFilmBase<3>::convertNormalVecIdToParamId(int normalVecIndex) 
   const
   {
      UTIL_CHECK((normalVecIndex >= 0) && (normalVecIndex <= 2));
      if (normalVecIndex == 0) {
         return 0;
      } else if (normalVecIndex == 1) { 
         UnitCell<3>::LatticeSystem lattice = systemLatticeSystem();
         UTIL_CHECK(lattice != UnitCell<3>::Null);
         if ((lattice == UnitCell<3>::Orthorhombic) ||
             (lattice == UnitCell<3>::Monoclinic) ||
             (lattice == UnitCell<3>::Triclinic)) {
            return 1;
         } else {
            return 0;
         }
      } else { // normalVecIndex == 2
         UnitCell<3>::LatticeSystem lattice = systemLatticeSystem();
         if ((lattice == UnitCell<3>::Cubic) ||
             (lattice == UnitCell<3>::Rhombohedral)) {
            return 0;
         } else if ((lattice == UnitCell<3>::Tetragonal) ||
                    (lattice == UnitCell<3>::Hexagonal)) {
            return 1;
         } else { // lattice is Orthorhombic, Monoclinic, or Triclinic
            return 2;
         }
      }
   }

   // Class declarations
   template class MaskGenFilmBase<1>;
   template class MaskGenFilmBase<2>;
   template class MaskGenFilmBase<3>;
}
}