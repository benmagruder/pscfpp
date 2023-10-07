#ifndef PSPG_WAVE_LIST_H
#define PSPG_WAVE_LIST_H
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/gpu/RFieldDft.h>
#include <prdc/gpu/RField.h>

#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <util/containers/DArray.h>
#include <util/containers/GArray.h>
#include <util/containers/DMatrix.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Container for wavevector data.
   */
   template <int D>
   class WaveList{
   public:

      /**
      * Constructor.
      */
      WaveList();

      /**
      * Destructor
      */
      ~WaveList();

      /**
      * Allocate memory for all arrays. Call in readParameters.
      *
      * \param mesh  spatial discretization mesh (input)
      * \param unitCell  crystallographic unit cell (input)
      */
      void allocate(Mesh<D> const & mesh, 
                    UnitCell<D> const & unitCell);

      /**
      * Compute minimum images of wavevectors.
      *
      * This is only done once, at beginning, using mesh and the initial 
      * unit cell. This can also be called in the readParameters function.
      *
      * \param mesh  spatial discretization Mesh<D> object
      * \param unitCell  crystallographic UnitCell<D> object
      */
      void computeMinimumImages(Mesh<D> const & mesh, 
                                UnitCell<D> const & unitCell);

      /**
      * Compute square norm |k|^2 for all wavevectors.
      *
      * Called once per iteration with unit cell relaxation.
      * Implementation can copy geometric data (rBasis, kBasis, etc.)
      * into local data structures and implement the actual 
      * calculation of kSq or dKSq for each wavevector on the GPU.
      *
      * \param unitCell crystallographic UnitCell<D>
      */
      void computeKSq(UnitCell<D> const & unitCell);

      /**
      * Compute derivatives of |k|^2 w/ respect to unit cell parameters.
      *
      * Called once per iteration with unit cell relaxation.
      *
      * \param unitCell crystallographic UnitCell<D>
      */
      void computedKSq(UnitCell<D> const & unitCell);

      /**
      * Get the minimum image vector for a specified wavevector.
      *
      * \param i index for wavevector
      */
      const IntVec<D>& minImage(int i) const;

      /**
      * Get a pointer to the kSq array on device.
      */
      cudaReal* kSq() const;

      /**
      * Get a pointer to the dkSq array on device.
      */
      cudaReal* dkSq() const;

      /**
      * Get size of k-grid (number of wavewavectors).
      */
      int kSize() const;

      bool isAllocated() const
      {  return isAllocated_; }

      bool hasMinimumImages() const
      {  return hasMinimumImages_; }

   private:

      // Bare C array holding precomputed minimum images
      int* minImage_d;

      // Bare C array holding values of kSq_
      cudaReal*  kSq_;

      // Bare C array holding values of dkSq_
      cudaReal*  dkSq_;

      cudaReal* dkkBasis_d;
      cudaReal* dkkBasis;

      int* partnerIdTable;
      int* partnerIdTable_d;

      int* selfIdTable;
      int* selfIdTable_d;

      bool* implicit;
      bool* implicit_d;

      IntVec<D> dimensions_;
      int kSize_;
      int rSize_;
      int nParams_;

      DArray< IntVec<D> > minImage_;

      bool isAllocated_;

      bool hasMinimumImages_;

   };

   template <int D>
   inline const IntVec<D>& WaveList<D>::minImage(int i) const
   {  return minImage_[i]; }

   template <int D>
   inline cudaReal* WaveList<D>::kSq() const
   {  return kSq_; }

   template <int D>
   inline cudaReal* WaveList<D>::dkSq() const
   {  return dkSq_; }

   template <int D>
   inline int WaveList<D>::kSize() const
   { return kSize_; }

   #ifndef PSPG_WAVE_LIST_TPP
   // Suppress implicit instantiation
   extern template class WaveList<1>;
   extern template class WaveList<2>;
   extern template class WaveList<3>;
   #endif

}
}
//#include "WaveList.tpp"
#endif
