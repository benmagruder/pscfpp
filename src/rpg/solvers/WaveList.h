#ifndef RPG_WAVE_LIST_H
#define RPG_WAVE_LIST_H
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/RField.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>
#include <pscf/cuda/DeviceArray.h> 

#include <util/containers/DArray.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Class to calculate and store properties of wavevectors.
   * 
   * In particular, minimum images, square norms of wavevectors (kSq), and 
   * derivatives of the square norms of wavevectors with respect to the 
   * lattice parameters (dKSq) are calculated and stored by this class. 
   * Any time the lattice parameters change the clearUnitCellData() method 
   * should be called, which will effectively reset the WaveList object so
   * that the wavevector properties will need to be recalculated before 
   * being used.
   */
   template <int D>
   class WaveList
   {
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
      * Allocate memory and set association with a Mesh and UnitCell object. 
      *
      * \param m  spatial discretization mesh (input)
      * \param c  crystallographic unit cell (input)
      */
      void allocate(Mesh<D> const & m, UnitCell<D> const & c);

      /**
      * Clear all internal data that depends on lattice parameters.
      * 
      * Sets hasKSq_ and hasdKSq_ to false, and sets hasMinimumImages_ to
      * false if hasVariableAngle == true.
      */
      void clearUnitCellData();

      /**
      * Compute minimum images of wavevectors. (Also calculates kSq.)
      *
      * The minimum images may change if a lattice angle in the unit cell 
      * is changed, so this method should be called whenever such changes
      * occur. The method hasVariableAngle() identifies whether the
      * minimum images may change under changes in the lattice parameters.
      * 
      * In the process of computing the minimum images, the square norm
      * |k|^2 for all wavevectors is also calculated and stored, so it is
      * not necessary to call computeKSq after calling this method. 
      * computeKSq is provided to allow calculation of kSq without 
      * recalculating minimum images.
      */
      void computeMinimumImages();

      /**
      * Compute sq. norm |k|^2 for all wavevectors, using existing min images.
      */
      void computeKSq();

      /**
      * Compute derivatives of |k|^2 w/ respect to unit cell parameters.
      */
      void computedKSq();

      /**
      * Get the array of minimum images on the device by reference.
      * 
      * The array has size kSize * D, where kSize is the number of grid points 
      * in reciprocal space. The array is unwrapped into a linear array in an
      * index-by-index manner, in which the first kSize elements of the array 
      * contain the first index of each minimum image, and so on. 
      */
      DeviceArray<int> const & minImages() const;

      /**
      * Get the kSq array on the device by reference.
      */
      RField<D> const & kSq() const;

      /**
      * Get the full dKSq array on the device by reference.
      * 
      * The array has size kSize * nParams, where kSize is the number of grid 
      * points in reciprocal space and nParams is the number of lattice 
      * parameters. The array is unwrapped into a linear array in which the 
      * first kSize elements of the array contain dKSq for the first lattice
      * parameter, and so on. 
      */
      DeviceArray<cudaReal> const & dKSq() const;

      /**
      * Get the slice of the dKSq array for parameter i by reference.
      * 
      * \param i index of lattice parameter
      */
      RField<D> const & dKSq(int i) const;

      /**
      * Get the implicitInverse array by reference.
      * 
      * This array is defined on a k-grid mesh, with a boolean value for 
      * each gridpoint. The boolean represents whether the inverse of the
      * wave at the given gridpoint is an implicit wave. Implicit here is
      * used to mean any wave that is outside the bounds of the k-grid.
      */
      DeviceArray<bool> const & implicitInverse() const;

      /**
      * Does this unit cell have an angle that can change?
      * 
      * The minimum images can only change if one of the lattice parameters
      * is an angle that may vary. Therefore, this method checks the crystal
      * system and returns true if there are any angles that may vary.
      */ 
      bool hasVariableAngle() const;

      /**
      * Has memory been allocated for arrays?
      */ 
      bool isAllocated() const
      {  return isAllocated_; }

      /**
      * Have minimum images been computed?
      */ 
      bool hasMinimumImages() const
      {  return hasMinimumImages_; }

      /**
      * Has the kSq array been computed?
      */ 
      bool hasKSq() const
      {  return hasKSq_; }

      /**
      * Has the dKSq array been computed?
      */ 
      bool hasdKSq() const
      {  return hasdKSq_; }

   private:

      /**
      * Array containing minimum images for each wave, stored on device.
      * 
      * The array has size kSize * D, where kSize is the number of grid points 
      * in reciprocal space. The array is unwrapped into a linear array in an
      * index-by-index manner, in which the first kSize elements of the array 
      * contain the first index of each minimum image, and so on. 
      */ 
      DeviceArray<int> minImages_;

      /// Array containing values of kSq_, stored on the device.
      RField<D> kSq_;

      /// Array containing all values of dKSq_, stored on the device.
      DeviceArray<cudaReal> dKSq_;

      /// Array of RFields, where each RField is a slice of the dKSq_ array.
      DArray<RField<D> > dKSqSlices_;

      /// Array indicating whether a given gridpoint has an implicit partner
      DeviceArray<bool> implicitInverse_;

      /// Dimensions of the mesh in reciprocal space.
      IntVec<D> kMeshDimensions_;

      /// Number of grid points in reciprocal space.
      int kSize_;

      /// Has memory been allocated for arrays?
      bool isAllocated_;

      /// Have minimum images been computed?
      bool hasMinimumImages_;

      /// Has the kSq array been computed?
      bool hasKSq_;

      /// Has the dKSq array been computed?
      bool hasdKSq_;

      /// Pointer to associated UnitCell<D> object
      UnitCell<D> const * unitCellPtr_;

      /// Pointer to associated Mesh<D> object
      Mesh<D> const * meshPtr_;

      /// Access associated UnitCell<D> by reference.
      UnitCell<D> const & unitCell() const
      {  return *unitCellPtr_; }

      /// Access associated Mesh<D> by reference.
      Mesh<D> const & mesh() const
      {  return *meshPtr_; }

   };

   // Get the array of minimum images on the device by reference.
   template <int D>
   inline DeviceArray<int> const & WaveList<D>::minImages() const
   {  
      UTIL_CHECK(hasMinimumImages_);
      return minImages_; 
   }
   
   // Get the kSq array on the device by reference.
   template <int D>
   inline RField<D> const & WaveList<D>::kSq() const
   {  
      UTIL_CHECK(hasKSq_);
      return kSq_; 
   }

   // Get the full dKSq array on the device by reference.
   template <int D>
   inline DeviceArray<cudaReal> const & WaveList<D>::dKSq() const
   {  
      UTIL_CHECK(hasdKSq_);
      return dKSq_; 
   }

   // Get a slice of the dKSq array on the device by reference.
   template <int D>
   inline RField<D> const & WaveList<D>::dKSq(int i) const
   {  
      UTIL_CHECK(hasdKSq_);
      return dKSqSlices_[i];
   }

   // Get the implicitInverse array by reference.
   template <int D>
   inline DeviceArray<bool> const & WaveList<D>::implicitInverse() const
   {  
      UTIL_CHECK(isAllocated_);
      return implicitInverse_;
   }

   #ifndef RPG_WAVE_LIST_TPP
   // Suppress implicit instantiation
   extern template class WaveList<1>;
   extern template class WaveList<2>;
   extern template class WaveList<3>;
   #endif

}
}
//#include "WaveList.tpp"
#endif
