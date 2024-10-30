#ifndef PRDC_CUDA_R_FIELD_H
#define PRDC_CUDA_R_FIELD_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/Field.h>
#include <pscf/cuda/GpuResources.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   template <typename Data> class HostField;

   using namespace Util;

   /**
   * Field of real values on a regular mesh, allocated on a GPU device.
   *
   * Type cudaReal is float or double, depending on preprocessor macro.
   *
   * \ingroup Prdc_Cuda_Module 
   */
   template <int D>
   class RField : public Field<cudaReal>
   {

   public:

      /**
      * Default constructor.
      */
      RField();

      /**
      * Allocating constructor.
      *
      * Allocates memory by calling allocate(meshDimensions) internally.
      *  
      * \param meshDimensions dimensions of the associated Mesh<D>
      */
      RField(IntVec<D> const & meshDimensions);

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the RField to be copied.
      */
      RField(RField<D> const& other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~RField();

      /**
      * Allocate the underlying C array for an FFT grid.
      *
      * \throw Exception if the RField is already allocated.
      *
      * \param meshDimensions number of grid points in each direction
      */
      void allocate(IntVec<D> const & meshDimensions);

      /**
      * Assignment operator, assignment from another RField<D>.
      *
      * Performs a deep copy, by copying all elements of the RHS RField<D>
      * from device memory to device memory, and copying the meshDimensions.
      *
      * The RHS RField<D> must be allocated on entry. If this LHS object is 
      * not allocated, allocate with the required capacity.  If the LHS and
      * RHS arrays are both allocated, capacity values must be equal.
      * 
      * \param other the RHS RField
      */
      RField<D>& operator = (const RField<D>& other);

      /**
      * Assignment operator, assignment from a HostField<cudaReal>.
      *
      * Performs a deep copy, by copying all elements of the RHS RField<D>
      * from host memory to device memory.
      *
      * The RHS HostField<cudaReal> and LHS RField<D> must both be allocated
      * with equal capacity values on entry. 
      * 
      * \param other the RHS HostField<cudaReal>
      */
      RField<D>& operator = (const HostField<cudaReal>& other);

      /**
      * Return mesh dimensions by constant reference.
      */
      IntVec<D> const & meshDimensions() const;

      /**
      * Serialize a Field to/from an Archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      using Field<cudaReal>::allocate;
      using Field<cudaReal>::operator=;

   private:

      // Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

   };

   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D>
   inline const IntVec<D>& RField<D>::meshDimensions() const
   {  return meshDimensions_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void RField<D>::serialize(Archive& ar, const unsigned int version)
   {
      int capacity;
      if (Archive::is_saving()) {
         capacity = capacity_;
      }
      ar & capacity;
      if (Archive::is_loading()) {
         if (!isAllocated()) {
            if (capacity > 0) {
               allocate(capacity);
            }
         } else {
            if (capacity != capacity_) {
               UTIL_THROW("Inconsistent Field capacities");
            }
         }
      }

      if (isAllocated()) {
         double* tempData = new double[capacity];
         cudaMemcpy(tempData, data_, capacity * sizeof(cudaReal), 
                    cudaMemcpyDeviceToHost);
         for (int i = 0; i < capacity_; ++i) {
            ar & tempData[i];
         }
         delete[] tempData;
      }
      ar & meshDimensions_;
   }

   #ifndef PRDC_CUDA_R_FIELD_TPP
   extern template class RField<1>;
   extern template class RField<2>;
   extern template class RField<3>;
   #endif

}
}
}
#endif
