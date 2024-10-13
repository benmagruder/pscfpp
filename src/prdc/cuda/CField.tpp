#ifndef PRDC_CUDA_C_FIELD_TPP
#define PRDC_CUDA_C_FIELD_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/CField.h>
#include <prdc/cuda/HostField.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;
   using namespace Pscf;

   /**
   * Default constructor.
   */
   template <int D>
   CField<D>::CField()
    : Field<cudaComplex>()
   {}

   /*
   * Destructor.
   */
   template <int D>
   CField<D>::~CField()
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the Field to be copied.
   */
   template <int D>
   CField<D>::CField(const CField<D>& other)
    : Field<cudaComplex>(other),
      meshDimensions_(0)
   {
      meshDimensions_ = other.meshDimensions_;
   }

   /*
   * Assignment from another RField<D>.
   */
   template <int D>
   CField<D>& CField<D>::operator = (const CField<D>& other)
   {
      Field<cudaComplex>::operator = (other);
      meshDimensions_ = other.meshDimensions_;

      return *this;
   }

   /*
   * Assignment from RHS HostField<Data> host array.
   */
   template <int D>
   CField<D>& CField<D>::operator = (const HostField<cudaComplex>& other)
   {
      // Preconditions: both arrays must be allocated with equal capacities
      if (!other.isAllocated()) {
         UTIL_THROW("Error: RHS HostField<cudaComplex> is not allocated.");
      }
      if (!isAllocated()) {
         UTIL_THROW("Error: LHS CField<D> is not allocated.");
      }
      if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign Fields of unequal capacity");
      }

      // Use base class assignment operator to copy elements
      Field<cudaComplex>::operator = (other);

      return *this;
   }

   /*
   * Allocate the underlying C array sized for an associated mesh.
   */
   template <int D>
   void CField<D>::allocate(const IntVec<D>& meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         size *= meshDimensions[i];
      }
      Field<cudaComplex>::allocate(size);
   }

}
}
}
#endif
