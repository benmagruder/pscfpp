#ifndef PRDC_CUDA_HOST_FIELD_TPP
#define PRDC_CUDA_HOST_FIELD_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HostField.h"
#include <util/misc/Memory.h>

#include <fftw3.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <typename Data>
   HostField<Data>::HostField()
    : data_(0),
      capacity_(0)
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the HostField to be copied.
   */
   template <typename Data>
   HostField<Data>::HostField(const HostField<Data>& other)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other HostField must be allocated.");
      }

      allocate(other.capacity_);
      cudaMemcpy(data_, other.cArray(), 
                 capacity_ * sizeof(Data), cudaMemcpyHostToHost);

   }

   /*
   * Destructor.
   */
   template <typename Data>
   HostField<Data>::~HostField()
   {
      if (isAllocated()) {
         cudaFreeHost(data_);
         capacity_ = 0;
      }
   }

   /*
   * Allocate the underlying C array.
   *
   * Throw an Exception if the HostField has already allocated.
   *
   * \param capacity number of elements to allocate.
   */
   template <typename Data>
   void HostField<Data>::allocate(int capacity)
   {
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate a HostField");
      }
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      gpuErrchk(cudaMallocHost((void**) &data_, capacity * sizeof(Data)));
      capacity_ = capacity;
   }

   /*
   * Deallocate the underlying C array.
   *
   * Throw an Exception if this HostField is not allocated.
   */
   template <typename Data>
   void HostField<Data>::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Attempt to deallocate unallocated HostField");
      }
      cudaFreeHost(data_);
      capacity_ = 0;
   }

   /*
   * Assignment from another HostField<Data> host array.
   */
   template <typename Data>
   HostField<Data>& 
   HostField<Data>::operator = (const HostField<Data>& other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other HostField must be allocated.");
      }

      // Allocate this if necessary 
      if (!isAllocated()) {
         allocate(other.capacity());
      } 

      // Require equal capacities
      if (capacity_ != other.capacity_) {
         UTIL_THROW("Cannot assign HostFields of unequal capacity");
      }

      // Copy elements
      cudaMemcpy(data_, other.cArray(), 
                 capacity_ * sizeof(Data), cudaMemcpyHostToHost);

      return *this;
   }

   /*
   * Assignment from Field<Data> RHS device array.
   */
   template <typename Data>
   HostField<Data>& 
   HostField<Data>::operator = (const Field<Data>& other)
   {
      // Precondition - RHS must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("RHS Field<Data> must be allocated.");
      }

      // Allocate this if necessary 
      if (!isAllocated()) {
         allocate(other.capacity());
      } 

      // Require equal capacities
      if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign fields of unequal capacity");
      }

      // Copy all elements
      cudaMemcpy(data_, other.cArray(), 
                 capacity_ * sizeof(Data), cudaMemcpyDeviceToHost);

      return *this;
   }

}
}
}
#endif