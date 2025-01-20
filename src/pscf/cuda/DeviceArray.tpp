#ifndef PSCF_DEVICE_ARRAY_TPP
#define PSCF_DEVICE_ARRAY_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DeviceArray.h"
#include "HostDArray.h"
#include <cuda_runtime.h>

namespace Pscf {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray()
    : dataPtr_(nullptr),
      capacity_(0),
      ownerPtr_(nullptr),
      ownerCapacity_(0),
      ownerDataPtr_(nullptr)
   {}

   /*
   * Allocating constructor.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray(int capacity)
    : dataPtr_(nullptr),
      capacity_(0),
      ownerPtr_(nullptr),
      ownerCapacity_(0),
      ownerDataPtr_(nullptr)
   {  allocate(capacity); }

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray(const DeviceArray<Data>& other)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other array must be allocated.");
      }

      allocate(other.capacity_);
      cudaMemcpy(dataPtr_, other.cArray(), 
                 capacity_ * sizeof(Data), cudaMemcpyDeviceToDevice);
   }

   /*
   * Destructor.
   */
   template <typename Data>
   DeviceArray<Data>::~DeviceArray()
   {
      if (isAllocated() && isOwner()) {
         cudaFree(dataPtr_);
         capacity_ = 0;
      }
   }

   /*
   * Allocate the underlying C array.
   */
   template <typename Data>
   void DeviceArray<Data>::allocate(int capacity)
   {
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate an array");
      }
      if (ownerPtr_) {
         UTIL_THROW("Attempt to allocate array already associated with data");
      }
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      gpuErrChk(cudaMalloc((void**) &dataPtr_, capacity * sizeof(Data)));
      capacity_ = capacity;

      ownerPtr_ = nullptr;
   }

   /*
   * Associate this object with a slice of a different DeviceArray.
   */
   template <typename Data>
   void DeviceArray<Data>::associate(DeviceArray<Data>& arr, int beginId, 
                                     int capacity)
   {
      UTIL_CHECK(arr.isAllocated());
      UTIL_CHECK(arr.isOwner());
      UTIL_CHECK(beginId >= 0);
      UTIL_CHECK(capacity > 0);
      UTIL_CHECK(beginId + capacity <= arr.capacity());
      
      if (isAllocated()) {
         UTIL_THROW("Attempt to associate an already-allocated array.");
      }
      
      dataPtr_ = arr.cArray() + beginId;
      capacity_ = capacity;
      ownerPtr_ = &arr;
      ownerCapacity_ = arr.capacity();
      ownerDataPtr_ = arr.cArray();
   }

   /*
   * Deallocate the underlying C array.
   *
   * Throw an Exception if this array is not allocated.
   */
   template <typename Data>
   void DeviceArray<Data>::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Array is not allocated");
      }
      if (!isOwner()) {
         UTIL_THROW("Cannot deallocate, data not owned by this object.");
      }
      cudaFree(dataPtr_);
      dataPtr_ = nullptr; // reset to null
      capacity_ = 0;
   }

   /*
   * Dissociate this object from the associated array.
   *
   * Throw Exception if this object is not associated with another array.
   */
   template <typename Data>
   void DeviceArray<Data>::dissociate()
   {
      if (isOwner()) {
         UTIL_THROW("Cannot dissociate: this object owns the array.");
      }
      if (!isAllocated()) {
         UTIL_THROW("Cannot dissociate: no associated data found.");
      }
      dataPtr_ = nullptr; // reset to null
      capacity_ = 0;
      ownerPtr_ = nullptr; // reset to null
      ownerCapacity_ = 0;
      ownerDataPtr_ = nullptr;
   }

   /*
   * Assignment from another DeviceArray<Data>.
   */
   template <typename Data>
   DeviceArray<Data>& 
   DeviceArray<Data>::operator = (const DeviceArray<Data>& other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition - RHS DeviceArray<Data> must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("Other DeviceArray<Data> must be allocated.");
      }

      // If this is not allocated, then allocate
      if (!isAllocated()) {
         allocate(other.capacity());
      } 

      // Require equal capacity values 
      if (capacity_ != other.capacity_) {
         UTIL_THROW("Cannot assign arrays of unequal capacity");
      }

      // Copy elements
      cudaMemcpy(dataPtr_, other.cArray(), 
                 capacity_ * sizeof(Data), cudaMemcpyDeviceToDevice);

      ownerPtr_ = nullptr;

      return *this;
   }

   /*
   * Assignment of LHS DeviceArray<Data> from RHS HostDArray<Data>.
   */
   template <typename Data>
   DeviceArray<Data>& 
   DeviceArray<Data>::operator = (const HostDArray<Data>& other)
   {
      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("RHS HostDArray<Data> must be allocated.");
      }

      // Allocate this if necessary
      if (!isAllocated()) {
         allocate(other.capacity());
      } 

      // Require equal capacity values
      if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign arrays of unequal capacity");
      }

      // Copy elements
      cudaMemcpy(dataPtr_, other.cArray(), 
                 capacity_ * sizeof(Data), cudaMemcpyHostToDevice);

      ownerPtr_ = nullptr;

      return *this;
   }

   /*
   * Get a pointer to the underlying C array.
   */
   template <typename Data>
   Data* DeviceArray<Data>::cArray()
   {  
      // Make sure that owner array has not been deallocated / reallocated
      if (ownerPtr_) {
         UTIL_CHECK(ownerPtr_->isAllocated());
         UTIL_CHECK(ownerPtr_->capacity() == ownerCapacity_);
         UTIL_CHECK(ownerPtr_->cArray() == ownerDataPtr_);
      }

      return dataPtr_; 
   }

   /*
   * Get a pointer to const to the underlying C array.
   */
   template <typename Data>
   const Data* DeviceArray<Data>::cArray() const
   {  
      // Make sure that owner array has not been deallocated / reallocated
      if (ownerPtr_) {
         UTIL_CHECK(ownerPtr_->isAllocated());
         UTIL_CHECK(ownerPtr_->capacity() == ownerCapacity_);
         UTIL_CHECK(ownerPtr_->cArray() == ownerDataPtr_);
      }

      return dataPtr_; 
   }

   /*
   * Return true if the array has been allocated, false otherwise.
   */
   template <typename Data>
   bool DeviceArray<Data>::isAllocated() const
   {  
      // Make sure that owner array has not been deallocated / reallocated
      if (ownerPtr_) {
         UTIL_CHECK(ownerPtr_->isAllocated());
         UTIL_CHECK(ownerPtr_->capacity() == ownerCapacity_);
         UTIL_CHECK(ownerPtr_->cArray() == ownerDataPtr_);
      }

      return (bool)dataPtr_; 
   }

} // namespace Pscf
#endif
