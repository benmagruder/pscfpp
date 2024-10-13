#ifndef PRDC_CUDA_FIELD_H
#define PRDC_CUDA_FIELD_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/GpuResources.h>
#include <util/global.h>
#include <cufft.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   // Forward declaration of analogous container for data on host.
   template <typename Data> class HostField;

   /**
   * Dynamic array on the GPU device with aligned data.
   *
   * This class wraps an aligned C array with elements of type Data that is
   * allocated in GPU device global memory.  All member functions may be 
   * called from the CPU host, but the class does not offer access to 
   * individual elements via the subscript operator, operator[]
   *
   * \ingroup Prdc_Cuda_Module
   */
   template <typename Data>
   class Field
   {

   public:

      /**
      * Default constructor.
      */
      Field();

      /**
      * Copy constructor.
      * 
      * \param other Field<Data> to be copied (input)
      */
      Field(Field<Data> const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~Field();

      /**
      * Allocate the underlying C array on the device.
      *
      * \throw Exception if the Field is already allocated.
      *
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the Field is not allocated.
      */
      void deallocate();

      /**
      * Assignment operator, assignment from another Field<Data> array.
      *  
      * Performs a deep copy, by copying values of all elements from device
      * memory to device memory.
      *
      * This function will allocate memory if this (LHS) Field is not 
      * allocated.  If this is allocated, it must have the same dimensions 
      * as the RHS Field<Data>.
      *
      * \param other Field<Data> on rhs of assignent (input)
      */
      virtual Field<Data>& operator = (const Field<Data>& other);

      /**
      * Assignment operator, assignment from HostField<Data> host array.
      *
      * Performs a deep copy from a RHS HostField<Data> host array to this 
      * LHS Field<D> device array, by copying underlying C array from host
      * memory to device memory.
      *
      * This function will allocate memory if this (LHS) Field<D> is not 
      * allocated.  If this is allocated, it must have the same dimensions 
      * as the RHS HostField<D>.
      *
      * \param other HostField<Data> on RHS of assignent (input)
      */
      virtual Field<Data>& operator = (const HostField<Data>& other);

      /**
      * Return allocated capacity.
      *
      * \return Number of elements allocated in array
      */
      int capacity() const;

      /**
      * Return true if the Field has been allocated, false otherwise.
      */
      bool isAllocated() const;

      /**
      * Return pointer to underlying C array.
      */
      Data* cField();

      /**
      * Return pointer to const to underlying C array.
      */
      Data const * cField() const;

   protected:

      /// Pointer to a C array of Data elements on the GPU device.
      Data* data_;

      /// Allocated size (capacity) of the data_ array.
      int capacity_;

   };

   /*
   * Return allocated capacity.
   */
   template <typename Data>
   inline int Field<Data>::capacity() const
   {  return capacity_; }

   /*
   * Get a pointer to the underlying C array.
   */
   template <typename Data>
   inline Data* Field<Data>::cField()
   {  return data_; }

   /*
   * Get a pointer to const to the underlying C array.
   */
   template <typename Data>
   inline const Data* Field<Data>::cField() const
   {  return data_; }

   /*
   * Return true if the Field has been allocated, false otherwise.
   */
   template <typename Data>
   inline bool Field<Data>::isAllocated() const
   {  return (bool)data_; }

   #ifndef PRDC_CUDA_FIELD_TPP
   extern template class Field<cudaReal>;
   extern template class Field<cudaComplex>;
   #endif

} // namespace Cuda
} // namespace Prdc
} // namespace Pscf
#endif
