#ifndef PRDC_CUDA_R_FIELD_DFT_H
#define PRDC_CUDA_R_FIELD_DFT_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Field.h"
#include <pscf/cuda/GpuResources.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>
#include <cufft.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;
   using namespace Pscf;

   /**
   * Discrete Fourier Transform (DFT) of a real field on an FFT mesh.
   *
   * The DFT is stored internally as a C array of cudaComplex elements 
   * located in global GPU memory. All member functions are C++ functions 
   * that can be called from the host CPU. 
   *
   * \ingroup Prdc_Cuda_Module
   */
   template <int D>
   class RFieldDft : public Field<cudaComplex>
   {

   public:

      /**
      * Default constructor.
      */
      RFieldDft();

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the RFieldDft to be copied.
      */
      RFieldDft(const RFieldDft<D>& other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~RFieldDft();

      /**
      * Assignment operator.
      *
      * If this Field is not allocated, allocates and copies all elements.
      *
      * If this and the other Field are both allocated, the capacities must
      * be exactly equal. If so, this method copies all elements.
      *
      * \param other the RHS Field
      */
      RFieldDft<D>& operator = (const RFieldDft<D>& other);

      /**
      * Allocate the underlying C array for an FFT grid.
      *
      * \throw Exception if the RFieldDft is already allocated.
      *
      * \param meshDimensions vector of mesh dimensions
      */
      void allocate(const IntVec<D>& meshDimensions);

      /**
      * Return vector of mesh dimensions by constant reference.
      */
      const IntVec<D>& meshDimensions() const;

      /**
      * Return vector of dft (Fourier) grid dimensions by const reference.
      *  
      * The last element of dftDimensions() and meshDimensions() differ by
      * about a factor of two: dftDimension()[D-1] = meshDimensions()/2 + 1.
      * For D > 1, other elements are equal. 
      */
      const IntVec<D>& dftDimensions() const;

      /**
      * Serialize a Field to/from an Archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      using Field<cudaComplex>::allocate;
      using Field<cudaComplex>::operator =;

   private:

      // Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

      // Vector containing dimensions of dft (Fourier) grid.
      IntVec<D> dftDimensions_;

   };

   /*
   * Allocate the underlying C array for an FFT grid.
   */
   template <int D>
   void RFieldDft<D>::allocate(const IntVec<D>& meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         if (i < D - 1) {
            dftDimensions_[i] = meshDimensions[i];
            size *= meshDimensions[i];
         } else {
            dftDimensions_[i] = (meshDimensions[i]/2 + 1); 
            size *= (meshDimensions[i]/2 + 1);
         }
      }
      Field<cudaComplex>::allocate(size);
   }

   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D>
   inline const IntVec<D>& RFieldDft<D>::meshDimensions() const
   {  return meshDimensions_; }

   /*  
   * Return dimensions of dft grid by constant reference. 
   */
   template <int D>
   inline const IntVec<D>& RFieldDft<D>::dftDimensions() const
   {  return dftDimensions_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void RFieldDft<D>::serialize(Archive& ar, const unsigned int version)
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
         cudaComplex* tempData = new cudaComplex[capacity];
         cudaMemcpy(tempData, data_, capacity * sizeof(cudaComplex), 
                    cudaMemcpyDeviceToHost);
         for (int i = 0; i < capacity_; ++i) {
            ar & tempData[i].x;
            ar & tempData[i].y;
         }
         delete[] tempData;
      }
      ar & meshDimensions_;
   }

   #ifndef PRDC_CUDA_R_FIELD_DFT_TPP
   extern template class RFieldDft<1>;
   extern template class RFieldDft<2>;
   extern template class RFieldDft<3>;
   #endif

}
}
}
//#include "RFieldDft.tpp"
#endif
