#ifndef PSPG_FFT_H
#define PSPG_FFT_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspg/field/RDField.h>
#include <pspg/field/RDFieldDft.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <cufft.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <pspg/math/GpuResources.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;
   using namespace Pscf;

   /**
   * Fourier transform wrapper for real data.
   *
   * \ingroup Pspg_Field_Module
   */
   template <int D>
   class FFT 
   {

   public:

      /**
      * Default constructor.
      */
      FFT();

      /**
      * Destructor.
      */
      virtual ~FFT();

      /**
      * Setup grid dimensions, plans and work space.
      *
      * \param meshDimensions Dimensions of real-space grid.
      */
      void setup(IntVec<D> const & meshDimensions);

      /**
      * Check and setup grid dimensions if necessary.
      *
      * \param rDField  real data on r-space grid (device mem)
      * \param kDField  complex data on k-space grid (device mem)
      */
      void setup(RDField<D>& rDField, RDFieldDft<D>& kDField);

      /**
      * Compute forward (real-to-complex) discrete Fourier transform.
      *
      * \param rField  real values on r-space grid (input, gpu mem)
      * \param kField  complex values on k-space grid (output, gpu mem)
      */
      void forwardTransform(RDField<D> & rField, RDFieldDft<D>& kField) 
      const;

      /**
      * Compute forward Fourier transform without destroying input.
      *
      * \param rField  real values on r-space grid (input, gpu mem)
      * \param kField  complex values on k-space grid (output, gpu mem)
      */
      void forwardTransformSafe(RDField<D> const & rField, 
                                RDFieldDft<D>& kField) const;

      /**
      * Compute inverse (complex-to-real) discrete Fourier transform.
      *
      * \param kField  complex values on k-space grid (input, gpu mem)
      * \param rField  real values on r-space grid (output, gpu mem)
      */
      void inverseTransform(RDFieldDft<D> & kField, RDField<D>& rField) 
      const;

      /**
      * Compute inverse (complex to real) DFT without destroying input.
      *
      * \param kField  complex values on k-space grid (input, gpu mem)
      * \param rField  real values on r-space grid (output, gpu mem)
      */
      void inverseTransformSafe(RDFieldDft<D> const & kField, 
                                RDField<D>& rField) const;

      /**
      * Return the dimensions of the grid for which this was allocated.
      */
      const IntVec<D>& meshDimensions() const;

      /**
      *  Has this FFT object been setup?  
      */
      bool isSetup() const;

      /**
      *  Get the plan for the forward DFT.
      */
      cufftHandle& fPlan();

      /**
      *  Get the plan for the inverse DFT.
      */
      cufftHandle& iPlan();

   private:

      // Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

      // Private r-space array for performing safe transforms.
      mutable RDField<D> rFieldCopy_;

      // Private k-space array for performing safe transforms.
      mutable RDFieldDft<D> kFieldCopy_;

      // Number of points in r-space grid
      int rSize_;

      // Number of points in k-space grid
      int kSize_;

      // Pointer to a plan for a forward transform.
      cufftHandle fPlan_;

      // Pointer to a plan for an inverse transform.
      cufftHandle iPlan_;

      // Have array dimension and plan been initialized?
      bool isSetup_;

      /**
      * Make FFTW plans for transform and inverse transform.
      */
      void makePlans(RDField<D>& rDField, RDFieldDft<D>& kDField);

   };

   // Declarations of explicit specializations

   template <>
   void FFT<1>::makePlans(RDField<1>& rField, RDFieldDft<1>& kField);

   template <>
   void FFT<2>::makePlans(RDField<2>& rField, RDFieldDft<2>& kField);

   template <>
   void FFT<3>::makePlans(RDField<3>& rField, RDFieldDft<3>& kField);

   /*
   * Return the dimensions of the grid for which this was allocated.
   */
   template <int D>
   inline const IntVec<D>& FFT<D>::meshDimensions() const
   {  return meshDimensions_; }

   template <int D>
   inline bool FFT<D>::isSetup() const
   { return isSetup_; }

   template <int D>
   inline cufftHandle& FFT<D>::fPlan() 
   { return fPlan_; }

   template <int D>
   inline cufftHandle& FFT<D>::iPlan()
   { return iPlan_; }


   #ifndef PSPG_FFT_TPP
   // Suppress implicit instantiation
   extern template class FFT<1>;
   extern template class FFT<2>;
   extern template class FFT<3>;
   #endif

static __global__ 
void scaleRealData(cudaReal* data, cudaReal scale, int size) {
   //write code that will scale
   int nThreads = blockDim.x * gridDim.x;
   int startId = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startId; i < size; i += nThreads ) {
      data[i] *= scale;
   }
}

}
}



//#include "FFT.tpp"
#endif
