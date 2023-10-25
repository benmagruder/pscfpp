/*
* PSCF Package
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.tpp"
#ifdef PSCF_OPENMP
#include <omp.h>
#endif

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   // Explicit class instantiations

   template class FFT<1>;
   template class FFT<2>;
   template class FFT<3>;


   // Planning functions, explicit specializations.

   template<>
   void FFT<1>::makePlans(RField<1>& rField, RFieldDft<1>& kField)
   {
      #ifdef PSCF_OPENMP
      int nThread = omp_get_max_threads();
      if (nThread > 1) {
         // std::cout << "Planning 1D FFT with " << nThread << " threads\n";
         fftw_plan_with_nthreads(nThread);
      }
      #endif
      unsigned int flags = FFTW_ESTIMATE;
      fPlan_ = fftw_plan_dft_r2c_1d(rSize_, &rField[0], &kField[0], flags);
      iPlan_ = fftw_plan_dft_c2r_1d(rSize_, &kField[0], &rField[0], flags);
   }

   template <>
   void FFT<2>::makePlans(RField<2>& rField, RFieldDft<2>& kField)
   {
      #ifdef PSCF_OPENMP
      int nThread = omp_get_max_threads();
      if (nThread > 1) {
         // std::cout << "Planning 2D FFT with " << nThread << " threads\n";
         fftw_plan_with_nthreads(nThread);
      }
      #endif
      unsigned int flags = FFTW_ESTIMATE;
      fPlan_ = fftw_plan_dft_r2c_2d(meshDimensions_[0], meshDimensions_[1],
      	                           &rField[0], &kField[0], flags);
      iPlan_ = fftw_plan_dft_c2r_2d(meshDimensions_[0], meshDimensions_[1],
                                    &kField[0], &rField[0], flags);
   }

   template <>
   void FFT<3>::makePlans(RField<3>& rField, RFieldDft<3>& kField)
   {
      #ifdef PSCF_OPENMP
      int nThread = omp_get_max_threads();
      if (nThread > 1) {
         // std::cout << "Planning 3D FFT with " << nThread << " threads\n";
         fftw_plan_with_nthreads(nThread);
      }
      #endif
      unsigned int flags = FFTW_ESTIMATE;
      fPlan_ = fftw_plan_dft_r2c_3d(meshDimensions_[0], meshDimensions_[1],
      	                           meshDimensions_[2], &rField[0], &kField[0],
      	                           flags);
      iPlan_ = fftw_plan_dft_c2r_3d(meshDimensions_[0], meshDimensions_[1],
                                    meshDimensions_[2], &kField[0], &rField[0],
                                    flags);
   }

}
}
}
