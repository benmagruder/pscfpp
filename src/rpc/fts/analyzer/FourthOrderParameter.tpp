#ifndef RPC_FOURTH_ORDER_PARAMETER_TPP
#define RPC_FOURTH_ORDER_PARAMETER_TPP

#include "FourthOrderParameter.h"

#include <rpc/fts/simulator/Simulator.h>
#include <rpc/System.h>

#include <prdc/cpu/RField.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/containers/DArray.h>
#include <util/param/ParamComposite.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/global.h>

#include <fftw3.h>

#include <iostream>
#include <complex>
#include <vector>
#include <numeric>
#include <cmath>
#include <set>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   FourthOrderParameter<D>::FourthOrderParameter(Simulator<D>& simulator, System<D>& system) 
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      kSize_(1),
      hasAverage_(true),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("FourthOrderParameter"); }


   /*
   * Read parameters from file, and allocate memory.
   */
   template <int D>
   void FourthOrderParameter<D>::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      readOptional(in, "hasAverage", hasAverage_);
      readOptional(in,"nSamplePerBlock", nSamplePerBlock_);
      
      system().fileMaster().openOutputFile(outputFileName(), outputFile_);
      outputFile_ << "    chi       " << "FourthOrderParameter" << "\n";
   }
   
   /*
   * FourthOrderParameter setup
   */
   template <int D>
   void FourthOrderParameter<D>::setup() 
   {
      //Check if the system is AB diblock copolymer
      const int nMonomer = system().mixture().nMonomer();
      if (nMonomer != 2) {
         UTIL_THROW("The FourthOrderParameter Analyzer is designed specifically for diblock copolymer system. Please verify the number of monomer types in your system.");
      }
      
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();

      // Compute Fourier space dimension
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = dimensions[i];
            kSize_ *= dimensions[i];
         } else {
            kMeshDimensions_[i] = dimensions[i]/2 + 1;
            kSize_ *= (dimensions[i]/2 + 1);
         }
      }

      // Allocate variables
      if (!isInitialized_){
         wK_.allocate(dimensions);
         prefactor_.allocate(kSize_);
         for (int i = 0; i < kSize_; ++i){
            prefactor_[i] = 0;
         }
      }
         
      isInitialized_ = true;
      
      // Clear accumulators
      if (hasAverage_){
         accumulator_.clear();
      }
      
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
      
      computePrefactor();
   }

   /* 
   * Increment structure factors for all wavevectors and modes.
   */
   template <int D>
   void FourthOrderParameter<D>::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) return;
      computeFourthOrderParameter();
      
      if (hasAverage_){
         accumulator_.sample(FourthOrderParameter_);
      }
      
      double chi =  system().interaction().chi(0,1);
      UTIL_CHECK(outputFile_.is_open());
      outputFile_ << Dbl(chi);
      outputFile_ << Dbl(FourthOrderParameter_);
      outputFile_<< "\n";
   }
   
   template <int D>
   void FourthOrderParameter<D>::computeFourthOrderParameter()
   {
      UTIL_CHECK(system().w().hasData());
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      
      MeshIterator<D> itr;
      itr.setDimensions(kMeshDimensions_);
      std::vector<double> psi(kSize_);
      
      // Conver W_(r) to fourier mode W_(k)
      system().domain().fft().forwardTransform(simulator().wc(0), wK_);
      
      for (itr.begin(); !itr.atEnd(); ++itr) {
         std::complex<double> wK(wK_[itr.rank()][0], wK_[itr.rank()][1]);
         psi[itr.rank()] = std::norm(wK) * std::norm(wK);
         psi[itr.rank()] *= prefactor_[itr.rank()];
      }
      
      // Get sum over all wavevectors
      FourthOrderParameter_ = std::accumulate(psi.begin(), psi.end(), 0.0);
      FourthOrderParameter_ = std::pow(FourthOrderParameter_, 0.25);
      
      #if 0
      // Debugging output
      IntVec<D> meshDimensions = system().domain().mesh().dimensions(); 
      UnitCell<D> const & unitCell = system().domain().unitCell();
      IntVec<D> G; 
      IntVec<D> Gmin;
      IntVec<D> nGmin;
      double kSq;
      std::vector<double> k(kSize_);
      
      // Calculate GminList
      for (itr.begin(); !itr.atEnd(); ++itr){
         G = itr.position();
         Gmin = shiftToMinimum(G, meshDimensions, unitCell);
         kSq = unitCell.ksq(Gmin);
         k[itr.rank()] = kSq;
      }
      
      auto maxIt = std::max_element(psi.begin(), psi.end());
      
      // Calculate the index of the maximum element
      size_t maxIndex = std::distance(psi.begin(), maxIt);
      double kmax = k[maxIndex];
      
      Log::file() << std::endl;
      for (itr.begin(); !itr.atEnd(); ++itr){
         if (k[itr.rank()] == kmax){
            G = itr.position();
            Gmin = shiftToMinimum(G, meshDimensions, unitCell);
            Log::file() << "ksq: " << k[itr.rank()] << std::endl;
            Log::file() << " G: " << G<< std::endl;
            Log::file() << " Gmin: " << Gmin<< std::endl;
            Log::file() << " prefactor: " <<  prefactor_[itr.rank()]<< std::endl;
            Log::file() << " psi: " <<  psi[itr.rank()]<< std::endl;
         }
      
      }
      #endif
      
   }
   
   template <int D>
   void FourthOrderParameter<D>::computePrefactor()
   {
      IntVec<D> meshDimensions = system().domain().mesh().dimensions(); 
      UnitCell<D> const & unitCell = system().domain().unitCell();
      IntVec<D> G; 
      IntVec<D> Gmin;
      IntVec<D> nGmin;
      DArray<IntVec<D>> GminList;
      GminList.allocate(kSize_);
      MeshIterator<D> itr(kMeshDimensions_);
      MeshIterator<D> searchItr(kMeshDimensions_);
      
      // Calculate GminList
      for (itr.begin(); !itr.atEnd(); ++itr){
         G = itr.position();
         Gmin = shiftToMinimum(G, meshDimensions, unitCell);
         GminList[itr.rank()] = Gmin;
      }
      
      // Compute prefactor for each G wavevector
      for (itr.begin(); !itr.atEnd(); ++itr){
         bool inverseFound = false;

         // If prefactor of current wavevector has not been assigned
         if (prefactor_[itr.rank()] == 0){
            Gmin = GminList[itr.rank()];
            
            // Compute inverse of wavevector
            nGmin.negate(Gmin);
            
            // Search for inverse of wavevector
            searchItr = itr;
            for (; !searchItr.atEnd(); ++searchItr){
               if (nGmin == GminList[searchItr.rank()]){
                  prefactor_[itr.rank()] = 1.0/2.0;
                  prefactor_[searchItr.rank()] = 1.0/2.0;
                  inverseFound = true;
               }
            }
            
            if (inverseFound == false){
               prefactor_[itr.rank()]  = 1.0;
            }
            
         }
         
      }
   }
   
   /*
   * Output final results to output file.
   */
   template <int D>  
   void FourthOrderParameter<D>::output() 
   {
      if (hasAverage_){
         Log::file() << std::endl;
         Log::file() << "At chi = " << system().interaction().chi(0,1) << "\n";
         Log::file() << "Time average of the FourthOrderParameter is: "
                     << Dbl(accumulator_.average())
                     << " +- " << Dbl(accumulator_.blockingError(), 9, 2) 
                     << "\n";
      } 
      
   }

}
}
#endif 