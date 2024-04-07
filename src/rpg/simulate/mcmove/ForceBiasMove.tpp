#ifndef RPG_FORCE_BIAS_MOVE_TPP
#define RPG_FORCE_BIAS_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ForceBiasMove.h"
#include "McMove.h" 
#include <rpg/simulate/mcmove/McSimulator.h>
#include <rpg/compressor/Compressor.h>
#include <rpg/System.h>
#include <util/param/ParamComposite.h>
#include <util/random/Random.h>
#include <pscf/cuda/CudaRandom.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ForceBiasMove<D>::ForceBiasMove(McSimulator<D>& simulator) 
    : McMove<D>(simulator),
      w_(),
      dwc_(),
      mobility_(0.0)
   {  setClassName("ForceBiasMove"); }

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   ForceBiasMove<D>::~ForceBiasMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void ForceBiasMove<D>::readParameters(std::istream &in)
   {
      //Read the probability
      readProbability(in);

      // Read Brownian dynamics mobility parameter
      read(in, "mobility", mobility_);

      // Allocate memory for private containers
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      w_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         w_[i].allocate(meshDimensions);
      }
      dc_.allocate(nMonomer-1);
      dwc_.allocate(nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         dc_[i].allocate(meshDimensions);
         dwc_[i].allocate(meshDimensions);
      }

   }
   
   template <int D>
   void ForceBiasMove<D>::setup()
   {  
      // Check array capacities
      int meshSize = system().domain().mesh().size();
      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(w_.capacity() == nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         UTIL_CHECK(w_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dc_.capacity() == nMonomer-1);
      UTIL_CHECK(dwc_.capacity() == nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         UTIL_CHECK(dc_[i].capacity() == meshSize);
         UTIL_CHECK(dwc_[i].capacity() == meshSize);
      }

      McMove<D>::setup();
      biasField_.allocate(meshSize);
      dwd_.allocate(meshSize);
      gaussianField_.allocate(meshSize);
   }
 
   /*
   * Attempt unconstrained move
   */
   template <int D>
   bool ForceBiasMove<D>::move()
   {
      totalTimer_.start();
      incrementNAttempt();

      // Preconditions
      UTIL_CHECK(simulator().hasWc());
      UTIL_CHECK(simulator().hasDc());
      UTIL_CHECK(simulator().hasHamiltonian());

      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j;
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);

      // Get current Hamiltonian
      double oldHamiltonian = simulator().hamiltonian();

      // Save current state
      simulator().saveState();

      // Clear both eigen-components of the fields and hamiltonian
      simulator().clearData();

      attemptMoveTimer_.start();
      
      // Copy current W fields from parent system into wc_
      DArray<RField<D>> const * Wr = &system().w().rgrid();
      for (i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks, nThreads>>>
            (w_[i].cField(), (*Wr)[i].cField(), meshSize);
      }

      // Copy current derivative fields from into member variable dc_
      DArray<RField<D>> const * Dc = &simulator().dc();
      for (i = 0; i < nMonomer - 1; ++i) {
         assignReal<<<nBlocks, nThreads>>>
            (dc_[i].cField(), (*Dc)[i].cField(), meshSize);
      }

      // Constants for dynamics
      const double vSystem = system().domain().unitCell().volume();
      const double vNode = vSystem/double(meshSize);
      double a = -1.0*mobility_;
      double b = sqrt(2.0*mobility_/vNode);
      
      // Constants for normal distribution
      double stddev = 1.0;
      double mean = 0;

      // Modify local variables dwc_ and wc_
      // Loop over eigenvectors of projected chi matrix
      double evec;
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> & dwc = dwc_[j];
         
         // Generagte normal distributed random floating point numbers
         cudaRandom().normal(gaussianField_.cField(), meshSize, (cudaReal)stddev, (cudaReal)mean);
         
         // dwr
         scaleReal<<<nBlocks, nThreads>>>(gaussianField_.cField(), b, meshSize);
         
         // dwd
         assignReal<<<nBlocks, nThreads>>>(dwd_.cField(), dc_[j].cField(), meshSize);
         scaleReal<<<nBlocks, nThreads>>>(dwd_.cField(), a, meshSize);
         
         // dwc
         pointWiseBinaryAdd<<<nBlocks, nThreads>>>(dwd_.cField(), gaussianField_.cField(), dwc.cField(), meshSize);
         
         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            RField<D> const & dwc = dwc_[j];
            RField<D> & wn = w_[i];
            evec = simulator().chiEvecs(j,i);
            pointWiseAddScale<<<nBlocks, nThreads>>>(wn.cField(), dwc.cField(), evec, meshSize);
         }
         
      }

      // Set modified fields in parent system
      system().setWRGrid(w_);
      simulator().clearData();

      attemptMoveTimer_.stop();

      // Call compressor
      compressorTimer_.start();
      int compress = simulator().compressor().compress();
      compressorTimer_.stop();
      
      bool isConverged = false;
      if (compress != 0){
         incrementNFail();
         simulator().restoreState();
      } else {
         isConverged = true;
         
         // Compute eigenvector components of current fields
         computeWcTimer_.start();
         simulator().computeWc();
         simulator().computeCc();
         simulator().computeDc();
         computeWcTimer_.stop();

         // Evaluate new Hamiltonian
         computeHamiltonianTimer_.start();
         simulator().computeHamiltonian();
         double newHamiltonian = simulator().hamiltonian();
         double dH = newHamiltonian - oldHamiltonian;

         // Compute force bias 
         double bias = 0.0;
         for (j = 0; j < nMonomer - 1; ++j) {
            RField<D> const & di = dc_[j];
            RField<D> const & df = simulator().dc(j);
            RField<D> const & dwc = dwc_[j];
         
            computeForceBias<<<nBlocks, nThreads>>>
               (biasField_.cField(), di.cField(), df.cField(), dwc.cField(), mobility_, meshSize);
            bias += (double) gpuSum(biasField_.cField(), meshSize);
         }
         bias *= vNode;
         computeHamiltonianTimer_.stop();

         // Accept or reject move
         decisionTimer_.start();
         bool accept = false;
         double weight = exp(bias - dH);
         accept = random().metropolis(weight);
         if (accept) {
            incrementNAccept();
            simulator().clearState();
         } else {
            simulator().restoreState();
         }
         decisionTimer_.stop();
      }
      
      totalTimer_.stop();
      return isConverged;
   }

   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void ForceBiasMove<D>::output()
   {}

   template<int D>
   void ForceBiasMove<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "Real Move times contributions:\n";
      McMove<D>::outputTimers(out);
   }

}
}
#endif
