#ifndef PSPC_REAL_MOVE_TPP
#define PSPC_REAL_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RealMove.h"
#include "McMove.h" 
#include <util/param/ParamComposite.h>
#include <pspc/System.h>
#include <util/archives/Serializable_includes.h>
#include <util/random/Random.h>


namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   RealMove<D>::RealMove(McSimulator<D>& simulator) 
    : McMove<D>(simulator),
      isAllocated_(false)
   { setClassName("RealMove"); }

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   RealMove<D>::~RealMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void RealMove<D>::readParameters(std::istream &in)
   {
      //Read the probability
      readProbability(in);
      // attampt move range [A, -A]
      read(in, "A", A_);
   }
   
   template <int D>
   void RealMove<D>::setup()
   {  
      McMove<D>::setup();
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      if (!isAllocated_){
         wFieldTmp_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            wFieldTmp_[i].allocate(meshSize);
         }
         isAllocated_ = true;
      }
   }
   
   
   
   /*
   * Attempt unconstrained move
   */
   template <int D>
   void RealMove<D>::attemptMove()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      #if 0
      //attampt move with rigid field change randomly chosen from uniform distribution [-A, A] 
      for (int i = 0; i < nMonomer; i++){
         for (int k = 0; k < meshSize; k++){
            //Random number generator
            double r = random().uniform(-A_,A_);
            //Log::file() << "random change " << r << "\n";
            wFieldTmp_[i][k] = system().w().rgrid()[i][k] + r;
         }
      }
      #endif
      
      for (int k = 0; k < meshSize; k++){
         Util::Random::SeedType seed = 0;
         random().setSeed(seed);
         double r = random().uniform(-A_,A_);
         wFieldTmp_[0][k] = system().w().rgrid()[0][k] - r;
         wFieldTmp_[1][k] = system().w().rgrid()[1][k] + r;
      }
      
      
      system().setWRGrid(wFieldTmp_);

   }


   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void RealMove<D>::output()
   {}

}
}
#endif
