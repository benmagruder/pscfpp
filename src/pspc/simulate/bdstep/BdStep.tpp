#ifndef PSPC_BD_STEP_TPP
#define PSPC_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdStep.h"

#include <pspc/simulate/BdSimulator.h>
#include <pspc/System.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   BdStep<D>::BdStep(BdSimulator<D>& bdSimulator)
    : bdSimulatorPtr_(&bdSimulator),
      systemPtr_(&(bdSimulator.system())),
      randomPtr_(&(bdSimulator.random()))
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   BdStep<D>::~BdStep()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void BdStep<D>::readParameters(std::istream &in)
   {}

   /*
   * Setup at beginning of loop.
   */
   template <int D>
   void BdStep<D>::setup()
   {}

   #if 0
   template <int D>
   void BdStep<D>::step()
   {}
   #endif

   template <int D>
   void BdStep<D>::output()
   {}

   template<int D>
   void BdStep<D>::outputTimers(std::ostream& out)
   {}

   template<int D>
   void BdStep<D>::clearTimers()
   {}

}
}
#endif
