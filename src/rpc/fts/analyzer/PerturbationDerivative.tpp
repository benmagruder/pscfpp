#ifndef RPC_PERTURBATION_DERIVATIVE_TPP
#define RPC_PERTURBATION_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PerturbationDerivative.h"

#include <rpc/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/perturbation/Perturbation.h>

namespace Pscf {
namespace Rpc 
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   PerturbationDerivative<D>::PerturbationDerivative(Simulator<D>& simulator, 
                                                     System<D>& system) 
    : ThermoDerivativeAnalyzer<D>(simulator, system)
   { setClassName("PerturbationDerivative"); }

   /*
   * Destructor.
   */
   template <int D>
   PerturbationDerivative<D>::~PerturbationDerivative() 
   {}

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void PerturbationDerivative<D>::readParameters(std::istream& in) 
   {
      ThermoDerivativeAnalyzer<D>::readParameters(in);
   }
   
   /*
   * Setup before simulation loop.
   */ 
   template <int D>
   void PerturbationDerivative<D>::setup()
   {}
   
   template <int D>
   double PerturbationDerivative<D>::computeDerivative()
   {
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(simulator().hasPerturbation());
      if (!system().hasCFields()) {
         system().compute();
      }
      
      // Obteain Hamiltonian per monomer
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }
      
      return simulator().perturbation().df(); 
   }
   
   template <int D>
   double PerturbationDerivative<D>::variable()
   { return simulator().perturbation().lambda(); }
   
   template <int D>
   std::string PerturbationDerivative<D>::parameterType()
   { return "Perturbation Derivative"; }

}
}
#endif