#ifndef RPC_PERTURBATION_TPP
#define RPC_PERTURBATION_TPP

#include "Perturbation.h"
#include <rpc/simulate/Simulator.h>
#include <prdc/cpu/RField.h>

#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /* 
   * Constructor.
   */
   template <int D>
   Perturbation<D>::Perturbation(Simulator<D>& simulator)
    : ParamComposite(),
      simulatorPtr_(&simulator),
      lambda_(1.0),
      savedLambda_(1.0)
   {}

   /* 
   * Destructor.
   */
   template <int D>
   Perturbation<D>::~Perturbation()
   {}

   /*
   * Read parameters from stream, empty default implementation.
   */
   template <int D>
   void Perturbation<D>::readParameters(std::istream& in)
   {}

   /*
   * Setup before simulation, empty default implementation.
   */
   template <int D>
   void Perturbation<D>::setup()
   {}

   /*
   * Compute and return perturbation, default implementation.
   */
   template <int D>
   double Perturbation<D>::hamiltonian(double unperturbedHamiltonian)
   {  return 0.0; }

   /*
   * Modify functional derivatives, empty default implementation.
   */
   template <int D>
   void Perturbation<D>::incrementDc(DArray< RField<D> > & dc)
   {}

   /*
   * Save internal variables (including lambda).
   */
   template <int D>
   void Perturbation<D>::saveState()
   {  savedLambda_ = lambda_; }

   /*
   * Restored saved internal variables (including lambda).
   */
   template <int D>
   void Perturbation<D>::restoreState()
   {  lambda_ = savedLambda_; }

}
}
#endif 
