#ifndef RPC_RAMP_FACTORY_TPP
#define RPC_RAMP_FACTORY_TPP

#include "RampFactory.h"  

// Subclasses of Ramp 
#include "LinearRamp.h"

#include <rpc/fts/simulator/Simulator.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   RampFactory<D>::RampFactory(Simulator<D>& simulator)
    : simulatorPtr_(&simulator)
   {}

   /* 
   * Return a pointer to a instance of Ramp subclass className.
   */
   template <int D>
   Ramp<D>* RampFactory<D>::factory(const std::string & className) const
   {
      Ramp<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
       
      // Try to match classname
      if (className == "LinearRamp") {
         ptr = new LinearRamp<D>(*simulatorPtr_);
      }

      return ptr;
   }

}
}
#endif
