#ifndef PSPC_MC_MOVE_FACTORY_H
#define PSPC_MC_MOVE_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <pspc/simulate/mcmove/McMove.h>
#include <pspc/simulate/McSimulator.h>

#include <string>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /**
   * Factory for subclasses of McMove.
   *
   * \ingroup Pspc_McMove_Module
   */
   template <int D>
   class McMoveFactory : public Factory< McMove<D> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param mcSimulator parent McSimulator<D> object
      */
      McMoveFactory(McSimulator<D>& mcSimulator);

      /**
      * Method to create any McMove supplied with PSCF.
      *
      * \param className name of the McMove subclass
      * \return McMove* pointer to new instance of className
      */
      McMove<D>* factory(const std::string &className) const;

      using Factory< McMove<D> >::trySubfactories;

   private:

      /// Pointer to the parent mcSimulator.
      McSimulator<D>* mcSimulatorPtr_;

   };

   #ifndef PSPC_MC_MOVE_FACTORY_TPP
   // Suppress implicit instantiation
   extern template class McMoveFactory<1>;
   extern template class McMoveFactory<2>;
   extern template class McMoveFactory<3>;
   #endif

}
}
#endif
