#ifndef PSPC_ITERATOR_FACTORY_H
#define PSPC_ITERATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/iterator/Iterator.h>
#include <util/param/Factory.h>  

#include <string>

namespace Pscf {
namespace Pspc {

   template <int D> class System;

   using namespace Util;

   /**
   * Factory for subclasses of Iterator.
   *
   * \ingroup Pspc_Iterator_Module
   */

   template <int D>
   class IteratorFactory : public Factory< Iterator<D> > 
   {

   public:

      /// Constructor
      IteratorFactory(System<D>& system);

      /**
      * Method to create any Iterator supplied with PSCF.
      *
      * \param className name of the Iterator subclass
      * \return Iterator* pointer to new instance of className
      */
      Iterator<D>* factory(const std::string &className) const;

      using Factory< Iterator<D> >::trySubfactories;

   private:

      /// Pointer to the parent system.
      System<D>* sysPtr_;

   };

   #ifndef PSPC_ITERATOR_FACTORY_TPP
   // Suppress implicit instantiation
   extern template class IteratorFactory<1>;
   extern template class IteratorFactory<2>;
   extern template class IteratorFactory<3>;
   #endif

}
}
#endif
