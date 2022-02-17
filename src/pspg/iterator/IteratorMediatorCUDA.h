#ifndef PSPG_ITERATOR_MEDIATOR_CUDA_H
#define PSPG_ITERATOR_MEDIATOR_CUDA_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/containers/DArray.h>
#include <pscf/iterator/IteratorMediator.h>
#include <pspg/System.h>

namespace Pscf {
namespace Pspg{

   using namespace Util;

   // typedef DArray<double> FieldCUDA; WHAT IS THIS TYPE WOW IDK YET

   template <int D>
   class IteratorMediatorCUDA : public IteratorMediator<FieldCUDA>
   {
   public:

      /// Constructor
      IteratorMediatorCUDA(System<D>& sys);

      /// Destructor
      ~IteratorMediatorCUDA(); 

      /// Set iterator pointer
      void setIterator(Iterator<FieldCUDA>& iter); 

      /// Set up the iterator.
      void setup();

      /// Instructor the iterator to solve the fixed point problem.
      int solve();

      /// Checks if the system has an initial guess
      bool hasInitialGuess();
      
      /// Calculates and returns the number of elements in the
      /// array to be iterated
      int nElements();

      /// Gets a reference to the current state of the system
      void getCurrent(FieldCUDA& curr);

      /// Runs calculation to evaluate function for fixed point.
      void evaluate();

      /// Gets residual values from system
      void getResidual(FieldCUDA& resid);

      /// Updates the system with a passed in state of the iterator.
      void update(FieldCUDA& newGuess);

      

   private:

      // pointer to system
      System<D>* sys_;

      // pointer to iterator
      Pscf::Iterator<FieldCUDA>* iter_;

   };

   #ifndef PSPG_ITERATOR_MEDIATOR_CUDA_TPP
   // Suppress implicit instantiation
   extern template class IteratorMediatorCUDA<1>;
   extern template class IteratorMediatorCUDA<2>;
   extern template class IteratorMediatorCUDA<3>;
   #endif

}
}
#endif