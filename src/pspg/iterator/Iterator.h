#ifndef PSPG_ITERATOR_H
#define PSPG_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <pspg/field/DField.h>
#include <util/global.h>                  

namespace Pscf {
namespace Pspg
{

   template <int D>
   class System;

   typedef DField<cudaReal> FieldCUDA;


   using namespace Util;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Pspg_Iterator_Module
   */
   template <int D>
   class Iterator : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      Iterator();

      /**
      * Destructor.
      */
      ~Iterator();

      /**
      * Setup iterator.
      */
      virtual void setup() = 0;

      /**
      * Iterate to solution.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int solve() = 0;
      
   };

   template<int D>
   inline Iterator<D>::Iterator()
   {  setClassName("Iterator"); }

   template<int D>
   inline Iterator<D>::~Iterator()
   {}

} // namespace Pspg
} // namespace Pscf
#endif
