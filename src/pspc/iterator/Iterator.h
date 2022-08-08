#ifndef PSPC_ITERATOR_H
#define PSPC_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>
#include <util/global.h>                  

namespace Pscf {
namespace Pspc
{

   template <int D>
   class System;

   using namespace Util;

   typedef DArray<double> FieldCPU;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Pspc_Iterator_Module
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
      * Constructor.
      * 
      * \param system parent System object
      */
      Iterator(System<D>& system);

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

      /**
      * Return const reference to parent system.
      */
      System<D> const & system() const
      {  return *sysPtr_;}

      /**
      * Return true if unit cell is flexible, false if rigid.
      */
      bool const isFlexible() const 
      {  return isFlexible_;}

      /**
      * Sets the flexibleParams_ array.
      * 
      * \param flexParams array of indices of each flexible lattice parameter
      */
      void setFlexibleParams(FSArray<int, 6> const & flexParams)
      {  
         flexibleParams_ = flexParams; 
         if (flexibleParams_.size() == 0) {
            isFlexible_ = false;
         } else {
            isFlexible_ = true;
         }
      }

      /**
      * Return an array containing indices of each flexible lattice 
      * parameter, e.g. [0,1,2] if there are 3 lattice parameters and
      * all 3 are flexible.
      */
      FSArray<int, 6> const flexibleParams() const
      {  return flexibleParams_; }

   protected:

      /**
      * Return reference to parent system.
      */
      System<D>& system() 
      {  return *sysPtr_;}

      /// Is the unit cell flexible during iteration?
      bool isFlexible_;

      /// Array of indices of the lattice parameters that are flexible
      FSArray<int, 6> flexibleParams_;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;
      
   };

   template<int D>
   inline Iterator<D>::Iterator()
   {  setClassName("Iterator"); }

   template<int D>
   inline Iterator<D>::Iterator(System<D>& system)
   : sysPtr_(&system)
   {  setClassName("Iterator"); }

   template<int D>
   inline Iterator<D>::~Iterator()
   {}

} // namespace Pspc
} // namespace Pscf
#endif
