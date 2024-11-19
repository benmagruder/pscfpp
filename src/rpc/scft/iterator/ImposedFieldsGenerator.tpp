#ifndef RPC_IMPOSED_FIELDS_GENERATOR_TPP
#define RPC_IMPOSED_FIELDS_GENERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ImposedFieldsGenerator.h"
#include <rpc/scft/iterator/MaskGenFilm.h>
#include <rpc/scft/iterator/ExtGenFilm.h>

namespace Pscf {
namespace Rpc {

   // Constructor
   template <int D>
   ImposedFieldsGenerator<D>::ImposedFieldsGenerator(System<D>& sys)
    : ImposedFieldsTmpl::ImposedFieldsTmpl(),
      sysPtr_(&sys)
   {  setClassName("ImposedFieldsGenerator"); }

   // Destructor
   template <int D>
   ImposedFieldsGenerator<D>::~ImposedFieldsGenerator()
   {
      if (fieldGenPtr1_) {
         delete fieldGenPtr1_;
      }
      if (fieldGenPtr2_) {
         delete fieldGenPtr2_;
      }
   }

   // Modify the stress value if necessary.
   template <int D>
   double ImposedFieldsGenerator<D>::modifyStress(int paramId, double stress) 
   const
   {
      if (type() == "film") {
         return fieldGenPtr1_->modifyStress(paramId, stress);
      } else {
         return stress;
      }
   }

   // Create FieldGenerator objects for the mask & external field
   template <int D>
   void ImposedFieldsGenerator<D>::createGenerators()
   {
      if (type() == "film") {
         fieldGenPtr1_ = new MaskGenFilm<D>(*sysPtr_);
         fieldGenPtr2_ = new ExtGenFilm<D>(*sysPtr_);
      } else {
         UTIL_THROW(("Unrecognized type parameter: " + type()).c_str());
      }
   }
   
}
}
#endif