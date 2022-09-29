#ifndef PSPC_FIELD_CONTAINER_TPP
#define PSPC_FIELD_CONTAINER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldContainer.h"
#include <pspc/field/FieldIo.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   FieldContainer<D>::FieldContainer()
    : basis_(),
      rgrid_(),
      fieldIoPtr_(0),
      meshDimensions_(),
      meshSize_(0),
      nBasis_(0),
      isAllocated_(false),
      hasData_(false),
      isSymmetric_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   FieldContainer<D>::~FieldContainer()
   {}

   /*
   * Create an association with a FieldIo object.
   */
   template <int D>
   void FieldContainer<D>::setFieldIo(FieldIo<D> const & fieldIo)
   {  fieldIoPtr_ = &fieldIo; }

   /*
   * Allocate memory for fields.
   */
   template <int D>
   void FieldContainer<D>::allocate(int nMonomer, int nBasis, 
                                    IntVec<D> const & meshDimensions)
   {
      // Set mesh and basis dimensions
      nMonomer_ = nMonomer;
      meshDimensions_ = meshDimensions;
      meshSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshSize_ *= meshDimensions[i];
      }
  
      // Allocate field arrays 
      basis_.allocate(nMonomer);
      rgrid_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         basis_[i].allocate(nBasis);
         rgrid_[i].allocate(meshDimensions);
      }

      isAllocated_ = true;
   }

   /*
   * Set new w-field values.
   */
   template <int D>
   void FieldContainer<D>::setBasis(DArray< DArray<double> > const & fields)
   {
      // Update system wFields
      for (int i = 0; i < nMonomer_; ++i) {
         DArray<double> const & f = fields[i];
         DArray<double> &  w = basis_[i];
         for (int j = 0; j < nBasis_; ++j) {
            w[j] = f[j];
         }
      }

      // Update system wFieldsRGrid
      fieldIoPtr_->convertBasisToRGrid(basis_, rgrid_);

      hasData_ = true;
      isSymmetric_ = true;
   }

   /*
   * Set new w-field values, using r-grid fields as inputs.
   */
   template <int D>
   void FieldContainer<D>::setRGrid(DArray< RField<D> > const & fields,
                                    bool isSymmetric)
   {
      // Update system wFieldsRGrid
      for (int i = 0; i < nMonomer_; ++i) {
         RField<D> const & f = fields[i];
         RField<D>& w = rgrid_[i];
         for (int j = 0; j < meshSize_; ++j) {
            w[j] = f[j];
         }
      }

      if (isSymmetric) {
         fieldIoPtr_->convertRGridToBasis(rgrid_, basis_);
      }

      hasData_ = true;
      isSymmetric_ =  isSymmetric;
   }

} // namespace Pspc
} // namespace Pscf
#endif
