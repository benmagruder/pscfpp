#ifndef PSPC_C_FIELD_CONTAINER_H
#define PSPC_C_FIELD_CONTAINER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>        // member template
#include <pspc/field/RField.h>             // member template parameter

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /**
   * A list of c fields stored in both basis and r-grid format.
   *
   * A CFieldContainer<D> contains representations of a list of nMonomer
   * fields that are associated with different monomer types in two 
   * different related formats:
   * 
   *  - A DArray of DArray<double> containers holds components of each
   *    field in a symmetry-adapted Fourier expansion (i.e., in basis 
   *    format). This is accessed by the basis() and basis(int) 
   *    member functions.
   *
   *  - A DArray of RField<D> containers holds valus of each field on
   *    the nodes of a regular grid. This is accessed by the rgrid()
   *    and rgrid(int) member functions.
   *
   * \ingroup Pspc_Field_Module
   */
   template <int D>
   class CFieldContainer 
   {

   public:

      /**
      * Constructor.
      */
      CFieldContainer();

      /**
      * Destructor.
      */
      ~CFieldContainer();

      /**
      * Allocate memory for fields.
      *
      * A CFieldContainer<D> may only be allocated once. An Exception will
      * be thrown if this function is called more than once.
      *
      * \param nMonomer  number of monomer types
      * \param nBasis  number of basis functions 
      * \param dimensions  dimensions of spatial mesh
      */
      void allocate(int nMonomer, int nBasis, IntVec<D> const & dimensions);

      /**
      * Get array of all fields in basis format (non-const).
      */
      DArray< DArray<double> > & basis()
      {  return basis_; }

      /**
      * Get array of all fields in basis format (const)
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< DArray<double> > const & basis() const
      {  return basis_; }

      /**
      * Get the field for one monomer type in basis format (non-const).
      *
      * \param monomerId integer monomer type index (0, ... ,nMonomer-1)
      */
      DArray<double> & basis(int monomerId)
      {  return basis_[monomerId]; }

      /**
      * Get the field for one monomer type in basis format (const)
      *
      * \param monomerId integer monomer type index (0, ... ,nMonomer-1)
      */
      DArray<double> const & basis(int monomerId) const
      {  return basis_[monomerId]; }

      /**
      * Get array of all fields in r-grid format (non-const).
      */
      DArray< RField<D> > & rgrid()
      {  return rgrid_; }

      /**
      * Get array of all fields in r-grid format (const).
      */
      DArray< RField<D> > const & rgrid() const
      {  return rgrid_; }

      /**
      * Get field for one monomer type in r-grid format (non-const)
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      RField<D> & rgrid(int monomerId)
      {  return rgrid_[monomerId]; }

      /**
      * Get field for one monomer type in r-grid format (const).
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      RField<D> const & rgrid(int monomerId) const
      {  return rgrid_[monomerId]; }

      /**
      * Has memory been allocated?
      */
      bool isAllocated() const
      {  return isAllocated_; }

   private:

      /*
      * Array of fields in symmetry-adapted basis format
      *
      * Element basis_[i] is an array that contains the components
      * of the field associated with monomer i, in a symmetry-adapted
      * Fourier expansion. 
      */
      DArray< DArray<double> > basis_;

      /*
      * Array of fields in real-space grid (r-grid) format
      *
      * Element basis_[i] is an RField<D> that contains values of the 
      * field associated with monomer i on the nodes of a regular mesh.
      */
      DArray< RField<D> > rgrid_;

      /*
      * Has memory been allocated for fields?
      */
      bool isAllocated_;

   };

   #ifndef PSPC_FIELD_CONTAINER_TPP
   // Suppress implicit instantiation
   extern template class CFieldContainer<1>;
   extern template class CFieldContainer<2>;
   extern template class CFieldContainer<3>;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif
