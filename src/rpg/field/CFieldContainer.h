#ifndef RPG_C_FIELD_CONTAINER_H
#define RPG_C_FIELD_CONTAINER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/RField.h>             // member template parameter
#include <util/containers/DArray.h>       // member template

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

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
   * \ingroup Rpg_Field_Module
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
      * Set stored value of nMonomer.
      * 
      * May only be called once.
      *
      * \param nMonomer number of monomer types.
      */
      void setNMonomer(int nMonomer);

      /**
      * Allocate or re-allocate memory for fields in rgrid format.
      *
      * \param dimensions  dimensions of spatial mesh
      */
      void allocateRGrid(IntVec<D> const & dimensions);

      /**
      * De-allocate fields in rgrid format.
      */
      void deallocateRGrid();

      /**
      * Allocate or re-allocate memory for fields in basis format.
      *
      * \param nBasis  number of basis functions 
      */
      void allocateBasis(int nBasis);

      /**
      * De-allocate fields in basis format.
      */
      void deallocateBasis();

      /**
      * Allocate memory for both r-grid and basis field formats.
      *
      * This function may only be called once.
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
      {  
         UTIL_ASSERT(isAllocatedBasis_);
         return basis_; 
      }

      /**
      * Get the field for one monomer type in basis format (non-const).
      *
      * \param monomerId integer monomer type index (0, ... ,nMonomer-1)
      */
      DArray<double> & basis(int monomerId)
      {  
         UTIL_ASSERT(isAllocatedBasis_);
         return basis_[monomerId]; 
      }

      /**
      * Get the field for one monomer type in basis format (const)
      *
      * \param monomerId integer monomer type index (0, ... ,nMonomer-1)
      */
      DArray<double> const & basis(int monomerId) const
      {  
         UTIL_ASSERT(isAllocatedBasis_);
         return basis_[monomerId]; 
      }

      /**
      * Get array of all fields in r-grid format (non-const).
      */
      DArray< RField<D> > & rgrid()
      {
         UTIL_ASSERT(isAllocatedRGrid_);
         return rgrid_; 
      }

      /**
      * Get array of all fields in r-grid format (const).
      */
      DArray< RField<D> > const & rgrid() const
      {  
         UTIL_ASSERT(isAllocatedRGrid_);
         return rgrid_; 
      }

      /**
      * Get field for one monomer type in r-grid format (non-const)
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      RField<D> & rgrid(int monomerId)
      {  
         UTIL_ASSERT(isAllocatedRGrid_);
         return rgrid_[monomerId]; 
      }

      /**
      * Get field for one monomer type in r-grid format (const).
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      RField<D> const & rgrid(int monomerId) const
      {  
         UTIL_ASSERT(isAllocatedRGrid_);
         return rgrid_[monomerId]; 
      }

      /**
      * Has memory been allocated for fields in r-grid format?
      */
      bool isAllocatedRGrid() const
      {  return isAllocatedRGrid_; }

      /**
      * Has memory been allocated for fields in basis format?
      */
      bool isAllocatedBasis() const
      {  return isAllocatedBasis_; }

   private:

      /*
      * Array of fields in symmetry-adapted basis format
      *
      * Element basis_[i] is an array that contains the components
      * of the field associated with monomer i, in a symmetry-adapted
      * Fourier basis expansion. 
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
      * Number of monomer types.
      */
      int nMonomer_;

      /*
      * Has memory been allocated for fields in r-grid format?
      */
      bool isAllocatedRGrid_;

      /*
      * Has memory been allocated for fields in basis format?
      */
      bool isAllocatedBasis_;

   };

   #ifndef RPG_FIELD_CONTAINER_TPP
   // Suppress implicit instantiation
   extern template class CFieldContainer<1>;
   extern template class CFieldContainer<2>;
   extern template class CFieldContainer<3>;
   #endif

} // namespace Rpg
} // namespace Pscf
#endif
