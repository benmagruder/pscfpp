#ifndef RPC_DOMAIN_H
#define RPC_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class

#include <rpc/field/FieldIo.h>            // member

#include <prdc/cpu/FFT.h>                  // member
#include <prdc/crystal/Basis.h>            // member
#include <prdc/crystal/SpaceGroup.h>       // member
#include <prdc/crystal/UnitCell.h>         // member

#include <pscf/mesh/Mesh.h>                // member

#include <string>

namespace Pscf {
namespace Rpc
{

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Spatial domain and spatial discretization for a periodic structure.
   *
   * A Domain has (among other components):
   *
   *    - a Pscf::Mesh spatial discretization mesh
   *    - a Prdc::UnitCell crystallographic unit cell
   *    - a Prdc::SpaceGroup crystallographic space group
   *    - a Prdc::Basis symmetry-adapated Fourier basis
   *    - a Prdc::Cpu::FFT Fast Fourier Transform
   *    - a Rpc::FieldIo object for field IO and conversion operations
   *    - a lattice system enum (type Prdc::UnitCell\<D\>::LatticeSystem)
   *    - a groupName string
   *
   * Note: Class names Pscf::Mesh, Prdc::UnitCell, etc. mentioned above are
   * actually all class templates with a template parameter D. Actual class
   * names are Mesh \<D\>, Prdc::UnitCell \<D\>, etc. with D=1, 2, or 3.
   *
   * \see \ref rpc_Domain_page "Manual Page"
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class Domain : public ParamComposite
   {

   public:


      /**
      * Constructor.
      */
      Domain();

      /**
      * Destructor.
      */
      ~Domain();

      /// \name Initialization and Mutators
      ///@{

      /**
      * Create association with a FileMaster, needed by FieldIo.
      *
      * \param fileMaster associated FileMaster object.
      */
      void setFileMaster(FileMaster& fileMaster);

      /**
      * Read body of parameter block (without opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Read header of an r-grid field file to initialize this Domain.
      *
      * Useful for unit testing.
      *
      * \param in input parameter stream
      * \param nMonomer number of monomers in field file (output)
      */
      void readRGridFieldHeader(std::istream& in, int& nMonomer);

      /**
      * Set unit cell. 
      *
      * The lattice system in the unitCell must match any value that was
      * read from the parameter file. This function initializes the basis if 
      * a space group exists and the basis was not initialized previously.
      *
      * \param unitCell new unit cell
      */
      void setUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set unit cell state.
      *
      * The "lattice" lattice system enumeration match any non-null value that 
      * was read from the parameter file. This function initializes the basis 
      * if a space group exists and the basis was not initialized previously.
      *
      * \param lattice  lattice system
      * \param parameters array of unit cell parameters
      */
      void setUnitCell(typename UnitCell<D>::LatticeSystem lattice,
                       FSArray<double, 6> const & parameters);

      /**
      * Set unit cell parameters.
      *
      * The lattice system must already be set to non-null value on entry.
      * This function initializes the basis if a space group exists and the
      * basis was not initialized previously.
      *
      * \param parameters array of unit cell parameters
      */
      void setUnitCell(FSArray<double, 6> const & parameters);

      /**
      * Construct basis if not done already.
      */
      void makeBasis();

      ///@}
      /// \name Accessors (return by non-const or const reference)
      ///@{

      /**
      * Get UnitCell (i.e., lattice type and parameters) by reference.
      */
      UnitCell<D>& unitCell();

      /**
      * Get UnitCell (i.e., lattice type and parameters) by reference.
      */
      UnitCell<D> const & unitCell() const;

      /**
      * Get spatial discretization mesh by reference.
      */
      Mesh<D>& mesh();

      /**
      * Get spatial discretization mesh by const reference.
      */
      Mesh<D> const & mesh() const;

      /**
      * Get associated SpaceGroup object by const reference.
      */
      SpaceGroup<D> const & group() const ;

      /**
      * Get associated Basis object by reference.
      */
      Basis<D>& basis();

      /**
      * Get associated Basis object by const reference.
      */
      Basis<D> const & basis() const ;

      /**
      * Get associated FFT object.
      */
      FFT<D>& fft();

      /**
      * Get associated FFT object by const reference.
      */
      FFT<D> const & fft() const;

      /**
      * Get associated FieldIo object.
      */
      FieldIo<D>& fieldIo();

      /**
      * Get associated FieldIo object by const reference.
      */
      FieldIo<D> const & fieldIo() const;

      ///@}
      /// \name Accessors (return by value)
      ///@{

      /** 
      * Get lattice system.
      */  
      typename UnitCell<D>::LatticeSystem lattice() const;

      /** 
      * Get group name.
      */  
      std::string groupName() const;

      /**
      * Has a space group been identified?
      */ 
      bool hasGroup() const;

      /**
      * Has a symmetry-adapted Fourier basis been initialized?
      */ 
      bool hasBasis() const;

      ///@}

   private:

      // Private member variables

      /**
      * Crystallographic unit cell (crystal system and cell parameters).
      */
      UnitCell<D> unitCell_;

      /**
      * Spatial discretization mesh.
      */
      Mesh<D> mesh_;

      /**
      * SpaceGroup object
      */
      SpaceGroup<D> group_;

      /**
      * Basis object.
      */
      Basis<D> basis_;

      /**
      * FFT object to be used by solvers.
      */
      FFT<D> fft_;

      /**
      * FieldIo object for field input/output operations
      */
      FieldIo<D> fieldIo_;

      /**
      * Lattice system (enumeration value).
      */
      typename UnitCell<D>::LatticeSystem lattice_;

      /**
      * Group name.
      */
      std::string groupName_;

      /**
      * Has a space group been indentified?
      */
      bool hasGroup_;

      /**
      * Has a FileMaster been set?
      */
      bool hasFileMaster_;

      /**
      * Has this Domain object been initialized?
      */
      bool isInitialized_;

      // members of parent class with non-dependent names
      using ParamComposite::read;
      using ParamComposite::readOptional;

   };

   // Inline member functions

   // Get the UnitCell<D> object by non-const reference.
   template <int D>
   inline UnitCell<D>& Domain<D>::unitCell()
   {  return unitCell_; }

   // Get the UnitCell<D> object by const reference.
   template <int D>
   inline UnitCell<D> const & Domain<D>::unitCell() const
   {  return unitCell_; }

   // Get the Mesh<D> object by reference.
   template <int D>
   inline Mesh<D>& Domain<D>::mesh()
   {  return mesh_; }

   // Get the Mesh<D> object by const reference.
   template <int D>
   inline Mesh<D> const & Domain<D>::mesh() const
   {  return mesh_; }

   // Get the SpaceGroup<D> object by const reference.
   template <int D>
   inline SpaceGroup<D> const & Domain<D>::group() const
   {  return group_; }

   // Get the Basis<D> object by non-const reference.
   template <int D>
   inline Basis<D>& Domain<D>::basis()
   {  return basis_; }

   // Get the Basis<D> object by const reference.
   template <int D>
   inline Basis<D> const & Domain<D>::basis() const
   {  return basis_; }

   // Get the FFT<D> object.
   template <int D>
   inline FFT<D>& Domain<D>::fft()
   {  return fft_; }

   // Get the FFT<D> object.
   template <int D>
   inline FFT<D> const & Domain<D>::fft() const
   {  return fft_; }

   // Get the FieldIo<D> object.
   template <int D>
   inline FieldIo<D>& Domain<D>::fieldIo()
   {  return fieldIo_; }

   // Get the FieldIo<D> object by const reference.
   template <int D>
   inline FieldIo<D> const & Domain<D>::fieldIo() const
   {  return fieldIo_; }

   // Get the lattice system enumeration value
   template <int D>
   inline 
   typename UnitCell<D>::LatticeSystem Domain<D>::lattice() 
   const
   {  return lattice_; }

   // Get the groupName string.
   template <int D>
   inline std::string Domain<D>::groupName() const
   {  return groupName_; }

   // Has a space group been identified?
   template <int D>
   inline bool Domain<D>::hasGroup() const
   {  return hasGroup_; }

   // Has a symmetry-adapted Fourier basis been initialized ?
   template <int D>
   inline bool Domain<D>::hasBasis() const
   {  return basis_.isInitialized(); }

   #ifndef RPC_DOMAIN_TPP
   // Suppress implicit instantiation
   extern template class Domain<1>;
   extern template class Domain<2>;
   extern template class Domain<3>;
   #endif

} // namespace Rpc
} // namespace Pscf
#endif
