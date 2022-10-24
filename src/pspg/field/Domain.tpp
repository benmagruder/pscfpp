#ifndef PSPG_DOMAIN_TPP
#define PSPG_DOMAIN_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"

namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Domain<D>::Domain()
    : unitCell_(),
      mesh_(),
      basis_(),
      fft_(),
      fieldIo_(),
      lattice_(UnitCell<D>::Null),
      groupName_(),
      hasFileMaster_(false),
      isInitialized_(false)
   {  setClassName("Domain"); }

   /*
   * Destructor.
   */
   template <int D>
   Domain<D>::~Domain()
   {}

   template <int D>
   void Domain<D>::setFileMaster(FileMaster& fileMaster)
   {
      fieldIo_.associate(mesh_, fft_, 
                         lattice_, groupName_, group_, basis_, 
                         fileMaster);
      hasFileMaster_ = true;
   }

   /*
   * Read parameters and initialize.
   */
   template <int D>
   void Domain<D>::readParameters(std::istream& in)
   {
      UTIL_CHECK(hasFileMaster_);

      read(in, "unitCell", unitCell_);
      read(in, "mesh", mesh_);
      read(in, "groupName", groupName_);

      lattice_ = unitCell_.lattice();

      fft_.setup(mesh_.dimensions());

      // Initialize space group
      readGroup(groupName_, group_);

      // Make symmetry-adapted basis
      basis().makeBasis(mesh(), unitCell(), group_);

      isInitialized_ = true;
   }
   
 
   template <int D> 
   void Domain<D>::readFieldHeader(std::istream& in, int& nMonomer)
   {
      // Read common section of standard field header
      int ver1, ver2;
      Pscf::readFieldHeader(in, ver1, ver2, 
                            unitCell_, groupName_, nMonomer);
 
      // Read grid dimensions
      std::string label;
      in >> label;
      UTIL_CHECK(label == "ngrid");
      IntVec<D> nGrid;
      in >> nGrid;

      // Initialize mesh, fft and basis
      mesh_.setDimensions(nGrid);
      fft_.setup(mesh_.dimensions());
      basis_.makeBasis(mesh_, unitCell_, groupName_);
      
      isInitialized_ = true;
   }

} // namespace Pspg
} // namespace Pscf
#endif
