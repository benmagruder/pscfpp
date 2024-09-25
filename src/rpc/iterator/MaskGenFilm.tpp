#ifndef RPC_MASK_GEN_FILM_TPP
#define RPC_MASK_GEN_FILM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaskGenFilm.h"
#include <rpc/field/FieldIo.h>
#include <prdc/cpu/RField.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/paramIdConversions.h>
#include <pscf/math/RealVec.h>
#include <pscf/math/IntVec.h>
#include <util/containers/FArray.h>
#include <cmath>

namespace Pscf {
namespace Rpc
{

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Default constructor
   */
   template <int D>
   MaskGenFilm<D>::MaskGenFilm()
    : MaskGenFilmBase<D>::MaskGenFilmBase(),
      sysPtr_(0)
   {  setClassName("MaskGenFilm"); }
   
   /*
   * Constructor
   */
   template <int D>
   MaskGenFilm<D>::MaskGenFilm(System<D>& sys)
    : MaskGenFilmBase<D>::MaskGenFilmBase(),
      sysPtr_(&sys)
   {  setClassName("MaskGenFilm"); }

   /*
   * Destructor
   */
   template <int D>
   MaskGenFilm<D>::~MaskGenFilm()
   {}

   /*
   * Get contribution to the stress from this mask
   */
   template <int D>
   double MaskGenFilm<D>::stressTerm(int paramId) const
   {
      int normalVecParamId = convertFullParamIdToReduced<D>(normalVecId(),
                                                   system().domain().lattice());

      if (normalVecParamId == paramId) {
         UTIL_CHECK(system().hasMask());
         
         // Get the length L of the lattice basis vector normal to the walls
         double L = system().domain().unitCell().parameter(paramId);

         // Create a 3 element vector 'dim' that contains the grid dimensions.
         // If system is 2D (1D), then the z (y & z) dimensions are set to 1.
         IntVec<3> dim;
         for (int ind = 0; ind < 3; ind++) {
            if (ind < D) {
               dim[ind] = system().domain().mesh().dimensions()[ind];
            } else {
               dim[ind] = 1;
            }
         }

         // Create the field d\phi_m / dL in rgrid format
         RField<D> rGrid;
         FArray<int,3> coords;
         int x, y, z;
         int counter = 0;
         double d, maskVal;

         rGrid.allocate(system().domain().mesh().dimensions());
         RField<D> const & maskRGrid = system().mask().rgrid();
         
         for (x = 0; x < dim[0]; x++) {
            coords[0] = x;
            for (y = 0; y < dim[1]; y++) {
               coords[1] = y;
               for (z = 0; z < dim[2]; z++) {
                  coords[2] = z;
                  maskVal = maskRGrid[counter];

                  // Get the distance 'd' traveled along the lattice basis 
                  // vector normal to the walls, in reduced coordinates
                  d = (double)coords[normalVecId()] / 
                      (double)dim[normalVecId()];

                  rGrid[counter] = maskVal * (maskVal - 1) * 8
                                   * (std::abs(d - 0.5) - 0.5)
                                   / interfaceThickness();
                  counter++;
               }
            }
         }

         // Convert above field into basis format
         DArray<double> basis;
         int nBasis = system().basis().nBasis();
         int nMonomer = system().mixture().nMonomer();
         basis.allocate(nBasis);
         system().fieldIo().convertRGridToBasis(rGrid, basis);

         // Get the integral term in the stress
         double intTerm = 0;
         DArray<double> xi;
         xi.allocate(nBasis);
         DArray<double> wVals;
         wVals.allocate(nMonomer);

         if (system().hasExternalFields()) {
            for (int i = 0; i < nBasis; i++) {
               xi[i] = system().w().basis(0)[i] - system().h().basis(0)[i];
               for (int j = 1; j < nMonomer; j++) {
                  xi[i] -= system().c().basis(j)[i] * 
                        system().interaction().chi(0,j);
               }
               intTerm += xi[i] * basis[i];
            }
         } else {
            for (int i = 0; i < nBasis; i++) {
               xi[i] = system().w().basis(0)[i];
               for (int j = 1; j < nMonomer; j++) {
                  xi[i] -= system().c().basis(j)[i] * 
                        system().interaction().chi(0,j);
               }
               intTerm += xi[i] * basis[i];
            }
         }

         intTerm /= system().mask().phiTot();

         // Get the free energy term in the stress
         double fQTerm = 0;
         double phi, mu;
         int nPolymer = system().mixture().nPolymer();
         int nSolvent = system().mixture().nSolvent();
         
         if (nPolymer > 0) { // Subtract out polymer ideal gas contributions
            Polymer<D> const * polymerPtr;
            double length;
            for (int i = 0; i < nPolymer; ++i) {
               polymerPtr = &system().mixture().polymer(i);
               phi = polymerPtr->phi();
               mu = polymerPtr->mu();
               length = polymerPtr->length();
               // Recall: mu = ln(phi/q)
               if (phi > 1.0E-08) {
                  fQTerm += phi*( mu - 1.0 )/length;
               }
            }
         }
         if (nSolvent > 0) { // Subtract out solvent ideal gas contributions
            Solvent<D> const * solventPtr;
            double size;
            for (int i = 0; i < nSolvent; ++i) {
               solventPtr = &system().mixture().solvent(i);
               phi = solventPtr->phi();
               mu = solventPtr->mu();
               size = solventPtr->size();
               if (phi > 1.0E-08) {
                  fQTerm += phi*( mu - 1.0 )/size;
               }
            }
         }

         // Get system free energy
         if (!sysPtr_->hasFreeEnergy()) {
            sysPtr_->computeFreeEnergy();
         }
         double fSys = sysPtr_->fHelmholtz();

         double fTerm = (fSys - fQTerm) * excludedThickness() / 
                        (system().mask().phiTot() * L * L);

         double term = - intTerm - fTerm;
         return term;

      } else {

         return 0.0;

      }
   }

   /*
   * Allocate container necessary to generate and store field
   */
   template <int D>
   void MaskGenFilm<D>::allocate()
   {
      UTIL_CHECK(system().basis().isInitialized());
      UTIL_CHECK(system().unitCell().isInitialized());

      system().mask().setFieldIo(system().fieldIo());

      // Allocate the mask containers if needed
      if (!system().mask().isAllocated()) {
         system().mask().allocate(system().basis().nBasis(), 
                                  system().mesh().dimensions());
      }
   }

   /*
   * Generate the field and store where the Iterator can access
   */
   template <int D>
   void MaskGenFilm<D>::generate()
   {
      UTIL_CHECK(interfaceThickness() > 0);
      UTIL_CHECK(excludedThickness() > interfaceThickness());
      UTIL_CHECK(system().mask().isAllocated());

      // Get the length L of the lattice basis vector normal to the walls
      int paramId = convertFullParamIdToReduced<D>(normalVecId(),
                                                system().domain().lattice());
      double L = system().domain().unitCell().parameter(paramId);

      // Create a 3 element vector 'dim' that contains the grid dimensions.
      // If system is 2D (1D), then the z (y & z) dimensions are set to 1.
      IntVec<3> dim;
      for (int ind = 0; ind < 3; ind++) {
         if (ind < D) {
            dim[ind] = system().domain().mesh().dimensions()[ind];
         } else {
            dim[ind] = 1;
         }
      }

      // Generate an r-grid representation of the walls
      RField<D> rGrid;
      rGrid.allocate(system().domain().mesh().dimensions());
      int x, y, z;
      int counter = 0;
      FArray<int,3> coords;
      double d, rhoW;

      for (x = 0; x < dim[0]; x++) {
         coords[0] = x;
         for (y = 0; y < dim[1]; y++) {
            coords[1] = y;
            for (z = 0; z < dim[2]; z++) {
               coords[2] = z;

               // Get the distance 'd' traveled along the lattice basis 
               // vector that is normal to the walls
               d = coords[normalVecId()] * L / dim[normalVecId()];

               // Calculate wall volume fraction (rhoW) at gridpoint (x,y,z)
               rhoW = 0.5 * (1 + tanh(4 * (((.5 * (excludedThickness()-L)) + 
                                  fabs(d - (L/2))) / interfaceThickness())));
               rGrid[counter++] = 1-rhoW;
            }
         }
      }

      // Store this mask in System
      system().mask().setRGrid(rGrid,true);

      // Store lattice vector normal to film used to construct this mask
      normalVecCurrent_ = systemLatticeVector(normalVecId());

   }

   template <int D>
   double MaskGenFilm<D>::modifyStressValue(int paramId, double stress) 
   const
   {
      int nvParamId = convertFullParamIdToReduced<D>(normalVecId(),
                                                 system().domain().lattice());

      if (nvParamId == paramId) {
         UTIL_CHECK(system().hasMask());

         // Get system free energy
         if (!sysPtr_->hasFreeEnergy()) {
            sysPtr_->computeFreeEnergy();
         }
         double fSys = sysPtr_->fHelmholtz();

         // Get the length L of the lattice basis vector normal to the walls
         double L = system().domain().unitCell().parameter(paramId);

         double modifiedStress = (stress * (L - excludedThickness())) 
                                 + (fSys - fBulk_);
         return modifiedStress;
      } else {
         return stress;
      }
   }

   // Explicit Specializations for setFlexibleParams are in MaskGenFilm.cpp
}
}
#endif