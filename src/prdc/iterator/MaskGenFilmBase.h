#ifndef PRDC_MASK_GEN_FILM_BASE_H
#define PRDC_MASK_GEN_FILM_BASE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/FieldGenerator.h>  // Base class
#include <pscf/math/RealVec.h>        // container
#include <prdc/crystal/UnitCell.h>
#include <util/containers/FSArray.h>  // container
#include <iostream>
#include <string>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Base class defining mask that imposes thin film confinement.
   * 
   * This is a base class for MaskGenFilm that defines all traits of a 
   * MaskGenFilm that do not require access to the System (System access is
   * needed, for example, to get the space group and set the mask field).
   * 
   * If the user chooses a MaskGenFilm object to construct the mask, then 
   * the system will contain two parallel hard surfaces ("walls"), confining
   * the polymers/solvents to a "thin film" region of the unit cell. The
   * shape of the mask is defined by three input parameters: normalVecId,
   * excludedThickness, and interfaceThickness. See \ref 
   * scft_thin_films_page for more information. 
   * 
   * \ingroup Prdc_Iterator_Module
   */
   template <int D>
   class MaskGenFilmBase : public FieldGenerator
   {

   public:

      /**
      * Constructor.
      */
      MaskGenFilmBase();

      /**
      * Destructor.
      */
      ~MaskGenFilmBase();

      /**
      * Read parameter file block and initialize.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Check that the system is compatible with this field.
      * 
      * This method calls setFlexibleParams, checkLatticeVectors, and
      * checkSpaceGroup.
      */
      void checkCompatibility();

      /**
      * Check whether system has changed such that the field needs updating.
      */
      bool updateNeeded() const;

      /**
      * Get value of normalVecId.
      */
      int normalVecId() const;

      /**
      * Get value of interfaceThickness.
      */
      double interfaceThickness() const;

      /**
      * Get value of excludedThickness.
      */
      double excludedThickness() const;

      /**
      * Check whether a value of fBulk was provided.
      */
      bool hasFBulk() const;

      /**
      * Check whether the field has been generated.
      */
      bool isGenerated() const = 0;

   protected:

      /**
      * Check that space group is compatible with the mask.
      */
      void checkSpaceGroup() const;

      /**
      * Check that lattice vectors are compatible with thin film constraint.
      * 
      * Check that user-defined lattice basis vectors (stored in the
      * Domain member of the parent System object) are compatible with 
      * thin film confinement. The lattice basis vector with index 
      * normalVecId should be normal to the walls, while any other lattice
      * basis vectors must be parallel to the walls.
      */
      virtual void checkLatticeVectors() const;

      /**
      * Allocate container necessary to generate and store field.
      */ 
      virtual void allocate() = 0;

      /**
      * Generate the field and store where the Iterator can access.
      */
      virtual void generate() = 0;

      /**
      * Modify stress value.
      * 
      * If the lattice parameter corresponds to normalVecId (and 
      * therefore defines the film thickness), then the stress is 
      * calculated using the excess free energy per unit area, rather 
      * than the absolute free energy. This requires a slight
      * modification of the equation used to calculate stress, as
      * implemented in this method. For all other lattice parameters, 
      * the stress returned by this method is the same value that was 
      * passed in as an input.
      * 
      * Note that, by default, the lattice parameter corresponding to 
      * normalVecId should not be flexible, so this method should not be
      * called for that lattice parameter in most instances. In order to 
      * calculate the excess free energy per unit area, a reference free 
      * energy must be provided equal to the free energy of the bulk 
      * phase corresponding to this thin film morphology. Users can allow
      * the film thickness to be flexible and optimized by providing a 
      * reference free energy via the optional input parameter fBulk. 
      * If fBulk is provided, the lattice parameter corresponding to 
      * normalVecId will be flexible unless the user specified that it is
      * rigid in the Iterator parameters.
      * 
      * This method is left virtual here because it requires a significant
      * amount of information from the system. Therefore, it is 
      * implemented by subclasses that have System access.
      * 
      * \param paramId  index of the lattice parameter with this stress
      * \param stress  stress value calculated by Mixture object
      */
      virtual double modifyStressValue(int paramId, double stress) const = 0;

      /** Which lattice parameter corresponds to normalVecId?
      * 
      * The index normalVecId refers to one lattice basis vector (e.g.,
      * in 3D the value is 0, 1, or 2). However, for basis vector i, the 
      * vector length is not always the i-th lattice parameter in the 
      * list of lattice parameters. For example, if normalVecId = 2 but
      * the lattice system is tetragonal, then the length of the vector
      * denoted by this normalVecId has index 1 in the list of lattice
      * parameters. This method converts a normalVecId index into the
      * corresponding lattice parameter index.
      * 
      * \param normalVecIndex the value of normalVecId to convert
      */
      int convertNormalVecIdToParamId(int normalVecIndex) const;

      /**
      * Sets flexible lattice parameters to be compatible with the mask.
      * 
      * An iterator for a thin film SCFT calculation should allow for
      * some lattice parameters to be fixed, while others are held 
      * constant. Subclasses should define this method to set the 
      * flexibility of these lattice parameters appropriately so that 
      * the film thickness is held constant, as are all but one of the 
      * angles between basis vectors (the angle in the plane of the film
      * may vary). The other lattice parameters are allowed to be 
      * flexible if the user specified that they are flexible in the 
      * Iterator parameters. 
      * 
      * Note that the lattice parameter that defines the film thickness
      * may be flexible if the optional input parameter fBulk is provided.
      * This parameter is necessary to compute the stress in the direction
      * normal to the film.
      */
      virtual void setFlexibleParams() = 0;

      /**
      * Get the space group name for this system.
      */
      virtual std::string systemSpaceGroup() const = 0;

      /**
      * Get the lattice system for this system.
      */
      virtual 
      typename UnitCell<D>::LatticeSystem systemLatticeSystem() const = 0;

      /**
      * Get one of the lattice vectors for this system.
      * 
      * \param id  index of the desired lattice vector
      */
      virtual RealVec<D> systemLatticeVector(int id) const = 0;

      /**
      * The lattice vector normal to the film used to generate these fields.
      * 
      * This vector is set to be equal to the system's lattice vector with
      * index normalVecId_ each time the external fields are generated. The 
      * system's lattice vectors may then change, and this normalVecCurrent_
      * vector is used to detect whether they have changed. This is used to 
      * decide whether a new set of external fields needs to be generated.
      */
      RealVec<D> normalVecCurrent_;

      /// Reference free energy used to calculate stress normal to the film
      double fBulk_;

      using ParamComposite::read;
      using ParamComposite::readOptional;

   private:

      /// Lattice basis vector that is normal to the walls
      int normalVecId_;

      /// Interface thickness
      double interfaceThickness_;

      /// Excluded (wall) thickness
      double excludedThickness_;

      /// Does this object have a value of fBulk from the parameter file?
      bool hasFBulk_;

      using FieldGenerator::type_;

   };

   // Inline member functions

   // Get value of normalVecId.
   template <int D> 
   inline int MaskGenFilmBase<D>::normalVecId() const
   {  return normalVecId_; }

   // Get value of interfaceThickness.
   template <int D> 
   inline double MaskGenFilmBase<D>::interfaceThickness() const
   {  return interfaceThickness_; }

   // Get value of excludedThickness.
   template <int D> 
   inline double MaskGenFilmBase<D>::excludedThickness() const
   {  return excludedThickness_; }

   // Check whether a value of fBulk was provided.
   template <int D> 
   inline bool MaskGenFilmBase<D>::hasFBulk() const
   {  return hasFBulk_; }

}
}
#endif
