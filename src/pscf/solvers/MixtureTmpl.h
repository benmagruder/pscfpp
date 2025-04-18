#ifndef PSCF_MIXTURE_TMPL_H
#define PSCF_MIXTURE_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Monomer.h>
#include <util/param/ParamComposite.h>
#include <util/containers/DArray.h>

namespace Pscf
{

   using namespace Util;

   /**
   * A mixture of polymer and solvent species.
   *
   * \ingroup Pscf_Solver_Module
   */
   template <class TP, class TS>
   class MixtureTmpl : public ParamComposite
   {
   public:

      // Public typedefs

      /**
      * Polymer species solver typename.
      */
      typedef TP Polymer;

      /**
      * Solvent species solver typename.
      */
      typedef TS Solvent;

      // Public member functions

      /**
      * Constructor.
      */
      MixtureTmpl();

      /**
      * Destructor.
      */
      ~MixtureTmpl();

      /**
      * Read parameters from file and initialize.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      /// \name Accessors (by non-const reference)
      ///@{
 
      /**
      * Get a Monomer type descriptor (const reference).
      *
      * \param id integer monomer type index (0 <= id < nMonomer)
      */
      Monomer const & monomer(int id) const;

      /**
      * Get a polymer object.
      *
      * \param id integer polymer species index (0 <= id < nPolymer)
      */
      Polymer& polymer(int id);

      /**
      * Get a polymer object by const reference.
      *
      * \param id integer polymer species index (0 <= id < nPolymer)
      */
      Polymer const & polymer(int id) const;

      /**
      * Set a solvent solver object.
      *
      * \param id integer solvent species index (0 <= id < nSolvent)
      */
      Solvent& solvent(int id);

      /**
      * Set a solvent solver object.
      *
      * \param id integer solvent species index (0 <= id < nSolvent)
      */
      Solvent const & solvent(int id) const;

      ///@}
      /// \name Accessors (by value)
      ///@{
 
      /**
      * Get number of monomer types.
      */
      int nMonomer() const;

      /**
      * Get number of polymer species.
      */
      int nPolymer() const;

      /** 
      * Get number of total blocks in the mixture across all polymers.
      */
      int nBlock() const;

      /**
      * Get number of solvent (point particle) species.
      */
      int nSolvent() const;

      /**
      * Get monomer reference volume (set to 1.0 by default).
      */
      double vMonomer() const;
      
      /**
      * Set value of monomer reference volume.
      * 
      * An initial value is normally read from a parameter file. 
      * This function is provided for use by a ramp in which 
      * monomer reference volume is modified after initialization.
      * 
      * \param vMonomer new monomer reference volume.
      */
      void setVmonomer(double vMonomer);

      ///@}

   protected:

      /**
      * Get a Monomer type descriptor (non-const reference).
      *
      * \param id integer monomer type index (0 <= id < nMonomer)
      */
      Monomer& monomer(int id);

   private:

      /**
      * Array of monomer type descriptors.
      */
      DArray<Monomer> monomers_;

      /**
      * Array of polymer species solver objects.
      *
      * Array capacity = nPolymer.
      */
      DArray<Polymer> polymers_;

      /**
      * Array of solvent species objects.
      *
      * Array capacity = nSolvent.
      */
      DArray<Solvent> solvents_;

      /**
      * Number of monomer types.
      */
      int nMonomer_; 

      /**
      * Number of polymer species.
      */
      int nPolymer_;

      /**
      * Number of solvent species.
      */
      int nSolvent_;

      /**
      * Number of blocks total, across all polymers.
      */
      int nBlock_;

      /**
      * Monomer reference volume (set to 1.0 by default).
      */
      double vMonomer_;

   };

   // Inline member functions

   template <class TP, class TS>
   inline int MixtureTmpl<TP,TS>::nMonomer() const
   {  return nMonomer_; }

   template <class TP, class TS>
   inline int MixtureTmpl<TP,TS>::nPolymer() const
   {  return nPolymer_; }

   template <class TP, class TS>
   inline int MixtureTmpl<TP,TS>::nSolvent() const
   {  return nSolvent_; }

   template <class TP, class TS>
   inline int MixtureTmpl<TP,TS>::nBlock() const
   {  return nBlock_; }

   template <class TP, class TS>
   inline Monomer const & MixtureTmpl<TP,TS>::monomer(int id) const
   {  
      UTIL_CHECK(id < nMonomer_);
      return monomers_[id]; 
   }

   template <class TP, class TS>
   inline Monomer& MixtureTmpl<TP,TS>::monomer(int id)
   {  
      UTIL_CHECK(id < nMonomer_);
      return monomers_[id]; 
   }

   template <class TP, class TS>
   inline TP& MixtureTmpl<TP,TS>::polymer(int id)
   {  
      UTIL_CHECK(id < nPolymer_);
      return polymers_[id];
   }

   template <class TP, class TS>
   inline TP const & MixtureTmpl<TP,TS>::polymer(int id) const
   {  
      UTIL_CHECK(id < nPolymer_);
      return polymers_[id];
   }

   template <class TP, class TS>
   inline TS& MixtureTmpl<TP,TS>::solvent(int id)
   {  
      UTIL_CHECK(id < nSolvent_);
      return solvents_[id]; 
   }

   template <class TP, class TS>
   inline TS const & MixtureTmpl<TP,TS>::solvent(int id) const
   {  
      UTIL_CHECK(id < nSolvent_);
      return solvents_[id]; 
   }

   template <class TP, class TS>
   inline double MixtureTmpl<TP,TS>::vMonomer() const
   {  return vMonomer_; }
   
   template <class TP, class TS>
   inline void MixtureTmpl<TP,TS>::setVmonomer(double vMonomer)
   { vMonomer_ = vMonomer; }
   
   // Non-inline member functions

   /*
   * Constructor.
   */
   template <class TP, class TS>
   MixtureTmpl<TP,TS>::MixtureTmpl()
    : ParamComposite(),
      monomers_(),
      polymers_(),
      solvents_(),
      nMonomer_(0), 
      nPolymer_(0),
      nSolvent_(0),
      nBlock_(0),
      vMonomer_(1.0)
   {}

   /*
   * Destructor.
   */
   template <class TP, class TS>
   MixtureTmpl<TP,TS>::~MixtureTmpl()
   {}

   /*
   * Read all parameters and initialize.
   */
   template <class TP, class TS>
   void MixtureTmpl<TP,TS>::readParameters(std::istream& in)
   {
      // Read nMonomer and monomers array
      read<int>(in, "nMonomer", nMonomer_);
      monomers_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         monomers_[i].setId(i);
      }
      readDArray< Monomer >(in, "monomers", monomers_, nMonomer_);

      /*
      * The input format for a single monomer is defined in the istream
      * extraction operation (operator >>) for a Pscf::Monomer, in file
      * pscf/chem/Monomer.cpp. The text representation contains only the
      * value for the monomer statistical segment Monomer::kuhn.
      */

      // Read nPolymer
      read<int>(in, "nPolymer", nPolymer_);
      UTIL_CHECK(nPolymer_ > 0);

      // Optionally read nSolvent, with nSolvent=0 by default
      nSolvent_ = 0;
      readOptional<int>(in, "nSolvent", nSolvent_);

      // Read polymers and compute nBlock
      nBlock_ = 0;
      if (nPolymer_ > 0) {

         polymers_.allocate(nPolymer_);
         for (int i = 0; i < nPolymer_; ++i) {
            readParamComposite(in, polymer(i));
            nBlock_ = nBlock_ + polymer(i).nBlock();
         }
   
         // Set statistical segment lengths for all blocks
         double kuhn;
         int monomerId;
         for (int i = 0; i < nPolymer_; ++i) {
            for (int j = 0; j < polymer(i).nBlock(); ++j) {
               monomerId = polymer(i).block(j).monomerId();
               kuhn = monomer(monomerId).kuhn();
               polymer(i).block(j).setKuhn(kuhn);
            }
         }

      }

      // Read solvents
      if (nSolvent_ > 0) {

         solvents_.allocate(nSolvent_);
         for (int i = 0; i < nSolvent_; ++i) {
            readParamComposite(in, solvent(i));
         }

      }

      // Optionally read monomer reference value
      vMonomer_ = 1.0; // Default value
      readOptional(in, "vMonomer", vMonomer_);

   }

}
#endif
