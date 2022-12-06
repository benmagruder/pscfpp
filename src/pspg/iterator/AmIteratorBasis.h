#ifndef PSPG_AM_ITERATOR_BASIS_H
#define PSPG_AM_ITERATOR_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <pscf/iterator/AmIteratorTmpl.h>                 
#include <pscf/iterator/AmbdInteraction.h>   // member variable

namespace Pscf {
namespace Pspg
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Pspg implementation of the Anderson Mixing iterator.
   *
   * \ingroup Pspg_Iterator_Module
   */
   template <int D>
   class AmIteratorBasis : public AmIteratorTmpl<Iterator<D>, DArray<double> >
   {

   public:

      /**
      * Constructor.
      *   
      * \param system parent system object
      */
      AmIteratorBasis(System<D>& system);

      /**
      * Destructor.
      */ 
      ~AmIteratorBasis();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      void readParameters(std::istream& in);

      using Iterator<D>::isFlexible;
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::solve;

   protected:

      using ParamComposite::readOptional;
      using Iterator<D>::system;
      using Iterator<D>::isFlexible_;
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::setClassName;
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::verbose;

      /**
      * Setup iterator just before entering iteration loop.
      *
      * \param isContinuation Is this a continuation within a sweep?
      */
      void setup(bool isContinuation);

   private:

      // Local copy of interaction, adapted for use AMBD residual definition
      AmbdInteraction interaction_;

      /// How are stress residuals scaled in error calculation?
      double scaleStress_;

      /**
      * Set a vector equal to another (assign a = b)
      * 
      * \param a the field to be set (LHS, result)
      * \param b the field for it to be set to (RHS, input)
      */
      void setEqual(DArray<double>& a, DArray<double> const & b);

      /**
      * Compute the inner product of two real vectors.
      */
      double dotProduct(DArray<double> const & a, DArray<double> const & b);

      /**
      * Find the maximum magnitude element of a vector.
      */
      double maxAbs(DArray<double> const & hist);

      /**
      * Update the series of residual vectors.
      * 
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of pase residual or field vectors
      */
      void updateBasis(RingBuffer< DArray<double> > & basis, 
                       RingBuffer< DArray<double> > const & hists);

      /**
      * Compute trial field so as to minimize L2 norm of residual.
      * 
      * \param trial resulting trial field (output)
      * \param basis RingBuffer of residual basis vectors.
      * \param coeffs coefficients of basis vectors
      * \param nHist number of prior states stored
      */
      void addHistories(DArray<double>& trial, 
                        RingBuffer<DArray<double> > const & basis, 
                        DArray<double> coeffs, int nHist);

      /**
      * Add predicted error to the trial field.
      * 
      * \param fieldTrial trial field (input/output)
      * \param resTrial predicted error for current trial field
      * \param lambda Anderson-Mixing mixing parameter 
      */
      void addPredictedError(DArray<double>& fieldTrial, 
                             DArray<double> const & resTrial, 
                             double lambda);

      /// Checks if the system has an initial guess
      bool hasInitialGuess();
     
      /** 
      * Compute the number of elements in field or residual.
      */
      int nElements();

      /*
      * Get the current state of the system.
      *
      * \param curr current field vector (output)
      */
      void getCurrent(DArray<double>& curr);

      /**
      * Solve MDE for current state of system.
      */
      void evaluate();

      /**
      * Gets the residual vector from system.
      *  
      * \param curr current residual vector (output)
      */
      void getResidual(DArray<double>& resid);

      /**
      * Update the system with a new trial field vector.
      *
      * \param newGuess trial field configuration
      */
      void update(DArray<double>& newGuess);

      /**
      * Output relevant system details to the iteration log file.
      */
      void outputToLog();

      // --- Private member functions specific to this implementation --- 
      
      cudaReal findAverage(cudaReal * const field, int n);

   };

} // namespace Pspg
} // namespace Pscf
#endif
