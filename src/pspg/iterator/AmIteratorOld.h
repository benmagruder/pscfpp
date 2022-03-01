#ifndef PSPG_AM_ITERATOR_OLD_H
#define PSPG_AM_ITERATOR_OLD_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspg/iterator/Iterator.h> // base class
#include <pspg/solvers/Mixture.h>
#include <pscf/math/LuSolver.h>
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>
#include <util/containers/RingBuffer.h>
#include <pspg/iterator/HistMatOld.h>
#include <pspg/field/RDField.h>
#include <util/containers/FArray.h>

// Always defined. Locks out legacy code that does nothing. 
// Eventually remove all legacy code that references memory in cpu
//#define GPU_OUTER

// Define this if you want to use scft
//#define GPU_SCFT

namespace Pscf {
namespace Pspg {

   using namespace Util;

   /**
   * Anderson mixing iterator for the pseudo spectral method
   *
   * \ingroup Pspg_Iterator_Module
   */
   template <int D>
   class AmIteratorOld : public Iterator<D>
   {
   public:

      typedef RDField<D> WField;
      typedef RDField<D> CField;

      /**
      * Constructor
      *
      * \param system pointer to a system object
      */
      AmIteratorOld(System<D>& system);

      /**
      * Destructor
      */
      ~AmIteratorOld();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      void readParameters(std::istream& in);

      /**
      * Allocate all arrays
      *
      */
      void setup();

      /**
      * Iterate to a solution
      *
      */
      int solve();

      /**
      * Getter for epsilon
      */
      double epsilon();

      /**
      * Getter for the maximum number of field histories to
      * convolute into a new field
      */
      int maxHist();

      /**
      * Getter for the maximum number of iteration before convergence
      */
      int maxItr();

      /**
      * Compute the deviation of wFields from a mean field solution
      */
      void computeDeviation();

      /**
      * Compute the error from deviations of wFields and compare with epsilon_
      * \return true for error < epsilon and false for error >= epsilon
      */
      bool isConverged();

      /**
      * Determine the coefficients that would minimize invertMatrix_ Umn
      * 
      */
      int minimizeCoeff(int itr);

      /**
      * Rebuild wFields for the next iteration from minimized coefficients
      */
      void buildOmega(int itr);

   private:

      ///error tolerance
      double epsilon_;

      /// Type of error checked for convergence.
      /// Either maxResid or normResid.
      std::string errorType_;

      /// Flexible unit cell (1) or rigid cell (0), default value = 0
      bool isFlexible_;

      ///free parameter for minimization
      double lambda_;

      ///number of histories to convolute into a new solution [0,maxHist_]
      int nHist_;

      //maximum number of histories to convolute into a new solution
      //AKA size of matrix
      int maxHist_;

      ///number of maximum iteration to perform
      int maxItr_;

      /// holds histories of deviation for each monomer
      /// 1st index = history, 2nd index = monomer, 3rd index = ngrid
      // The ringbuffer used is now slightly modified to return by reference

      RingBuffer< DArray < RDField<D> > > d_resHists_;
      RingBuffer< DArray < RDField<D> > > d_omHists_;

      /// holds histories of deviation for each cell parameter
      /// 1st index = history, 2nd index = cell parameter
      // The ringbuffer used is now slightly modified to return by reference
      RingBuffer< FArray <double, 6> > devCpHists_;
      RingBuffer< FSArray<double, 6> > CpHists_;

      FSArray<double, 6> cellParameters_;

      /// Umn, matrix to be minimized
      DMatrix<double> invertMatrix_;

      /// Cn, coefficient to convolute previous histories with
      DArray<double> coeffs_;

      DArray<double> vM_;

      /// bigW, blended omega fields
      DArray<RDField<D> > wArrays_;

      /// bigWcP, blended parameter
      FArray <double, 6> wCpArrays_;

      /// bigDCp, blened deviation parameter. new wParameter = bigWCp + lambda * bigDCp
      FArray <double, 6> dCpArrays_;

      /// bigD, blened deviation fields. new wFields = bigW + lambda * bigD
      DArray<RDField<D> > dArrays_;

      DArray<RDField<D> > tempDev;

      HistMatOld <cudaReal> histMat_;

      cudaReal* d_temp_;
      cudaReal* temp_;

      using Iterator<D>::setClassName;
      using Iterator<D>::sys_;
      using ParamComposite::read;
      using ParamComposite::readOptional;

      //friend:
      //for testing purposes

   };

   template<int D>
   inline double AmIteratorOld<D>::epsilon()
   {
      return epsilon_;
   }

   template<int D>
   inline int AmIteratorOld<D>::maxHist()
   {
      return maxHist_;
   }

   template<int D>
   inline int AmIteratorOld<D>::maxItr()
   {
      return maxItr_;
   }

   #ifndef PSPG_AM_ITERATOR_OLD_TPP
   // Suppress implicit instantiation
   extern template class AmIteratorOld<1>;
   extern template class AmIteratorOld<2>;
   extern template class AmIteratorOld<3>;
   #endif

}
}
#endif
