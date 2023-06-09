#ifndef PSPC_MC_SIMULATOR_H
#define PSPC_MC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McState.h"                     // member
#include "McMoveManager.h"               // member

#include <util/param/Manager.h>          // base class template
#include <util/param/ParamComposite.h>   // base class
#include <util/random/Random.h>          // member
#include <util/containers/DArray.h>      // member template
#include <util/containers/DMatrix.h>     // member template

namespace Pscf {
namespace Pspc {

   using namespace Util;

   template <int D> class System;
   template <int D> class McMove;

   /**
   * Resources for a Monte-Carlo simulation of system.
   */
   template <int D>
   class McSimulator : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System
      */
      McSimulator(System<D>& system);

      /**
      * Destructor.
      */
      ~McSimulator();

      /**
      * Read parameters for a MC simulation.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /**
      * Perform a field theoretic Monte-Carlo simulation.
      *
      * Perform a field theoretic Monte-Carlo simulation using the 
      * partial saddle-point approximation. 
      * 
      * \param nStep  number of Monte-Carlo steps
      */
      void simulate(int nStep);

      /**
      * Compute the Hamiltonian used in Monte-Carlo simulations.
      */
      void computeMcHamiltonian();

      /**
      * Get the Hamiltonian used in Monte-Carlo simulations.
      */
      double mcHamiltonian() const;

      /**
      * Save a copy of the Monte-Carlo state.
      *
      * This function and restoreMcState() are intended for use in 
      * the implementation of field theoretic Monte Carlo moves. This
      * function stores the current w fields and the corresponding 
      * Hamiltonian value.  This is normally the first step of a MC 
      * move, prior to an attempted modification of the fields stored
      * in the system w field container.
      */
      void saveMcState();

      /**
      * Restore the saved copy of the Monte-Carlo state.
      *
      * This function  and saveMcState() are intended to be used 
      * together in the implementation of Monte-Carlo moves. If an 
      * attempted move is rejected, restoreMcState() is called to
      * restore the fields ahd Hamiltonian value that were saved 
      * by a previous call to the function saveMcState().
      */
      void restoreMcState();

      /**
      * Clear the saved copy of the Monte-Carlo state.
      *
      * This function, restoreMcState(), and saveMcState() are intended to be used 
      * together in the implementation of Monte-Carlo moves. If an 
      * attempted move is accepted, clearMcState() is called to
      * clear mcState_.hasData
      */
      void clearMcState();
      
      /**
      * Get an array of the eigenvalues of the projected chi matrix.
      *
      * The projected chi matrix is given by the matrix product P*chi*P,
      * where P is the symmetric projection matrix that projects onto the
      * subspace perpendicular to the vector (1,1,...,1). The projected
      * chi matrix is singular, and always has one zero eigenvalue, with 
      * eigenvector (1,1, ... ,1). By convention, this zero value is the 
      * last eigenvalue, with index nMonomer -  1, is zero.
      */
      DArray<double> const & chiEvals() const
      {  return chiEvals_; }

      /**
      * Get an single eigenvalue of the projected chi matrix.
      *
      * \int i index of eigenvalue (0, ... , nMonomer - 1)
      */
      double chiEval(int i ) const
      {  return chiEvals_[i]; }

      /**
      * Get a matrix of eigenvectors of the projected chi matrix.
      *
      * The first index of the matrix indexes the eigenvector, while 
      * the second index indexes vector components. All eigenvectors
      * have a Euclidean norm of unity. The sign of each vector is 
      * chosen so as to make the first (0) component of each vector
      * positive.
      */
      DMatrix<double> const & chiEvecs() const
      {  return chiEvecs_; }

      /**
      * Compute and store the eigenvector components of the current w fields.
      */
      void computeWC();

      /**
      * Get an eigenvector component of the w fields.
      *
      * Each component is a point-wise projection of the w fields onto
      * a corresponding eigenvector of the projected chi matrix.
      *
      * \int i eigenvector / eigenvalue index
      */
      RField<D> const & wc(int i) const
      {   return wc_[i]; }

      /**
      * Get parent system by reference.
      */
      System<D>& system();

      /**
      * Get McMoveManger
      */
      McMoveManager<D>& mcMoveManager();
      
      /**
      * Get random number generator by reference.
      */
      Random& random();

   private:
      
      // Private data members

      /**
      * Manger for Monte Carlo Move.
      */
      McMoveManager<D> mcMoveManager_;  
      
      /**
      * Random number generator
      */
      Random random_;
      
      /**
      * State saved during MC simulation.
      */
      mutable McState<D> mcState_;

      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_;  

      /**
      * Eigenvector components of w on a real space grid.
      *
      * Each field component corresponds to a point-wise projection of w 
      * onto an eigenvector of the projected chi matrix.
      */
      DArray< RField<D> > wc_;

      /**
      * Projected chi matrix
      *
      * Projected matrix chiP_ = P*chi*P, where P is projection matrix that 
      * projects onto the subspace orthogonal to the vector e = [1, ... , 1].
      */
      DMatrix<double> chiP_;

      /**
      * Eigenvectors of the projected chi matrix.
      */
      DMatrix<double> chiEvecs_;

      /**
      * Eigenvalues of the projected chi matrix.
      */
      DArray<double>  chiEvals_;

      /**
      * Monte-Carlo System Hamiltonian (extensive value).
      */
      double mcHamiltonian_;

      /**
      * Has the MC Hamiltonian been computed for the current w and c fields?
      */ 
      bool hasMcHamiltonian_;

      // Private member functions
       
      /**
      * Perform eigenvalue analysis of projected chi matrix.
      */
      void analyzeChi();

      /**
      * Called at the beginning of the simulation member function.
      */
      void setup();

   };

   // Inline functions

   // Get the Monte-Carlo move manager.
   template <int D>
   inline McMoveManager<D>& McSimulator<D>::mcMoveManager()
   {  return mcMoveManager_; }

   // Get the random number generator.
   template <int D>
   inline Random& McSimulator<D>::random()
   {  return random_; }
   
   // Get the parent system.
   template <int D>
   inline System<D>& McSimulator<D>::system()
   {  return *systemPtr_; }
   
   // Get the precomputed MC Hamiltonian
   template <int D>
   inline double McSimulator<D>::mcHamiltonian() const
   {  
      UTIL_CHECK(hasMcHamiltonian_);
      return mcHamiltonian_; 
   }

   #ifndef PSPC_MC_SIMULATOR_TPP
   // Suppress implicit instantiation
   extern template class McSimulator<1>;
   extern template class McSimulator<2>;
   extern template class McSimulator<3>;
   #endif

}
}
#endif
