#ifndef PSPC_MC_MOVE_MANAGER_TPP
#define PSPC_MC_MOVE_MANAGER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMoveManager.h"
#include <pspc/System.h>
#include <pspc/mcmove/McMoveFactory.h>

#include <util/random/Random.h>
#include <util/global.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McMoveManager<D>::McMoveManager(System<D>& system)
   : Manager< McMove<D> >(),
     systemPtr_(&system),
     randomPtr_(&system.random())
   {  setClassName("McMoveManager"); }

   /*
   * Destructor.
   */
   template <int D>
   McMoveManager<D>::~McMoveManager()
   {}

   /*
   * Return a pointer to a new McMoveFactory object.
   */
   template <int D>
   Factory< McMove<D> >* McMoveManager<D>::newDefaultFactory() const
   {  return new McMoveFactory<D>(*systemPtr_); }

   /* 
   * Read instructions for creating objects from file.
   */
   template <int D>
   void McMoveManager<D>::readParameters(std::istream &in)
   {
      Manager< McMove<D> >::readParameters(in);

      // Allocate and store probabilities
      probabilities_.allocate(size());
      double  totalProbability = 0.0;
      int     iMove;
      for (iMove = 0; iMove < size(); ++iMove) {
         probabilities_[iMove] = (*this)[iMove].probability();
         totalProbability += probabilities_[iMove];
      }

      // Allocate and store and normalize probabilities
      for (iMove = 0; iMove < size(); ++iMove) {
         probabilities_[iMove] = probabilities_[iMove]/totalProbability;
         (*this)[iMove].setProbability(probabilities_[iMove]);
      }
   }

   /*
   * Load internal state from an archive.
   */
   template <int D>
   void McMoveManager<D>::loadParameters(Serializable::IArchive &ar)
   {
      Manager< McMove<D> >::loadParameters(ar);
      ar & probabilities_;
   }

   /*
   * Load internal state from an archive.
   */
   template <int D>
   void McMoveManager<D>::save(Serializable::OArchive &ar)
   {
      Manager< McMove<D> >::save(ar);
      ar & probabilities_;
   }

   /*
   * Initialize all moves just prior to a run.
   */
   template <int D>
   void McMoveManager<D>::setup()
   {
      for (int iMove = 0; iMove < size(); ++iMove) {
         (*this)[iMove].setup();
      }
   }

   /*
   * Choose a McMove at random.
   */
   template <int D>
   McMove<D>& McMoveManager<D>::chooseMove()
   {
      int iMove;
      iMove = randomPtr_->drawFrom(&probabilities_[0], size());
      return (*this)[iMove];
   }

   /*
   * Output statistics for every move.
   */
   template <int D>
   void McMoveManager<D>::output()
   {
      for (int i=0; i< size(); i++) {
         (*this)[i].output();
      }
   }

}
}
#endif
