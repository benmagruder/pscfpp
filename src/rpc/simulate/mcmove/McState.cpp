#ifndef RPC_MC_STATE_CPP
#define RPC_MC_STATE_CPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McState.tpp"

namespace Pscf {
namespace Rpc {

   template struct McState<1>;
   template struct McState<2>;
   template struct McState<3>;

}
}
#endif
