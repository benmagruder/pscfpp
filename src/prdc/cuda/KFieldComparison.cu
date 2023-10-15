/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "KFieldComparison.tpp"

namespace Pscf {
namespace Prdc {
namespace Cuda {

   template class KFieldComparison<1>;
   template class KFieldComparison<2>;
   template class KFieldComparison<3>;

} // namespace Pscf::Prdc::Cuda
} // namespace Pscf::Prdc
} // namespace Pscf
