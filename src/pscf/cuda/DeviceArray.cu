/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DeviceArray.tpp"

namespace Pscf {

   template class DeviceArray<cudaReal>;
   template class DeviceArray<cudaComplex>;
   template class DeviceArray<int>;
   template class DeviceArray<bool>;

}