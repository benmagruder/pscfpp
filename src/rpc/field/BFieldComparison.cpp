/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BFieldComparison.h"

namespace Pscf {
namespace Rpc {

   // Constructor
   BFieldComparison::BFieldComparison(int begin)
    : FieldComparison< DArray<double> > (begin)
   {};

}
}
