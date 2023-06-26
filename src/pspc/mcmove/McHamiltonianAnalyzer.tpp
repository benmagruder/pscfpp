#ifndef PSPC_MC_HAMILTONIAN_ANALYZER_TPP
#define PSPC_MC_HAMILTONIAN_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McHamiltonianAnalyzer.h"
#include "McSimulator.h"
#include <pspc/System.h>
#include <util/accumulators/Average.h>

namespace Pscf {
namespace Pspc 
{
   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McHamiltonianAnalyzer<D>::McHamiltonianAnalyzer(McSimulator<D>& mcSimulator, System<D>& system)
    : AverageListAnalyzer<D>(system),
      mcSimulatorPtr_(&mcSimulator),
      systemPtr_(&(mcSimulator.system())),
      idealId_(-1),
      fieldId_(-1),
      totalId_(-1)
   {  setClassName("McHamiltonianAnalyzer"); }

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void McHamiltonianAnalyzer<D>::readParameters(std::istream& in) 
   {
      AverageListAnalyzer<D>::readParameters(in);

      // Count number of values
      int id = 0;
      idealId_ = id;
      ++id; 
      fieldId_ = id;
      ++id;
      totalId_ = id;
      ++id; 
      AverageListAnalyzer<D>::initializeAccumulators(id);
 
      setName(idealId_, "ideal");
      setName(fieldId_, "field");
      setName(totalId_, "total");
   }

   /*
   * Output energy to file
   */
   template <int D>
   void McHamiltonianAnalyzer<D>::compute() 
   {
      double ideal = mcSimulator().mcIdealHamiltonian();
      //Log::file() << "ideal" << ideal<< "\n";
      // outputFile_ << Dbl(ideal, 20)
      setValue(idealId_, ideal);
      
      double field = mcSimulator().mcFieldHamiltonian();
      // outputFile_ << Dbl(field, 20)
      setValue(fieldId_, field);
      //Log::file() << "field" << field<< "\n";
      
      double total = mcSimulator().mcHamiltonian();
      // outputFile_ << Dbl(total, 20) << std::endl;
      setValue(totalId_, total);

   }

}
}
#endif