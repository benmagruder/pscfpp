#ifndef RPG_THERMO_DERIVATIVE_ANALYZER_TPP
#define RPG_THERMO_DERIVATIVE_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ThermoDerivativeAnalyzer.h"

#include <rpg/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/fts/ramp/Ramp.h>

#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

#include <string>
#include <algorithm>

namespace Pscf {
namespace Rpg 
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ThermoDerivativeAnalyzer<D>::ThermoDerivativeAnalyzer(
                                Simulator<D>& simulator, System<D>& system) 
    : Analyzer<D>(),
      outputFileName_(""),
      hasAverage_(true),
      hasOutputFile_(false),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system()))
   { setClassName("ThermoDerivativeAnalyzer"); }

   /*
   * Destructor.
   */
   template <int D>
   ThermoDerivativeAnalyzer<D>::~ThermoDerivativeAnalyzer() 
   {}

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void ThermoDerivativeAnalyzer<D>::readParameters(std::istream& in) 
   {
      readInterval(in);
      read(in, "outputFileName", outputFileName_);
      readOptional(in, "hasAverage", hasAverage_);
      
      if (!outputFileName_.empty()) {
         hasOutputFile_ = true;
      }
      
      if (hasOutputFile_){
         system().fileMaster().openOutputFile(outputFileName_, outputFile_);
         outputFile_ << "    Variable       " << "Derivative" << "\n";
         
      }
   }
   
   /*
   * Setup before simulation loop.
   */ 
   template <int D>
   void ThermoDerivativeAnalyzer<D>::setup()
   {
      if (hasAverage_){
         accumulator_.clear();
      }
      
   }
   
   template <int D>
   void ThermoDerivativeAnalyzer<D>::sample(long iStep)
   {
      if (!isAtInterval(iStep)) return;
      
      double derivative;
      derivative = computeDerivative();
      if (hasAverage_){
         accumulator_.sample(derivative);
      }
      
      if (hasOutputFile_){
          UTIL_CHECK(outputFile_.is_open());
          outputFile_ << Dbl(variable());
          outputFile_ << Dbl(derivative);
          outputFile_<< "\n";
      }
   } 
   
   /*
   * Output results after a simulation is completed.
   */
   template <int D>
   void ThermoDerivativeAnalyzer<D>::output()
   {
      if (hasAverage_){
         Log::file() << std::endl;
         Log::file() << parameterType() << " : "<< std::endl;
         Log::file() << "Time average of the derivative is: "
                     << Dbl(accumulator_.average())
                     << " +- " << Dbl(accumulator_.blockingError(), 9, 2) 
                     << "\n";
      }
      
      // Close data file, if any
      if (outputFile_.is_open()) {
         outputFile_.close();
      }
      
      // Output average results to average file 
      if (hasAverage_){
         
         std::string fileName;
         std::string type;
         fileName = outputFileName_ +  ".ave";
         system().fileMaster().openOutputFile(fileName, outputFile_);
         
         double ave, err;
         ave = accumulator_.average();
         err = accumulator_.blockingError();
         outputFile_ << " " << std::left << std::setw(2) 
                      << parameterType() << "   ";
         outputFile_ << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
         outputFile_.close();
      }
      
   }

}
}
#endif