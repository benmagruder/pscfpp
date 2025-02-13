#ifndef RPG_IMPOSED_FIELDS_GENERATOR_TEST_H
#define RPG_IMPOSED_FIELDS_GENERATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/scft/iterator/ImposedFieldsGenerator.h>
#include <rpg/System.h>

#include <prdc/cuda/RField.h>
#include <prdc/cuda/RFieldComparison.h>

#include <pscf/iterator/FieldGenerator.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

class ImposedFieldsGeneratorTest : public UnitTest
{

public:

   std::ofstream logFile_;

   void setUp()
   {  
      setVerbose(0); 
   }

   void tearDown()
   {
      if (logFile_.is_open()) {
         logFile_.close();
      }
   }

   void openLogFile(char const * filename)
   {
      openOutputFile(filename, logFile_);
      Log::setFile(logFile_);
   }

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      System<1> system;
      ImposedFieldsGenerator<1> ext(system);
   }

   void testReadParameters() // test ExtGenFilmBase::readParameters()
   {
      printMethod(TEST_FUNC);

      // Set up external field generator from file
      System<1> system;
      createSystem(system, "in/system1DGen");
      ImposedFieldsGenerator<1> gen(system);

      std::ifstream in;
      openInputFile("in/generator1", in);
      gen.readParam(in);
      in.close();

      // Check that the parameters that are publicly accessible were 
      // read correctly
      TEST_ASSERT(gen.type() == "film");
      TEST_ASSERT(gen.fieldGenerator1().type() == FieldGenerator::Mask);
      TEST_ASSERT(gen.fieldGenerator2().type() == FieldGenerator::External);
      
      DArray<int> ids;
      ids.allocate(1);
      ids[0] = 0;
      TEST_ASSERT(eq(gen.fieldGenerator2().getParameter("chi_bottom",ids), 5.0));
      TEST_ASSERT(eq(gen.fieldGenerator2().getParameter("chi_top",ids), 2.0));
      ids[0] = 1;
      TEST_ASSERT(eq(gen.fieldGenerator2().getParameter("chi_bottom",ids), 0.0));
      TEST_ASSERT(eq(gen.fieldGenerator2().getParameter("chi_top",ids), 10.0));
   }

   void testSolve1D() // solve a 1D system with an ImposedFieldsGenerator
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/GeneratorTestSolve1D.log");
      
      // Set up system with some data
      System<1> system;
      createSystem(system, "in/system1DGen");

      // Read initial guess
      system.readWBasis("in/wIn1D.bf");

      // Iterate to a solution
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0539149468e-01));

      // Check converged field is correct by comparing to ref files in in/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRef1D.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << rComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(rComparison.maxDiff() < 1.0E-5);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free energy error = " 
                   << (system.fHelmholtz() - 3.89135835701) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 11.7328996670) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.89135835701) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 11.7328996670) < 1e-4);
   }

   void testSolve2D() // solve a 2D system with an ImposedFieldsGenerator
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/GeneratorTestSolve2D.log");
      
      // Set up system with some data
      System<2> system;
      createSystem(system, "in/system2DGen");

      // Read initial guess
      system.readWBasis("in/wIn2D.bf");

      // Solve
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 7.99990525324e-01));
      
      // Check that lattice parameters are correct
      double aErr = system.domain().unitCell().parameter(0) - 1.63536608507;
      TEST_ASSERT(abs(aErr) < 1e-5);
      TEST_ASSERT(eq(system.domain().unitCell().parameter(1), 2.0));

      // Check converged field is correct by comparing to reference
      UnitCell<2> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<2> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRef2D.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<2> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
     
      double epsilon = 1.0E-4; 
      double diff = rComparison.maxDiff();

      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "Max field error = " << diff << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free energy error = " 
                   << (system.fHelmholtz() - 3.91021919092) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 12.4986763317) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.91021919092) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 12.4986763317) < 1e-4);
   }

   void testSweep() // test sweep along chiBottom and lattice parameter
   {
      // NOTE: this also tests that the ParameterModifier methods work
      printMethod(TEST_FUNC);
      
      openLogFile("out/GeneratorTestSweep.log");
      
      // Set up system
      System<1> system;
      createSystem(system, "in/system1DGen");

      // Read initial guess
      system.readWBasis("in/wIn1D.bf");

      // Run the sweep function
      system.sweep();

      // Check converged field is correct by comparing to reference
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRefSweep.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      double diff = rComparison.maxDiff();

      double epsilon = 1.0E-5; 
      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "Max field error = " << diff << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free Energy error = " 
                   << (system.fHelmholtz() - 3.86654196623) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 11.4844864688) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.86654196623) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 11.4844864688) < 1e-4);
   }

   void testSolveWithFBulk() // solve a 1D system w/ flexible film thickness
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/GeneratorTestSolveWithFBulk.log");
      
      // Set up system with some data
      System<1> system;
      createSystem(system, "in/system1DGenFBulk");

      // Read initial guess
      system.readWBasis("in/wIn1DFBulk.bf");

      // Iterate to a solution
      system.iterate();
      
      // Check that the right film thickness was found
      double paramErr = system.domain().unitCell().parameter(0) - 2.05584596449;
      if (verbose() > 0) {
         std::cout << "\nFilm thickness error = " << paramErr << "\n";
      }
      TEST_ASSERT(abs(paramErr) < 1e-5);
      TEST_ASSERT(abs(system.mask().phiTot() - 0.80542424387) < 1e-5);

      // Check converged field is correct by comparing to ref files in in/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRef1DFBulk.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      double diff = rComparison.maxDiff();
      if (verbose() > 0) {
         std::cout << "Max field error = " << diff << "\n";
      }
      TEST_ASSERT(diff < 1.0E-5);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free Energy error = " 
                   << (system.fHelmholtz() - 3.89123697966) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 11.8931100854) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.89123697966) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 11.8931100854) < 1e-4);
   }

   void testSolve1DGrid() // solve a 1D system with an ImposedFieldsGenerator
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/GeneratorTestSolve1DGrid.log");
      
      // Set up system with some data
      System<1> system;
      createSystem(system, "in/system1DGenGrid");

      // Read initial guess
      system.readWBasis("in/wIn1D.bf");

      // Iterate to a solution
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0539149468e-01));

      // Check converged field is correct by comparing to ref files in in/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRef1D.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << rComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(rComparison.maxDiff() < 1.0E-5);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free Energy error = " 
                   << (system.fHelmholtz() - 3.89135835701) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 11.7328996670) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.89135835701) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 11.7328996670) < 1e-4);
   }

   void testSolve2DGrid() // solve a 2D system with an ImposedFieldsGenerator
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/GeneratorTestSolve2DGrid.log");
      
      // Set up system with some data
      System<2> system;
      createSystem(system, "in/system2DGenGrid");

      // Read initial guess
      system.readWBasis("in/wIn2D.bf");

      // Solve
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 7.99990525324e-01));
      
      // Check that lattice parameters are correct
      double aErr = system.domain().unitCell().parameter(0) - 1.63536608507;
      TEST_ASSERT(abs(aErr) < 1e-5);
      TEST_ASSERT(eq(system.domain().unitCell().parameter(1), 2.0));

      // Check converged field is correct by comparing to reference
      UnitCell<2> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<2> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRef2D.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<2> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
     
      double epsilon = 1.0E-4; 
      double diff = rComparison.maxDiff();

      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "Max field error = " << diff << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free Energy error = " 
                   << (system.fHelmholtz() - 3.91022221196) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 12.4995705392) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.91022221196) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 12.4995705392) < 1e-4);
   }

   void testSweepGrid() // test sweep along chiBottom and lattice parameter
   {
      // NOTE: this also tests that the ParameterModifier methods work
      printMethod(TEST_FUNC);
      
      openLogFile("out/GeneratorTestSweepGrid.log");
      
      // Set up system
      System<1> system;
      createSystem(system, "in/system1DGenGrid");

      // Read initial guess
      system.readWBasis("in/wIn1D.bf");

      // Run the sweep function
      system.sweep();

      // Check converged field is correct by comparing to reference
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRefSweep.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      double diff = rComparison.maxDiff();

      double epsilon = 1.0E-5; 
      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "Max field error = " << diff << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free Energy error = " 
                   << (system.fHelmholtz() - 3.86656498107) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 11.4994801703) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.86656498107) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 11.4994801703) < 1e-4);
   }

   void testSolveWithFBulkGrid() // solve a 1D system w flexible film thickness
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/GeneratorTestSolveWithFBulkGrid.log");
      
      // Set up system with some data
      System<1> system;
      createSystem(system, "in/system1DGenFBulkGrid");

      // Read initial guess
      system.readWBasis("in/wIn1DFBulk.bf");

      // Iterate to a solution
      system.iterate();
      
      // Check that the right film thickness was found
      double paramErr = system.domain().unitCell().parameter(0) - 2.05584596449;
      if (verbose() > 0) {
         std::cout << "\nFilm thickness error = " << paramErr << "\n";
      }
      TEST_ASSERT(abs(paramErr) < 1e-5);
      TEST_ASSERT(abs(system.mask().phiTot() - 0.80542424387) < 1e-5);

      // Check converged field is correct by comparing to ref files in in/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRef1DFBulk.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      double diff = rComparison.maxDiff();
      if (verbose() > 0) {
         std::cout << "Max field error = " << diff << "\n";
      }
      TEST_ASSERT(diff < 1.0E-5);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free Energy error = " 
                   << (system.fHelmholtz() - 3.89123697966) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 11.8931100854) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.89123697966) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 11.8931100854) < 1e-4);
   }

   // Read parameter file to create a System object
   template <int D>
   void createSystem(System<D>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

};

TEST_BEGIN(ImposedFieldsGeneratorTest)
TEST_ADD(ImposedFieldsGeneratorTest, testConstructor)
TEST_ADD(ImposedFieldsGeneratorTest, testReadParameters)
TEST_ADD(ImposedFieldsGeneratorTest, testSolve1D)
TEST_ADD(ImposedFieldsGeneratorTest, testSolve2D)
TEST_ADD(ImposedFieldsGeneratorTest, testSweep)
TEST_ADD(ImposedFieldsGeneratorTest, testSolveWithFBulk)
TEST_ADD(ImposedFieldsGeneratorTest, testSolve1DGrid)
TEST_ADD(ImposedFieldsGeneratorTest, testSolve2DGrid)
TEST_ADD(ImposedFieldsGeneratorTest, testSweepGrid)
TEST_ADD(ImposedFieldsGeneratorTest, testSolveWithFBulkGrid)
TEST_END(ExtGenFilmTest)

#endif
