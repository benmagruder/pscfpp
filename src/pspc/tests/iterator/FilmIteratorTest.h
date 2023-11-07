#ifndef PSPC_FILM_ITERATOR_TEST_H
#define PSPC_FILM_ITERATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/iterator/FilmIterator.h>
#include <pspc/iterator/AmIteratorBasis.h>
#include <pspc/field/FieldIo.h>
#include <pspc/System.h>

#include <prdc/cpu/RFieldComparison.h>
#include <prdc/crystal/BFieldComparison.h>
#include <prdc/crystal/UnitCell.h>

#include <util/misc/Exception.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;
using namespace Pscf::Pspc;

class FilmIteratorTest : public UnitTest
{

public:

   std::ofstream logFile_;

   void setUp()
   {  setVerbose(0); }

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
      System<3> system;
      FilmIterator<3, AmIteratorBasis<3> > iterator(system);
   }

   void testReadParameters() // test FilmIterator::readParameters()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestReadParameters.log");

      // Set up system with some data
      System<2> system;
      FilmIteratorTest::setUpSystem(system, "in/film/system2D");

      // Set up iterator from file
      FilmIterator<2, AmIteratorBasis<2> > iterator(system);
      FilmIteratorTest::setUpFilmIterator(iterator, "in/film/film2D");

      // Check that everything was read in correctly
      TEST_ASSERT(eq(iterator.normalVecId(),1));
      TEST_ASSERT(eq(iterator.interfaceThickness(),0.2));
      TEST_ASSERT(eq(iterator.wallThickness(),0.4));
      TEST_ASSERT(eq(iterator.chiBottom(0),3.0));
      TEST_ASSERT(eq(iterator.chiBottom(1),0.0));
      TEST_ASSERT(eq(iterator.chiTop(0),0.0));
      TEST_ASSERT(eq(iterator.chiTop(1),4.0));
      TEST_ASSERT(iterator.isFlexible());
   }

   void testGenerateWallFields() // testFilmIterator::generateWallFields()
   {      

      printMethod(TEST_FUNC);
      openLogFile("out/filmTestGenerateWallFields.log");
      
      // Set up system with some data
      System<1> system;
      FilmIteratorTest::setUpSystem(system, "in/film/system1D");

      // Set up iterator from file
      FilmIterator<1, AmIteratorBasis<1> > iterator(system);
      FilmIteratorTest::setUpFilmIterator(iterator, "in/film/film1D");

      system.readWBasis("in/film/w_1D_ref.bf");

      // Allocate mask and external field containers
      system.mask().allocate(system.basis().nBasis(), 
                             system.mesh().dimensions());
      system.h().allocateRGrid(system.mesh().dimensions());
      system.h().allocateBasis(system.basis().nBasis());

      // Run the generateWallFields function
      iterator.generateWallFields();

      // Check that the homogeneous components of the mask
      // and the blocks were adjusted correctly
      TEST_ASSERT(eq(system.mask().phiTot(),8.0951532073e-01));

      // output mask field files for reference
      system.fieldIo().writeFieldBasis("out/mask_1D.bf", system.mask().basis(),
                                       system.unitCell());
      system.fieldIo().writeFieldRGrid("out/mask_1D.rf", system.mask().rgrid(),
                                       system.unitCell());
      
      // output external field for reference
      system.fieldIo().writeFieldsBasis("out/h_1D.bf", system.h().basis(),
                                       system.unitCell());

      // Check that the mask field files were generated correctly by 
      // comparing them to the reference files in in/film
      UnitCell<1> unitCell; // UnitCell object to pass into FieldIo functions
      DArray<double> cFieldsCheck; // Copy of reference field
      system.fieldIo().readFieldBasis("in/film/mask_1D_ref.bf", 
                                       cFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.mask().basis(), cFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << bComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-7);

      RField<1> cRGridCheck; // Array to store reference field
      system.fieldIo().readFieldRGrid("in/film/mask_1D_ref.rf", 
                                      cRGridCheck, unitCell);
      RField<1> cRGridFromIterator;
      cRGridFromIterator.allocate(system.domain().mesh().dimensions());

      // Put iterator cField inside a DArray so it can be passed into 
      // convertBasisToRGrid
      DArray<double> cFieldFromIterator;   
      cFieldFromIterator = system.mask().basis(); 

      system.fieldIo().convertBasisToRGrid(cFieldFromIterator,
                                           cRGridFromIterator);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(cRGridFromIterator, cRGridCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << rComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(rComparison.maxDiff() < 1.0E-7);
   }

   void testCheckSpaceGroup1DA() // test FilmIterator::checkSpaceGroup
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckSpaceGroup1DA.log");

      // Set up 1D system with a correct space group and check it
      System<1> system1;
      FilmIteratorTest::setUpSystem(system1, "in/film/system1D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      double parameter = 2.9;
      parameters.append(parameter);
      system1.setUnitCell(UnitCell<1>::Lamellar, parameters);
      //system1.readWBasis("in/film/w_ref.bf");

      FilmIterator<1, AmIteratorBasis<1> > iterator1(system1);
      FilmIteratorTest::setUpFilmIterator(iterator1, "in/film/film1D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator1,false));
   }

   void testCheckSpaceGroup1DB() // test FilmIterator::checkSpaceGroup
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckSpaceGroup1DB.log");

      // Set up 1D system with an incorrect space group and check it
      System<1> system2;
      FilmIteratorTest::setUpSystem(system2, "in/film/system_bad_1D");
      FilmIterator<1, AmIteratorBasis<1> > iterator2(system2);

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      //double parameter = 2.2;
      parameters.append(2.2);
      system2.setUnitCell(UnitCell<1>::Lamellar, parameters);

      FilmIteratorTest::setUpFilmIterator(iterator2, "in/film/film1D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator2,true));
   }

   void testCheckSpaceGroup2D() // test FilmIterator::checkSpaceGroup
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckSpaceGroup2D.log");

      // Set up 2D system with an incorrect space group and check it
      System<2> system3;
      FilmIteratorTest::setUpSystem(system3, "in/film/system_bad_2D_2");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(2.0);
      system3.setUnitCell(UnitCell<2>::Rectangular, parameters);

      FilmIterator<2, AmIteratorBasis<2> > iterator3(system3);

      FilmIteratorTest::setUpFilmIterator(iterator3, "in/film/film2D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator3,true));
   }

   void testCheckSpaceGroup3DA() 
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckSpaceGroup3DA.log");

      // Set up 3D system with a correct space group and check it
      System<3> system4;
      FilmIteratorTest::setUpSystem(system4, "in/film/system3D");
      FilmIterator<3, AmIteratorBasis<3> > iterator4(system4);
      FilmIteratorTest::setUpFilmIterator(iterator4, "in/film/film3D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(4.2);
      system4.setUnitCell(UnitCell<3>::Tetragonal, parameters);

      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator4,false));
      TEST_ASSERT(iterator4.isFlexible()); // check that isFlexible works
   }

   void testCheckSpaceGroup3DB() 
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckSpaceGroup3DB.log");

      // Set up 3D system with an incorrect space group and check it
      System<3> system5;
      FilmIteratorTest::setUpSystem(system5, "in/film/system_bad_3D_1");
      FilmIterator<3, AmIteratorBasis<3> > iterator5(system5);
      FilmIteratorTest::setUpFilmIterator(iterator5, "in/film/film3D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(4.2);
      system5.setUnitCell(UnitCell<3>::Tetragonal, parameters);

      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator5,true));

      // Set up another 3D system with an incorrect space group and check it
      System<3> system6;
      FilmIteratorTest::setUpSystem(system6, "in/film/system_bad_3D_2");
      FilmIterator<3, AmIteratorBasis<3> > iterator6(system6);
      FilmIteratorTest::setUpFilmIterator(iterator6, "in/film/film3D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator6,true));
   }

   void testCheckLatticeVectors() // test FilmIterator::checkLatticeVectors()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckLatticeVectors.log");
      
      // Set up 2D system with incorrect lattice vectors and check it
      System<2> system1;
      FilmIteratorTest::setUpSystem(system1, "in/film/system_bad_2D_1");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(2.0);
      parameters.append(1.0);
      system1.setUnitCell(UnitCell<2>::Oblique, parameters);

      FilmIterator<2, AmIteratorBasis<2> > iterator1(system1);
      FilmIteratorTest::setUpFilmIterator(iterator1, "in/film/film2D");
      try {
         iterator1.checkLatticeVectors();
         // If above does not throw an error, then it failed this test
         TEST_ASSERT(1 == 2);
      } catch (Exception& e) {
         Log::file() << "EXCEPTION CAUGHT, expected behavior occurred" 
                     << std::endl;
      }

      // Set up 3D system with correct lattice vectors and check it
      System<3> system2;
      FilmIteratorTest::setUpSystem(system2, "in/film/system_bad_3D_1");
      parameters.clear();
      parameters.append(2.0);
      parameters.append(4.2);
      system2.setUnitCell(UnitCell<3>::Tetragonal, parameters);
      FilmIterator<3, AmIteratorBasis<3> > iterator2(system2);
      FilmIteratorTest::setUpFilmIterator(iterator2, "in/film/film3D");
      iterator2.checkLatticeVectors(); // this should not throw an error

      // Set up 3D system with incorrect lattice vectors and check it
      System<3> system3;
      FilmIteratorTest::setUpSystem(system3, "in/film/system_bad_3D_2");
      parameters[1] = 2.0;
      parameters.append(2.0); 
      parameters.append(1.0);
      system3.setUnitCell(UnitCell<3>::Monoclinic, parameters);
      FilmIterator<3, AmIteratorBasis<3> > iterator3(system3);
      FilmIteratorTest::setUpFilmIterator(iterator3, "in/film/film3D");
      try {
         iterator3.checkLatticeVectors();
         // If above doesn't throw an error, then it failed this test
         TEST_ASSERT(1 == 2);
      } catch (Exception& e) {
         Log::file() << "EXCEPTION CAUGHT, expected behavior occurred" 
                     << std::endl;
      }
   }
   
   void testFlexibleParams() // test FilmIterator::flexibleParams
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestFlexibleParams.log");
      
      // Set up 1D system and make sure flexibleParams is empty
      System<1> system1;
      FilmIteratorTest::setUpSystem(system1, "in/film/system1D");
      FilmIterator<1, AmIteratorBasis<1> > iterator1(system1);
      FilmIteratorTest::setUpFilmIterator(iterator1, "in/film/film1D");
      TEST_ASSERT(iterator1.nFlexibleParams() == 0);

      // Set up 2D system and make sure flexibleParams is correct
      System<2> system2;
      FilmIteratorTest::setUpSystem(system2, "in/film/system2D");
      FilmIterator<2, AmIteratorBasis<2> > iterator2(system2);
      FilmIteratorTest::setUpFilmIterator(iterator2, "in/film/film2D");
      TEST_ASSERT(iterator2.nFlexibleParams() == 1);
      TEST_ASSERT(iterator2.flexibleParams()[0]);

      // Set up 3D tetragonal system, validate flexibleParams 
      System<3> system3;
      FilmIteratorTest::setUpSystem(system3, "in/film/system3D");
      FilmIterator<3, AmIteratorBasis<3> > iterator3(system3);
      FilmIteratorTest::setUpFilmIterator(iterator3, "in/film/film3D");
      TEST_ASSERT(iterator3.nFlexibleParams() == 1);
      TEST_ASSERT(iterator3.flexibleParams()[0]);

      // Set up 3D monoclinic system (monoclinic), validate flexibleParams 
      System<3> system4;
      FilmIteratorTest::setUpSystem(system4, "in/film/system_bad_3D_2");
      FilmIterator<3, AmIteratorBasis<3> > iterator4(system4);
      FilmIteratorTest::setUpFilmIterator(iterator4, "in/film/film2D");
      // Using film2D here because it has normalVecId=1 which 
      // we want for this example
      TEST_ASSERT(iterator4.nFlexibleParams() == 3);
      TEST_ASSERT(iterator4.flexibleParams()[0]);
      TEST_ASSERT(iterator4.flexibleParams()[2]);
      TEST_ASSERT(iterator4.flexibleParams()[3]);
   }

   void testReadFlexibleParams() // test manual entry of flexibleParams
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/filmTestReadFlexibleParams.log");
      
      // Set up system
      System<2> system;
      FilmIteratorTest::setUpSystem(system, "in/film/system_bad_2D_1");

      // Check flexibleParams array
      TEST_ASSERT(system.iterator().nFlexibleParams() == 2);
      TEST_ASSERT(system.iterator().flexibleParams()[0]);
      TEST_ASSERT(system.iterator().flexibleParams()[2]);

      // Set up another system
      System<3> system2;
      FilmIteratorTest::setUpSystem(system2, "in/film/system_bad_3D_2");

      // Check flexibleParams array
      TEST_ASSERT(system2.iterator().nFlexibleParams() == 3);
      TEST_ASSERT(system2.iterator().flexibleParams()[0]);
      TEST_ASSERT(system2.iterator().flexibleParams()[2]);
      TEST_ASSERT(system2.iterator().flexibleParams()[3]);
   }

   void testSolve1D() // test FilmIterator::solve
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/filmTestSolve1D.log");
      
      // Set up system with some data
      System<1> system;
      FilmIteratorTest::setUpSystem(system, "in/film/system1D");

      // Read initial guess
      system.readWBasis("in/film/w_1D_in.bf");

      // Set up iterator from file
      FilmIterator<1, AmIteratorBasis<1> > iterator(system);
      FilmIteratorTest::setUpFilmIterator(iterator, "in/film/film1D");

      // Run the solve function
      iterator.solve();
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0951532073e-01));

      // Check converged field is correct by comparing to files in in/film
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.fieldIo().readFieldsBasis("in/film/w_1D_ref.bf", 
                                       wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << bComparison.maxDiff() << "\n";
      }
      system.fieldIo().writeFieldsBasis("out/w_1D.bf", system.w().basis(), 
                                        system.unitCell());
      system.fieldIo().writeFieldsRGrid("out/w_1D.rf", system.w().rgrid(), 
                                        system.unitCell());
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-5);
   }

   void testSolve2D() // test FilmIterator::solve
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/filmTestSolve2D.log");
      
      // Set up system with some data
      System<2> system;
      FilmIteratorTest::setUpSystem(system, "in/film/system2D");

      // Read initial guess
      system.readWBasis("in/film/w_2D_in.bf");

      // Solve
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 8.7096155661e-01));
      
      // Check that lattice parameters are correct
      TEST_ASSERT((system.unitCell().parameter(0) - 1.6418585139) < 1.0e-6);
      TEST_ASSERT(eq(system.unitCell().parameter(1), 3.1));

      // Check converged field is correct by comparing to reference
      UnitCell<2> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.fieldIo().readFieldsBasis("in/film/w_2D_ref.bf", 
                                       wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
      system.writeWBasis("out/w_2D.bf");
     
      double epsilon = 1.0E-5; 
      double diff = bComparison.maxDiff();

      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "diff    = " << diff << "\n";
         std::cout << "epsilon = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
   }

   void testSweep() // test sweep along chiBottom and lattice parameter
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/filmTestSweep.log");
      
      // Set up system
      System<1> system;
      FilmIteratorTest::setUpSystem(system, "in/film/system1D");

      // Read initial guess
      system.readWBasis("out/w_1D.bf");

      // Run the sweep function
      system.sweep();

      // Check converged field is correct by comparing to reference
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.fieldIo().readFieldsBasis("in/film/w_1D_ref_sweep.bf", 
                                       wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
      double diff = bComparison.maxDiff();

      double epsilon = 1.0E-5; 
      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "diff    = " << diff << "\n";
         std::cout << "epsilon = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
   }

   void testFreeEnergy() // test System::computeFreeEnergy with mask/h fields
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/filmTestFreeEnergy.log");
      
      // Set up system
      System<1> system;
      FilmIteratorTest::setUpSystem(system, "in/film/system1D");

      // Read a converged solution as initial guess
      system.readWBasis("out/w_1D.bf");

      // Set up iterator
      FilmIterator<1, AmIteratorBasis<1> > iterator(system);
      FilmIteratorTest::setUpFilmIterator(iterator, "in/film/film1D");

      // Solve (should only take a few iterations)
      iterator.solve();

      // Compute free energy
      system.computeFreeEnergy();
      system.writeThermo(Log::file());

      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.fHelmholtz() - 3.87784944222) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.pressure() + 12.1117881919) << "\n";
      }
      TEST_ASSERT(system.fHelmholtz() - 3.87784944222 < 1e-6);
      TEST_ASSERT(system.pressure() + 12.1117881919 < 1e-5);
   }

   void testMaskAndH() // test manual entry of mask and h fields
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/filmTestMaskAndH.log");
      
      // Set up system
      System<1> system;
      FilmIteratorTest::setUpSystem(system, "in/film/system1D_noFilm");

      // Read the same initial guess as testSolve
      system.readWBasis("in/film/w_1D_in.bf");

      // Read in the mask and external fields from file
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      unitCell = system.unitCell();
      system.mask().setFieldIo(system.fieldIo());
      system.mask().allocate(system.basis().nBasis(), 
                             system.mesh().dimensions());
      system.mask().readBasis("out/mask_1D.bf", unitCell);
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0951532073e-01));

      system.h().setFieldIo(system.fieldIo());
      system.h().allocateBasis(system.basis().nBasis());
      system.h().allocateRGrid(system.mesh().dimensions());
      system.h().readBasis("out/h_1D.bf", unitCell);

      // Run the solve function
      system.iterate();

      // Check converged field is correct by comparing to files in in/film
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.fieldIo().readFieldsBasis("in/film/w_1D_ref.bf", wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);

      double diff = bComparison.maxDiff();
      double epsilon = 1.0E-5; 
      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "diff    = " << diff << "\n";
         std::cout << "epsilon = " << epsilon << "\n";
      }
      system.fieldIo().writeFieldsBasis("out/w_1D_2.bf", system.w().basis(), 
                                        system.unitCell());
      TEST_ASSERT(diff < epsilon);
   }

   // Read parameter file to create a System object
   template <int D>
   void setUpSystem(System<D>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   // Read parameter file section to create a FilmIterator object
   template <int D>
   void setUpFilmIterator(FilmIterator<D, AmIteratorBasis<D>>& iterator, 
                          std::string fname)
   {
      std::ifstream in;
      openInputFile(fname, in);
      iterator.readParam(in);
      in.close();
   }

   // Determine if we get expected result when running checkSpaceGroup.
   // Function accepts a boolean indicating whether we expect it to throw 
   // an error or not, and returns a boolean indicating whether the 
   // function demonstrated the expected behavior.
   template <int D>
   bool checkCheckSpaceGroup(FilmIterator<D, AmIteratorBasis<D>>& iterator, 
                             bool expectError)
   {
      bool pass = true;
      if (expectError) {
         try {
            iterator.checkSpaceGroup();
            // This is expected to fail. If it succeeds, the test fails.
            pass = false;
         } catch (Exception& e) {
            Log::file() << "EXCEPTION CAUGHT, expected behavior occurred" 
                        << std::endl;
         }
      } else {
         try {
            iterator.checkSpaceGroup();
            // This should succeed. If not, the test fails. 
         } catch (Exception& e) {
            pass = false;
         }
      }
      return pass;
   }

};

TEST_BEGIN(FilmIteratorTest)
TEST_ADD(FilmIteratorTest, testConstructor)
TEST_ADD(FilmIteratorTest, testReadParameters)
TEST_ADD(FilmIteratorTest, testGenerateWallFields)
TEST_ADD(FilmIteratorTest, testCheckSpaceGroup1DA)
TEST_ADD(FilmIteratorTest, testCheckSpaceGroup1DB)
TEST_ADD(FilmIteratorTest, testCheckSpaceGroup2D)
TEST_ADD(FilmIteratorTest, testCheckSpaceGroup3DA)
TEST_ADD(FilmIteratorTest, testCheckSpaceGroup3DB)
TEST_ADD(FilmIteratorTest, testFlexibleParams)
TEST_ADD(FilmIteratorTest, testReadFlexibleParams)
TEST_ADD(FilmIteratorTest, testCheckLatticeVectors)
TEST_ADD(FilmIteratorTest, testSolve1D)
TEST_ADD(FilmIteratorTest, testSolve2D)
TEST_ADD(FilmIteratorTest, testSweep)
TEST_ADD(FilmIteratorTest, testFreeEnergy)
TEST_ADD(FilmIteratorTest, testMaskAndH)
TEST_END(FilmIteratorTest)

#endif
