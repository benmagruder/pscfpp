#ifndef PRDC_CUDA_FIELD_COMPARISON_TEST_H
#define PRDC_CUDA_FIELD_COMPARISON_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cuda/RField.h>
#include <prdc/cuda/HostField.h>
#include <prdc/cuda/RFieldComparison.h>

#include <util/format/Dbl.h>

using namespace Util;
using namespace Pscf::Prdc;

class CudaFieldComparisonTest : public UnitTest 
{

public:

   void setUp()
   {  setVerbose(0);  }

   void tearDown() 
   {}

   void testRFieldComparison_1D()
   {
      printMethod(TEST_FUNC);

      // Define dimensions
      int n = 10;
      IntVec<1> dimensions;
      dimensions[0] = n;

      // Allocate arrays on host CPU
      Cuda::HostField<cudaReal> ha, hb;
      ha.allocate(n);
      hb.allocate(n);
      int size = ha.capacity();
      TEST_ASSERT(size == n);

      // Initialize data for two slightly different fields
      for (int i = 0; i < n; ++i) {
         ha[i] = 2.0;
         hb[i] = 2.001;
      }

      // Copy data to fields on the device
      Cuda::RField<1> da, db;
      da.allocate(dimensions);
      db.allocate(dimensions);
      da = ha;
      db = hb;

      // Make comparison
      Cuda::RFieldComparison<1> comparison;
      comparison.compare(da, db);


      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "MaxDiff = " 
                   << Dbl(comparison.maxDiff(), 20, 12) << "\n";
         std::cout << "RmsDiff = " 
                   << Dbl(comparison.rmsDiff(), 20, 12) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 0.0011);
      TEST_ASSERT(comparison.maxDiff() > 0.0009);
      TEST_ASSERT(comparison.rmsDiff() < 0.0011);
      TEST_ASSERT(comparison.rmsDiff() > 0.0009);

   }

   void testRFieldComparison_2D()
   {
      printMethod(TEST_FUNC);

      int m = 5;
      int n = 10;
      int size = m*n;

      IntVec<2> dimensions;
      dimensions[0] = m;
      dimensions[1] = n;

      // Allocate fields ha and hb on CPU host
      Cuda::HostField<cudaReal> ha, hb;
      ha.allocate(size);
      hb.allocate(size);
      TEST_ASSERT(ha.capacity() == size);
      TEST_ASSERT(hb.capacity() == size);

      // Initialize data
      for (int i = 0; i < size; ++i) {
         ha[i] = 2.0;
         hb[i] = 2.001;
      }

      // Copy data to fields da and db on the GPU device
      Cuda::RField<2> da, db;
      da.allocate(dimensions);
      db.allocate(dimensions);
      da = ha;
      db = hb;

      // Make comparison
      Cuda::RFieldComparison<2> comparison;
      comparison.compare(da, db);

      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "MaxDiff = " 
                   << Dbl(comparison.maxDiff(), 20, 12) << "\n";
         std::cout << "RmsDiff = " 
                   << Dbl(comparison.rmsDiff(), 20, 12) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 0.0011);
      TEST_ASSERT(comparison.maxDiff() > 0.0009);
      TEST_ASSERT(comparison.rmsDiff() < 0.0011);
      TEST_ASSERT(comparison.rmsDiff() > 0.0009);

   }

   void testRFieldArrayComparison_2D()
   {
      printMethod(TEST_FUNC);

      int nMonomer = 2;

      int m = 5;
      int n = 10;
      int meshSize = m*n;
      IntVec<2> dimensions;
      dimensions[0] = m;
      dimensions[1] = n;

      // Allocate fields ha and hb on CPU host
      DArray< Cuda::HostField<cudaReal> > ha, hb;
      ha.allocate(nMonomer);
      hb.allocate(nMonomer);

      // Initialize data
      for (int i = 0; i < nMonomer; ++i) {
         ha[i].allocate(meshSize);
         hb[i].allocate(meshSize);
         TEST_ASSERT(ha[i].capacity() == meshSize);
         TEST_ASSERT(hb[i].capacity() == meshSize);
         for (int j = 0; j < meshSize; ++j) {
            ha[i][j] = double(i) + 2.0;
            hb[i][j] = double(i) + 2.001;
         }
      }

      // Copy data to fields da and db on the GPU device
      DArray< Cuda::RField<2> > da;
      DArray< Cuda::RField<2> > db;
      da.allocate(nMonomer);
      db.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         da[i].allocate(dimensions);
         db[i].allocate(dimensions);
         da[i] = ha[i];
         db[i] = hb[i];
      }

      // Make comparison
      Cuda::RFieldComparison<2> comparison;
      comparison.compare(da, db);

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "MaxDiff = " 
                   << Dbl(comparison.maxDiff(), 20, 12) << "\n";
         std::cout << "RmsDiff = " 
                   << Dbl(comparison.rmsDiff(), 20, 12) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 0.0011);
      TEST_ASSERT(comparison.maxDiff() > 0.0009);
      TEST_ASSERT(comparison.rmsDiff() < 0.0011);
      TEST_ASSERT(comparison.rmsDiff() > 0.0009);

   }
};


TEST_BEGIN(CudaFieldComparisonTest)
TEST_ADD(CudaFieldComparisonTest, testRFieldComparison_1D)
TEST_ADD(CudaFieldComparisonTest, testRFieldComparison_2D)
TEST_ADD(CudaFieldComparisonTest, testRFieldArrayComparison_2D)
TEST_END(CudaFieldComparisonTest)

#endif