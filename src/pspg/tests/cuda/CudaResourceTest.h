#ifndef PSPG_CUDA_RESOURCE_TEST_H
#define PSPG_CUDA_RESOURCE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspg/field/RDField.h>
#include <util/math/Constants.h>
#include <pspg/math/GpuResources.h>

#include <cstdlib>
#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspg;

class CudaResourceTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testReductionSum() 
   {
      printMethod(TEST_FUNC);

      // GPU Resources
      NUMBER_OF_BLOCKS = 32; // parallel reduction into 32 blocks.
      THREADS_PER_BLOCK = 32;

      // Create device and host arrays
      const int n = NUMBER_OF_BLOCKS*THREADS_PER_BLOCK*2;
      cudaReal sum = 0;
      cudaReal sumCheck = 0;
      cudaReal* num = new cudaReal[n];
      cudaReal* d_temp;
      cudaReal* d_num;
      cudaMalloc((void**) &d_num, n*sizeof(cudaReal));
      cudaMalloc((void**) &d_temp, NUMBER_OF_BLOCKS*sizeof(cudaReal));

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 10000);
      }

      // Host find max
      for (int i = 0; i < n; i++) {
         sumCheck+=num[i];
      }

      // Launch kernel twice and get output
      cudaMemcpy(d_num, num, n*sizeof(cudaReal), cudaMemcpyHostToDevice);
      reductionSum<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(cudaReal)>>>(d_temp, d_num, n);
      reductionSum<<<1, NUMBER_OF_BLOCKS/2, NUMBER_OF_BLOCKS/2*sizeof(cudaReal)>>>(d_temp, d_temp, NUMBER_OF_BLOCKS);
      cudaMemcpy(&sum, d_temp, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      TEST_ASSERT(sum == sumCheck);
   }

   void testReductionMaxSmall()
   {
      printMethod(TEST_FUNC);

      // GPU Resources
      NUMBER_OF_BLOCKS = 1; // parallel reduction into 1 block.
      THREADS_PER_BLOCK = 32;

      // Create device and host arrays
      const int n = 1*32*2;
      cudaReal max = -1;
      cudaReal maxCheck = -10;
      cudaReal* num = new cudaReal[n];
      cudaReal* d_max;
      cudaReal* d_num;
      cudaMalloc((void**) &d_max, 1*sizeof(cudaReal));
      cudaMalloc((void**) &d_num, n*sizeof(cudaReal));

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 100);
      }
      cudaMemcpy(d_num, num, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Host find max
      maxCheck = 0;
      for (int i = 0; i < n; i++) {
         if (num[i] > maxCheck) {
            maxCheck = num[i];
         }
      }

      // Launch kernel and get output
      reductionMax<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(cudaReal)>>>(d_max, d_num, n);
      cudaMemcpy(&max, d_max, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      TEST_ASSERT(max == maxCheck);

   }

   void testReductionMaxLarge()
   {
      printMethod(TEST_FUNC);

      // GPU Resources
      NUMBER_OF_BLOCKS = 8; // parallel reduction into 32 blocks.
      THREADS_PER_BLOCK = 32;

      // Create device and host arrays
      const int n = NUMBER_OF_BLOCKS*THREADS_PER_BLOCK*2;
      cudaReal max = -1;
      cudaReal maxCheck = -10;
      cudaReal* num = new cudaReal[n];
      cudaReal* d_temp;
      cudaReal* d_max;
      cudaReal* d_num;
      cudaMalloc((void**) &d_num, n*sizeof(cudaReal));
      cudaMalloc((void**) &d_temp, NUMBER_OF_BLOCKS*sizeof(cudaReal));
      cudaMalloc((void**) &d_max, 1*sizeof(cudaReal));

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 10000);
      }

      // Host find max
      for (int i = 0; i < n; i++) {
         if (num[i] > maxCheck) {
            maxCheck = num[i];
         }
      }

      // Launch kernel twice and get output
      cudaMemcpy(d_num, num, n*sizeof(cudaReal), cudaMemcpyHostToDevice);
      reductionMax<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(cudaReal)>>>(d_temp, d_num, n);
      reductionMax<<<1, NUMBER_OF_BLOCKS/2, NUMBER_OF_BLOCKS/2*sizeof(cudaReal)>>>(d_max, d_temp, NUMBER_OF_BLOCKS);
      cudaMemcpy(&max, d_max, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      TEST_ASSERT(max == maxCheck);

   }

   void testReductionMaxAbsLarge()
   {
      printMethod(TEST_FUNC);

      // GPU Resources
      NUMBER_OF_BLOCKS = 8; // parallel reduction into 32 blocks.
      THREADS_PER_BLOCK = 32;

      // Create device and host arrays
      const int n = NUMBER_OF_BLOCKS*THREADS_PER_BLOCK*2;
      cudaReal max = -1;
      cudaReal maxCheck = -10;
      cudaReal* num = new cudaReal[n];
      cudaReal* d_temp;
      cudaReal* d_max;
      cudaReal* d_num;
      cudaMalloc((void**) &d_num, n*sizeof(cudaReal));
      cudaMalloc((void**) &d_temp, NUMBER_OF_BLOCKS*sizeof(cudaReal));
      cudaMalloc((void**) &d_max, 1*sizeof(cudaReal));

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 10000 - 6000);
      }

      // Host find max
      for (int i = 0; i < n; i++) {
         if (fabs(num[i]) > maxCheck) {
            maxCheck = fabs(num[i]);
         }
      }

      // Launch kernel twice and get output
      cudaMemcpy(d_num, num, n*sizeof(cudaReal), cudaMemcpyHostToDevice);
      reductionMaxAbs<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(cudaReal)>>>(d_temp, d_num, n);
      reductionMaxAbs<<<1, NUMBER_OF_BLOCKS/2, NUMBER_OF_BLOCKS/2*sizeof(cudaReal)>>>(d_max, d_temp, NUMBER_OF_BLOCKS);
      cudaMemcpy(&max, d_max, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      TEST_ASSERT(max == maxCheck);

   }

   void testReductionMinLarge()
   {
      printMethod(TEST_FUNC);

      // GPU Resources
      NUMBER_OF_BLOCKS = 8; // parallel reduction into 32 blocks.
      THREADS_PER_BLOCK = 32;

      // Create device and host arrays
      const int n = NUMBER_OF_BLOCKS*THREADS_PER_BLOCK*2;
      cudaReal min = 100000;
      cudaReal minCheck = 100000;
      cudaReal* num = new cudaReal[n];
      cudaReal* d_temp;
      cudaReal* d_min;
      cudaReal* d_num;
      cudaMalloc((void**) &d_num, n*sizeof(cudaReal));
      cudaMalloc((void**) &d_temp, NUMBER_OF_BLOCKS*sizeof(cudaReal));
      cudaMalloc((void**) &d_min, 1*sizeof(cudaReal));

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 10000);
      }

      // Host find max
      for (int i = 0; i < n; i++) {
         if (num[i] < minCheck) {
            minCheck = num[i];
         }
      }

      // Launch kernel twice and get output
      cudaMemcpy(d_num, num, n*sizeof(cudaReal), cudaMemcpyHostToDevice);
      reductionMin<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(cudaReal)>>>(d_temp, d_num, n);
      reductionMin<<<1, NUMBER_OF_BLOCKS/2, NUMBER_OF_BLOCKS/2*sizeof(cudaReal)>>>(d_min, d_temp, NUMBER_OF_BLOCKS);
      cudaMemcpy(&min, d_min, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      TEST_ASSERT(min == minCheck);

   }

   void testGpuInnerProduct()
   {
      printMethod(TEST_FUNC);
      // GPU Resources
      MAX_THREADS_PER_BLOCK = 128;

      // Data size
      // Non-power-of-two to check performance in weird situations
      int n = 14*MAX_THREADS_PER_BLOCK+77;

      // Device arrays
      DField<cudaReal> d_a, d_b;
      d_a.allocate(n);
      d_b.allocate(n);

      // Host arrays
      cudaReal* a = new cudaReal[n];
      cudaReal* b = new cudaReal[n];
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++ ) {
         a[i] = (cudaReal)(std::rand() % 10000);
         b[i] = (cudaReal)(std::rand() % 10000);
      }
      cudaMemcpy(d_a.cDField(), a, n*sizeof(cudaReal), cudaMemcpyHostToDevice);
      cudaMemcpy(d_b.cDField(), b, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal prodCheck = 0;
      for (int i = 0; i < n; i++) {
         prodCheck += a[i]*b[i];
      }

      // Inner product on device
      cudaReal prod = gpuInnerProduct(d_a.cDField(),d_b.cDField(),n);

      TEST_ASSERT(prodCheck==prod);
   }

   void testGpuSum()
   {
      printMethod(TEST_FUNC);
      // GPU Resources
      MAX_THREADS_PER_BLOCK = 128;

      // Data size
      // Non-power-of-two to check performance in weird situations
      int n = 14*MAX_THREADS_PER_BLOCK+77;

      // Device arrays
      DField<cudaReal> d_data;
      d_data.allocate(n);

      // Host arrays
      cudaReal* data = new cudaReal[n];
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000);
      }
      cudaMemcpy(d_data.cDField(), data, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal prodCheck = 0;
      for (int i = 0; i < n; i++) {
         prodCheck += data[i];
      }

      // Inner product on device
      cudaReal prod = gpuSum(d_data.cDField(),n);

      TEST_ASSERT(prodCheck==prod);
   }

   void testGpuMax()
   {
      printMethod(TEST_FUNC);
      // GPU Resources
      MAX_THREADS_PER_BLOCK = 128;

      // Data size
      // Non-power-of-two to check performance in weird situations
      int n = 14*MAX_THREADS_PER_BLOCK+77;

      // Device arrays
      DField<cudaReal> d_data;
      d_data.allocate(n);

      // Host arrays
      cudaReal* data = new cudaReal[n];
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000);
      }
      cudaMemcpy(d_data.cDField(), data, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal maxCheck = 0;
      for (int i = 0; i < n; i++) {
         if (data[i] > maxCheck) maxCheck = data[i];
      }

      // Inner product on device
      cudaReal max = gpuMax(d_data.cDField(),n);

      TEST_ASSERT(max==maxCheck);
   }

   void testGpuMaxAbs()
   {
      printMethod(TEST_FUNC);
      // GPU Resources
      MAX_THREADS_PER_BLOCK = 128;

      // Data size
      // Non-power-of-two to check performance in weird situations
      int n = 14*MAX_THREADS_PER_BLOCK+77;

      // Device arrays
      DField<cudaReal> d_data;
      d_data.allocate(n);

      // Host arrays
      cudaReal* data = new cudaReal[n];
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000 - 6000);
      }
      cudaMemcpy(d_data.cDField(), data, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal maxCheck = 0;
      for (int i = 0; i < n; i++) {
         if (fabs(data[i]) > maxCheck) maxCheck = fabs(data[i]);
      }

      // Inner product on device
      cudaReal max = gpuMaxAbs(d_data.cDField(),n);

      TEST_ASSERT(max==maxCheck);
   }

   void testGpuMin()
   {
      printMethod(TEST_FUNC);
      // GPU Resources
      MAX_THREADS_PER_BLOCK = 128;

      // Data size
      // Non-power-of-two to check performance in weird situations
      int n = 14*MAX_THREADS_PER_BLOCK+77;

      // Device arrays
      DField<cudaReal> d_data;
      d_data.allocate(n);

      // Host arrays
      cudaReal* data = new cudaReal[n];
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000);
      }
      cudaMemcpy(d_data.cDField(), data, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal minCheck = 1000000;
      for (int i = 0; i < n; i++) {
         if (data[i] < minCheck) minCheck = data[i];
      }

      // Inner product on device
      cudaReal min = gpuMin(d_data.cDField(),n);

      TEST_ASSERT(min==minCheck);
   }

   void testGpuMinAbs()
   {
      printMethod(TEST_FUNC);
      // GPU Resources
      MAX_THREADS_PER_BLOCK = 128;

      // Data size
      // Non-power-of-two to check performance in weird situations
      int n = 14*MAX_THREADS_PER_BLOCK+77;

      // Device arrays
      DField<cudaReal> d_data;
      d_data.allocate(n);

      // Host arrays
      cudaReal* data = new cudaReal[n];
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000 - 6000);
      }
      cudaMemcpy(d_data.cDField(), data, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal minCheck = 1E300;
      for (int i = 0; i < n; i++) {
         if (fabs(data[i]) < minCheck) minCheck = fabs(data[i]);
      }

      // Inner product on device
      cudaReal min = gpuMinAbs(d_data.cDField(),n);

      TEST_ASSERT(min==minCheck);
   }

};

TEST_BEGIN(CudaResourceTest)
TEST_ADD(CudaResourceTest, testReductionSum)
TEST_ADD(CudaResourceTest, testReductionMaxSmall)
TEST_ADD(CudaResourceTest, testReductionMaxLarge)
TEST_ADD(CudaResourceTest, testReductionMaxAbsLarge)
TEST_ADD(CudaResourceTest, testReductionMinLarge)
TEST_ADD(CudaResourceTest, testGpuInnerProduct)
TEST_ADD(CudaResourceTest, testGpuSum)
TEST_ADD(CudaResourceTest, testGpuMax)
TEST_ADD(CudaResourceTest, testGpuMaxAbs)
TEST_ADD(CudaResourceTest, testGpuMin)
TEST_ADD(CudaResourceTest, testGpuMinAbs)
TEST_END(CudaResourceTest)

#endif
