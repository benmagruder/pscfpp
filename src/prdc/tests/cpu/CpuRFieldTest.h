#ifndef PRDC_CPU_R_FIELD_TEST_H
#define PRDC_CPU_R_FIELD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cpu/RField.h>

#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>

using namespace Util;
using namespace Pscf::Prdc;

class CpuRFieldTest : public UnitTest 
{

private:

   typedef double Data;

public:

   void setUp() {}

   void tearDown() {}

   void testConstructor();
   void testAllocate3();
   void testSubscript();
   void testCopyConstructor();
   void testAssignment();
   void testSerialize1Memory();
   void testSerialize2Memory();
   void testSerialize1File();
   void testSerialize2File();

};


void CpuRFieldTest::testConstructor()
{

   printMethod(TEST_FUNC);
   {
      Cpu::RField<3> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   }
} 

void CpuRFieldTest::testAllocate3()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 2;
      d[1] = 3;
      d[2] = 4;

      Cpu::RField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 24);
      TEST_ASSERT(v.meshDimensions() == d);
   }
}

void CpuRFieldTest::testSubscript()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 2;
      d[1] = 3;
      d[2] = 4;

      Cpu::RField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 24);
      TEST_ASSERT(v.meshDimensions() == d);

      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0 ;
      }
   
      TEST_ASSERT(v[0] == 10.0);
      TEST_ASSERT(v[2] == 30.0);
   }
}
 
void CpuRFieldTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);

   {
      IntVec<3> d;
      d[0] = 2;
      d[1] = 3;
      d[2] = 4;

      Cpu::RField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 24);
      TEST_ASSERT(v.meshDimensions() == d);
   
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }
   
      Cpu::RField<3> u(v);
      TEST_ASSERT(u.capacity() == capacity);
      TEST_ASSERT(u.isAllocated() );
   
      TEST_ASSERT(v[0] == 10.0);
      TEST_ASSERT(v[2] == 30.0);
      TEST_ASSERT(u[0] == 10.0);
      TEST_ASSERT(u[2] == 30.0);
   }
} 

void CpuRFieldTest::testAssignment()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 2;
      d[1] = 3;
      d[2] = 4;

      Cpu::RField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 24);
      TEST_ASSERT(v.meshDimensions() == d);
   
      Cpu::RField<3> u;
      u.allocate(d);
      TEST_ASSERT(u.isAllocated() );
      TEST_ASSERT(u.capacity() == capacity);
   
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }
   
      u  = v;
   
      TEST_ASSERT(u.capacity() == capacity );
      TEST_ASSERT(u.isAllocated() );
      TEST_ASSERT(v[0] == 10.0);
      TEST_ASSERT(v[2] == 30.0);
      TEST_ASSERT(u[0] == 10.0);
      TEST_ASSERT(u[2] == 30.0);
   }
} 

void CpuRFieldTest::testSerialize1Memory()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 2;
      d[1] = 3;
      d[2] = 4;

      Cpu::RField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 24);
      TEST_ASSERT(v.meshDimensions() == d);
   
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }
      int size = memorySize(v);
     
      int i1 = 13;
      int i2;
   
      MemoryOArchive oArchive;
      oArchive.allocate(size + 12);
   
      oArchive << v;
      TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size);
      oArchive << i1;
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1]==20.0);
      TEST_ASSERT(v.capacity() == capacity);
   
      Cpu::RField<3> u;
      u.allocate(d);
   
      MemoryIArchive iArchive;
      iArchive = oArchive;
      TEST_ASSERT(iArchive.begin()  == oArchive.begin());
      TEST_ASSERT(iArchive.cursor() == iArchive.begin());
   
      // Load into u and i2
      iArchive >> u;
      TEST_ASSERT(iArchive.begin() == oArchive.begin());
      TEST_ASSERT(iArchive.end() == oArchive.cursor());
      TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);
   
      iArchive >> i2;
      TEST_ASSERT(iArchive.cursor() == iArchive.end());
      TEST_ASSERT(iArchive.begin() == oArchive.begin());
      TEST_ASSERT(iArchive.end() == oArchive.cursor());
   
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == capacity);
   
      // Release
      iArchive.release();
      TEST_ASSERT(!iArchive.isAllocated());
      TEST_ASSERT(iArchive.begin() == 0);
      TEST_ASSERT(iArchive.cursor() == 0);
      TEST_ASSERT(iArchive.end() == 0);
      TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size + sizeof(int));
   
      // Clear values of u and i2
      for (int i=0; i < capacity; i++ ) {
         u[i] = 0.0;
      }
      i2 = 0;
   
      // Reload into u and i2
      iArchive = oArchive;
      iArchive >> u;
      TEST_ASSERT(iArchive.begin() == oArchive.begin());
      TEST_ASSERT(iArchive.end() == oArchive.cursor());
      TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);
   
      iArchive >> i2;
      TEST_ASSERT(iArchive.cursor() == iArchive.end());
      TEST_ASSERT(iArchive.begin() == oArchive.begin());
      TEST_ASSERT(iArchive.end() == oArchive.cursor());
   
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == capacity);
   }

}

void CpuRFieldTest::testSerialize2Memory()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 2;
      d[1] = 3;
      d[2] = 4;

      Cpu::RField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 24);
      TEST_ASSERT(v.meshDimensions() == d);
   
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }
      int size = memorySize(v);
     
      MemoryOArchive oArchive;
      oArchive.allocate(size);
   
      oArchive << v;
      TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size);
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1] == 20.0);
      TEST_ASSERT(v.capacity() == capacity);
   
      Cpu::RField<3> u;
   
      // Note: We do not allocate RField u in this test.
      // This is the main difference from testSerialize1Memory()
   
      MemoryIArchive iArchive;
   
      iArchive = oArchive;
   
      TEST_ASSERT(iArchive.begin()  == oArchive.begin());
      TEST_ASSERT(iArchive.cursor() == iArchive.begin());
   
      iArchive >> u;
   
      TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(u.capacity() == capacity);
   }
}

void CpuRFieldTest::testSerialize1File()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 2;
      d[1] = 3;
      d[2] = 4;

      Cpu::RField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 24);
      TEST_ASSERT(v.meshDimensions() == d);
   
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }
     
      int i1 = 13;
      int i2;

      BinaryFileOArchive oArchive;
      openOutputFile("out/binary", oArchive.file());
      oArchive << v;
      oArchive << i1;
      oArchive.file().close();
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1]==20.0);
      TEST_ASSERT(v.capacity() == capacity);
   
      Cpu::RField<3> u;
      u.allocate(d);
   
      BinaryFileIArchive iArchive;
      openInputFile("out/binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
      iArchive.file().close();
   
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == capacity);
   
      // Clear values of u and i2
      for (int i=0; i < capacity; i++ ) {
         u[i] = 0.0;
      }
      i2 = 0;
   
      // Reload into u and i2
      openInputFile("out/binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
   
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == capacity);
   }
}

void CpuRFieldTest::testSerialize2File()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 2;
      d[1] = 3;
      d[2] = 4;

      Cpu::RField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 24);
      TEST_ASSERT(v.meshDimensions() == d);
   
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }
     
      int i1 = 13;
      int i2;
  
      BinaryFileOArchive oArchive;
      openOutputFile("out/binary", oArchive.file());
      oArchive << v;
      oArchive << i1;
      oArchive.file().close();
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1] == 20.0);
      TEST_ASSERT(v.capacity() == capacity);
   
      Cpu::RField<3> u;
   
      // u.allocate(d); -> 
      // Note: We do not allocate first. This is the difference 
      // from the previous test
   
      BinaryFileIArchive iArchive;
      openInputFile("out/binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
      iArchive.file().close();
   
      TEST_ASSERT(eq(u[1], 20.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == capacity);
   
      // Clear values of u and i2
      for (int i=0; i < capacity; i++ ) {
         u[i] = 0.0;
      }
      i2 = 0;
   
      // Reload into u and i2
      openInputFile("out/binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
   
      TEST_ASSERT(eq(u[1], 20.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == capacity);
   }
}

TEST_BEGIN(CpuRFieldTest)
TEST_ADD(CpuRFieldTest, testConstructor)
TEST_ADD(CpuRFieldTest, testAllocate3)
TEST_ADD(CpuRFieldTest, testSubscript)
TEST_ADD(CpuRFieldTest, testCopyConstructor)
TEST_ADD(CpuRFieldTest, testAssignment)
TEST_ADD(CpuRFieldTest, testSerialize1Memory)
TEST_ADD(CpuRFieldTest, testSerialize2Memory)
TEST_ADD(CpuRFieldTest, testSerialize1File)
TEST_ADD(CpuRFieldTest, testSerialize2File)
TEST_END(CpuRFieldTest)

#endif
