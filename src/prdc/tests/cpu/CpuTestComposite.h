#ifndef PRDC_CPU_TEST_COMPOSITE_H
#define PRDC_CPU_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CpuFieldTest.h"
#include "CpuRFieldTest.h"
#include "CpuRFieldDftTest.h"
#include "CpuFieldComparisonTest.h"
#include "CpuFftTest.h"
#include "CpuFieldBasisConverterTest.h"

TEST_COMPOSITE_BEGIN(CpuTestComposite)
TEST_COMPOSITE_ADD_UNIT(CpuFieldTest);
TEST_COMPOSITE_ADD_UNIT(CpuRFieldTest);
TEST_COMPOSITE_ADD_UNIT(CpuRFieldDftTest);
TEST_COMPOSITE_ADD_UNIT(CpuFieldComparisonTest);
TEST_COMPOSITE_ADD_UNIT(CpuFftTest);
TEST_COMPOSITE_ADD_UNIT(CpuFieldBasisConverterTest);
TEST_COMPOSITE_END

#endif
