#ifndef PSPG_FIELD_TEST_COMPOSITE_H
#define PSPG_FIELD_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "FftTest.h"
#include "FieldIoTest.h"
#include "DomainTest.h"

TEST_COMPOSITE_BEGIN(FieldTestComposite)
TEST_COMPOSITE_ADD_UNIT(FftTest);
TEST_COMPOSITE_ADD_UNIT(DomainTest);
TEST_COMPOSITE_ADD_UNIT(FieldIoTest);
TEST_COMPOSITE_END

#endif
