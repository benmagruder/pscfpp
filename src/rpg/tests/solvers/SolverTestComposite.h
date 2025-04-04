#ifndef RPG_TEST_SOLVER_TEST_COMPOSITE_H
#define RPG_TEST_SOLVER_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "PropagatorTest.h"
#include "BlockTest.h"
#include "MixtureTest.h"
#include "WaveListTest.h"

TEST_COMPOSITE_BEGIN(SolverTestComposite)
TEST_COMPOSITE_ADD_UNIT(WaveListTest);
TEST_COMPOSITE_ADD_UNIT(PropagatorTest);
TEST_COMPOSITE_ADD_UNIT(BlockTest);
TEST_COMPOSITE_ADD_UNIT(MixtureTest);
TEST_COMPOSITE_END

#endif
