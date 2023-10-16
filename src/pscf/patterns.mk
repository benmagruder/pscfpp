# ---------------------------------------------------------------------
# File: src/pscf/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/pscf directory, which
# contains all source code for the Util namespace. It is included by
# all "makefile" files in this directory tree. 
#
# This file should be included in other makefiles after inclusion of
# the files src/config.mk and src/pscf/config.mk because this file
# uses makefile variables defined in those files.
#-----------------------------------------------------------------------

# List of pscf-specific libraries needed in src/pscf (the order matters)
PSCF_LIBS=$(pscf_LIB) $(util_LIB) 

# List of all libraries needed for executables in src/pscf
LIBS=$(PSCF_LIBS)

# Add paths to Gnu Scientific Library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Conditionally add CUDA header include and library paths
ifdef PSCF_CUDA
  INCLUDES+=$(CUDA_INC)
  LIBS+=$(CUDA_LIB)
endif

# Preprocessor macro definitions needed in src/pscf
DEFINES=$(PSCF_DEFS) $(UTIL_DEFS)

# Dependencies on build configuration files
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pscf/config.mk

# Pattern rule to compile *.cpp C++ source files in src/pscf
# Note: Creates a *.d dependency file as a side effect of compilation
$(BLD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
   endif

# Pattern rule to compile *.cu CUDA source files in src/pscf
# Note: Creates a *.d dependency file as a side effect of compilation
$(BLD_DIR)/%.o: $(SRC_DIR)/%.cu
	$(NVXX) $(CPPFLAGS) $(NVXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP_CUDA
	$(MAKEDEP_CUDA) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
   endif

# Pattern rule to compile Test programs in src/pscf/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(PSCF_LIBS)
ifdef PSCF_CUDA
	$(NVXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $@ $< $(LIBS)
else
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $@ $< $(LIBS)
endif

