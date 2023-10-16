# ---------------------------------------------------------------------
# File: src/prdc/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/prdc directory, which
# contains all source code for the Pscf::Prdc namespace. It is included 
# by all "makefile" files in this directory tree. 
#
# This file must be included in other makefiles after inclusion of
# the root src/config.mk and relevant namespace level config.mk files 
# in the build directory, because this file uses makefile variables 
# defined in those configuration files.
#-----------------------------------------------------------------------

# Local pscf-specific libraries needed in src/prdc
PRDC_LIBS=$(prdc_LIB) $(pscf_LIB) $(util_LIB)

# All libraries needed in executables built in src/prdc
LIBS=$(PRDC_LIBS)

# Add paths to Gnu scientific library (GSL)
INCLUDES+=$(GSL_INC)
LIBS+=$(GSL_LIB) 

# Add paths to FFTW Fast Fourier transform library
INCLUDES+=$(FFTW_INC)
LIBS+=$(FFTW_LIB) 

# Add paths to CUDA FFT library
ifdef PSCF_CUDA
  INCLUDES+=$(CUDA_INC)
  LIBS+=$(CUDA_LIB)
endif

# Preprocessor macro definitions needed in src/prdc
DEFINES=$(UTIL_DEFS) $(PSCF_DEFS) $(PRDC_DEFS) 

# Dependencies on build configuration files
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/pscf/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/prdc/config.mk

# Pattern rule to compile *.cpp class source files in src/prdc
# Note: Creates a *.d dependency file as a side effect of compilation
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
   endif

# Pattern rule to compile *.cu class source files in src/prdc
# Note: Creates a *.d dependency file as a side effect of compilation
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cu
	$(NVXX) $(CPPFLAGS) $(NVXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
   ifdef MAKEDEP_CUDA
	$(MAKEDEP_CUDA) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
   endif

# Pattern rule to compile Test programs in src/prdc/tests
$(BLD_DIR)/%Test: $(BLD_DIR)/%Test.o $(PRDC_LIBS)
ifdef PSCF_CUDA
	$(NVXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $@ $< $(LIBS)
else
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $@ $< $(LIBS)
endif
