include $(SRC_DIR)/pspc/field/sources.mk
include $(SRC_DIR)/pspc/solvers/sources.mk
include $(SRC_DIR)/pspc/iterator/sources.mk
include $(SRC_DIR)/pspc/sweep/sources.mk
include $(SRC_DIR)/pspc/compressor/sources.mk
include $(SRC_DIR)/pspc/simulate/sources.mk

pspc_= \
  $(pspc_field_) \
  $(pspc_solvers_) \
  $(pspc_iterator_) \
  $(pspc_sweep_) \
  $(pspc_compressor_) \
  $(pspc_simulate_) \
  pspc/System.cpp 

pspc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_:.cpp=.o))

$(pspc_LIB): $(pspc_OBJS)
	$(AR) rcs $(pspc_LIB) $(pspc_OBJS)

