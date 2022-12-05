include $(SRC_DIR)/pscf/chem/sources.mk
include $(SRC_DIR)/pscf/inter/sources.mk
include $(SRC_DIR)/pscf/math/sources.mk
include $(SRC_DIR)/pscf/mesh/sources.mk
include $(SRC_DIR)/pscf/crystal/sources.mk
include $(SRC_DIR)/pscf/homogeneous/sources.mk
include $(SRC_DIR)/pscf/iterator/sources.mk

pscf_= \
  $(pscf_chem_) $(pscf_inter_) $(pscf_math_) \
  $(pscf_mesh_) $(pscf_crystal_) $(pscf_homogeneous_) \
  $(pscf_iterator)

pscf_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_))
pscf_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_:.cpp=.o))

$(pscf_LIB): $(pscf_OBJS)
	$(AR) rcs $(pscf_LIB) $(pscf_OBJS)

