pspc_simlate_trajectory_= \
  pspc/simulate/trajectory/TrajectoryReader.cpp \
  pspc/simulate/trajectory/TrajectoryReaderFactory.cpp \
  pspc/simulate/trajectory/FieldConfigReader.cpp 

pspc_simulate_trajectory_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_simulate_trajectory_))
pspc_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_simulate_trajectory_:.cpp=.o))

