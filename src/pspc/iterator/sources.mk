pspc_iterator_= \
  pspc/iterator/IteratorFactory.cpp \
  pspc/iterator/IteratorMediatorCPU.cpp \
  pspc/iterator/AmStrategyCPU.cpp \

pspc_iterator_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_iterator_))
pspc_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_iterator_:.cpp=.o))

