pspc_solvers_= \
  pspc/solvers/Propagator.cpp \
  pspc/solvers/Block.cpp \
  pspc/solvers/Polymer.cpp \
  pspc/solvers/Solvent.cpp \
  pspc/solvers/Mixture.cpp

pspc_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_solvers_:.cpp=.o))

