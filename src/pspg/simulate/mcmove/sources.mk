pspg_simulate_mcmove_= \
  pspg/simulate/mcmove/McMove.cu \
  pspg/simulate/mcmove/McMoveFactory.cu \
  pspg/simulate/mcmove/McMoveManager.cu \
  pspg/simulate/mcmove/RealMove.cu \
  pspg/simulate/mcmove/FourierMove.cu \
  
pspg_simulate_mcmove_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_simulate_mcmove_:.cu=.o))

