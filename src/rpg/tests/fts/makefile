BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/rpg/include.mk
include $(SRC_DIR)/rpg/tests/fts/sources.mk

TEST=rpg/tests/fts/Test

all: $(rpg_tests_fts_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(rpg_tests_fts_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpg/tests/fts/

run-quiet: $(rpg_tests_fts_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/rpg/tests/fts/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in rpg/tests/fts/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	rm -f out/*

clean:
	rm -f $(rpg_tests_fts_OBJS) $(rpg_tests_fts_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(rpg_tests_fts_OBJS:.o=.d)
-include $(rpg_tests_fts_OBJS:.o=.d)
