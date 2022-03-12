include src/config.mk
# =========================================================================
.PHONY:  html html-graphs clean clean-tests clean-bin clean-html veryclean 

# =========================================================================
# Main targets 
 
all-cpu:
	cd bld; $(MAKE) all-cpu

all:
	cd bld; $(MAKE) all

fd1d:
	cd bld; $(MAKE) fd1d

pspc:
	cd bld; $(MAKE) pspc

pspg:
	cd bld; $(MAKE) pspg


# =========================================================================
# HTML Documentation
 
html:
	cd docs; $(MAKE) html

html-dev:
	cd docs; $(MAKE) html-dev

# =========================================================================
# Clean targets

clean:
	cd bld; $(MAKE) clean
	cd src; $(MAKE) clean

clean-tests:
	cd src/; $(MAKE) clean-tests
	cd bld/; $(MAKE) clean-tests

clean-bin:
	rm -f $(BIN_DIR)/pscf*
 
clean-html:
	cd docs; $(MAKE) clean

veryclean:
	make clean-bin
	rm -f $(BIN_DIR)/makeDep*
	cd bld; $(MAKE) veryclean
	rm bld/makefile
	cd src; $(MAKE) veryclean
	cd docs; $(MAKE) clean
	cd examples; ./clean

# ==========================================================================
