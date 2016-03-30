#!/usr/bin/make

# defaults
STATIC = no
COMPILER = gnu

#main building variables
ifeq "$(STATIC)" "yes"
  DOBJ    = static/obj/
  DMOD    = static/mod/
  DEXE    = static/
  MAKELIB = ar -rcs $(DEXE)libwenoof.a $(DOBJ)*.o ; ranlib $(DEXE)libwenoof.a
  RULE    = WENOOF
else
  DOBJ = tests/obj/
  DMOD = tests/mod/
  DEXE = tests/
  RULE = $(DEXE)sin_reconstruction $(DEXE)modules_imports
endif
DSRC = src/
LIBS =
ifeq "$(COMPILER)" "gnu"
  FC    = gfortran
  OPTSC = -cpp -c -frealloc-lhs -O2  -J $(DMOD)
  OPTSL = -J $(DMOD)
endif
ifeq "$(COMPILER)" "ibm"
  FC    = bgxlf2008_r
  OPTSC = -c -O2 -qmoddir=$(DMOD) -I$(DMOD)
  OPTSL = -qmoddir=$(DMOD) -I$(DMOD)
endif
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"
RUTEXT = "Executed rule $@"

firsrule: $(RULE)

#building rules
$(DEXE)sin_reconstruction: $(MKDIRS) $(DOBJ)sin_reconstruction.o
	@rm -f $(filter-out $(DOBJ)sin_reconstruction.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) sin_reconstruction

$(DEXE)modules_imports: $(MKDIRS) $(DOBJ)modules_imports.o
	@rm -f $(filter-out $(DOBJ)modules_imports.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) modules_imports

WENOOF: $(MKDIRS) $(DOBJ)wenoof.o
	@echo $(LITEXT)
	@$(MAKELIB)

#compiling rules
$(DOBJ)wenoof.o: src/lib/wenoof.f90 \
	$(DOBJ)type_weno_interpolator.o \
	$(DOBJ)type_weno_interpolator_upwind.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)type_weno_interpolator.o: src/lib/type_weno_interpolator.f90 \
	$(DOBJ)penf.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)type_weno_interpolator_upwind.o: src/lib/type_weno_interpolator_upwind.f90 \
	$(DOBJ)penf.o \
	$(DOBJ)type_weno_interpolator.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)pyplot_module.o: src/third_party/pyplot-fortran/src/pyplot_module.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)penf.o: src/third_party/PENF/src/lib/penf.F90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sin_reconstruction.o: src/tests/sin_reconstruction.f90 \
	$(DOBJ)penf.o \
	$(DOBJ)wenoof.o \
	$(DOBJ)pyplot_module.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)modules_imports.o: src/tests/modules_imports.f90 \
	$(DOBJ)wenoof.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
