#*********************************************************
# Makefile for pywarpustruct
# C.A.Mader  July 10, 2014
#*********************************************************

MAIN_DIR = ../../
COMMON_FILE = $(MAIN_DIR)/Common.mk
RULES_FILE  = $(MAIN_DIR)/Rules.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

# Name of files to be compiled:
#MODS0 = 

OBJS = readUnstructuredCGNS.o\
       writeUnstructuredCGNS.o\

#FILES_TO_COMPLEXIFY =  enginemodules.F90\
		       kriging_data.F90\
	               weightmodules.F90\
		       krigingInputs.F90\
	               krigingModels.F90\
		       missionParameters.F90\
		       aircraftParameters.F90\
		       computationParameters.f90\
		       conversion.f90\
		       missionParameters.F90\

#FILES_TO_COPY = complexify.f90\
		diffSizes.f90\
		communication.f90\
	        diffSizes.f90\
	       	dvParameters.f90\
		modulesPETSc.F90\
		precision.f90\
		rbf_data.f90\

default:
#	$(MAKE) $(MODS0)
#	$(MAKE) $(MODS1)
#	cp -f *.mod $(MODDIR)
	$(MAKE) $(OBJS)

clean:
	@echo "        Making clean in src... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		if [ $$file -nt ../../src_cs/modules/$$file -o ! -e ../../src_cs/modules/$$file ]; then \
			echo complexifying file $$file; \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/modules/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/modules/$$file -o ! -e ../../src_cs/modules/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/modules/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/modules