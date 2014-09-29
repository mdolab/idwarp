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
FF90_OBJECTS_1 = precision.o\

FF90_OBJECTS_2 = constants.o\

FF90_OBJECTS_3 = communication.o\
		 gridData.o\
		 gridInput.o\
		 kd_tree.o\
		 structuredCGNS.o\

FILES_TO_COMPLEXIFY =	precision.F90\
			constants.f90\
			communication.f90\
			gridData.f90\
			gridInput.f90\
			kd_tree.f90\
			structuredCGNS.F90\

FILES_TO_COPY = 

default: all

all: 
	make $(FF90_OBJECTS_1)
	make $(FF90_OBJECTS_2)
	make $(FF90_OBJECTS_3)
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in modules ... "
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
