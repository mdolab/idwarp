#*********************************************************
# Makefile for IDWarp
# C.A.Mader       July 10, 2014
# G. K. W. Kenway Sept 30, 2014
#*********************************************************

MAIN_DIR = ../../
COMMON_FILE = $(MAIN_DIR)/Common.mk
RULES_FILE  = $(MAIN_DIR)/Rules.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

# Name of files to be compiled...Order here matters!
FF90_OBJECTS_1 = precision.o\

FF90_OBJECTS_2 = constants.o\

FF90_OBJECTS_3 = gridInput.o\
		 communication.o\

FF90_OBJECTS_4 = gridData.o\

FF90_OBJECTS_5 = kd_tree.o\
		 kd_tree2.o\
		 cgnsGrid.o\

FILES_TO_COMPLEXIFY =	precision.F90\
			constants.F90\
			gridInput.f90\
			communication.f90\
			gridData.F90\
			kd_tree.F90\
			cgnsGrid.F90\
			kd_tree2.f90\

default: all

all:
	make $(FF90_OBJECTS_1)
	make $(FF90_OBJECTS_2)
	make $(FF90_OBJECTS_3)
	make $(FF90_OBJECTS_4)
	make $(FF90_OBJECTS_5)
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in modules ... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		if [ $$file -nt ../../src_cs/modules/$$file -o ! -e ../../src_cs/modules/$$file ]; then \
			echo complexifying file $$file; \
			complexify $$file; \
			mv c_$$file ../../src_cs/modules/$$file; \
		fi;\
		done
	cp Makefile ../../src_cs/modules
