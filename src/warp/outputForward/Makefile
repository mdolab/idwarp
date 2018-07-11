#*********************************************************
# Makefile for IDWarp
# C.A.Mader       July 10, 2014
# G. K. W. Kenway Sept 30, 2014
#*********************************************************

MAIN_DIR = ../../..
COMMON_FILE = $(MAIN_DIR)/Common.mk
RULES_FILE  = $(MAIN_DIR)/Rules.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS =	cross_product_3d_d.o\
		getelementprops_d.o\
		getmag_d.o\
		getrotationmatrix3d_d.o

all:
	$(PMAKE)  src

clean:
	@echo "        Making clean in warp/outputForward... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src:	$(FF90_OBJECTS)
