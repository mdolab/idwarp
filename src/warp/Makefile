#*********************************************************
# Makefile for IDWarp
# C.A.Mader  July 10, 2014
#*********************************************************

MAIN_DIR = ../../
COMMON_FILE = $(MAIN_DIR)/Common.mk
RULES_FILE  = $(MAIN_DIR)/Rules.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FILES_TO_COMPLEXIFY = \
	getdXs.F90\
	getElementProps.F90\
	getSurfaceCoordinates.F90\
	getVolumeCoordinates.F90\
	getCommonVolumeCoordinates.F90\
	gridTransform.F90\
	initializeWarping.F90\
	setdXs.F90\
	setExternalMeshIndices.F90\
	setSurfaceCoordinates.F90\
	verifyWarpDeriv.F90\
	warpMesh.F90\
	warpDeriv.F90\
	warpDerivFwd.F90\

FILES_TO_COPY =

ALL_FILES = $(FILES_TO_COMPLEXIFY) $(FILES_TO_COPY)
OFILES := $(ALL_FILES:%.f90=%.o) $(ALL_FILES:%.F90=%.o) $(ALL_FILES:%.c=%.o)

default: all

all:
	$(PMAKE) src

clean:
	@echo "        Making clean in IO ... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src : $(OFILES)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		if [ $$file -nt ../../src_cs/warp/$$file -o ! -e ../../src_cs/warp/$$file ]; then \
			echo complexifying file $$file; \
			complexify $$file; \
			mv c_$$file ../../src_cs/warp/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/warp/$$file -o ! -e ../../src_cs/warp/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/warp/$$file; \
		fi;\
		done
	cp Makefile ../../src_cs/warp
