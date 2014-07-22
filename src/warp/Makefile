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

OBJS =  patchIO.o\
#	computeSi.o\
	computeWi.o\
	updateVolumeCoordinates.o\


#FILES_TO_COMPLEXIFY =  

#FILES_TO_COPY = 


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
		if [ $$file -nt ../../src_cs/warp/$$file -o ! -e ../../src_cs/warp/$$file ]; then \
			echo complexifying file $$file; \
			python ../complexify.py $$file; \
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