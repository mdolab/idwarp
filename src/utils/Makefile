MAIN_DIR = ../../
COMMON_FILE = $(MAIN_DIR)/Common.mk
RULES_FILE  = $(MAIN_DIR)/Rules.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FILES_TO_COMPLEXIFY = \
	initPETSc.F90\
	pointReduce.f90\
	vectorUtils.f90\

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
		if [ $$file -nt ../../src_cs/utils/$$file -o ! -e ../../src_cs/utils/$$file ]; then \
			echo complexifying file $$file; \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/utils/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/utils/$$file -o ! -e ../../src_cs/utils/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/utils/$$file; \
		fi;\
		done


