#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Gaetan Kenway                                   *
#      * Starting date: 05-23-2010                                      *
#      * Last modified: 12-23-2012                                      *
#      *                                                                *
#      ******************************************************************

SUBDIR_SRC    = src/modules	\
	        src/IO		\
	        src/utils	\
		src/ADFirstAidKit\
		src/warp	\
		src/warp/outputReverse\


WARP_SUBDIRS       = $(SUBDIR_SRC)
WARP_CLEAN_SUBDIRS = $(SUBDIR_SRC)

default:
# Check if the config.mk file is in the config dir.
	@if [ ! -f "config/config.mk" ]; then \
	echo "Before compiling, copy an existing config file from the "; \
	echo "config/defaults/ directory to the config/ directory and  "; \
	echo "rename to config.mk. For example:"; \
	echo " ";\
	echo "  cp config/defaults/config.LINUX_INTEL_OPENMPI.mk config/config.mk"; \
	echo " ";\
	echo "The modify this config file as required. Typically the CGNS directory "; \
	echo "will have to be modified. With the config file specified, rerun "; \
	echo "'make' and the build will start"; \
	else make warp;\
	fi;

clean:
	@echo " Making clean ... "

	@for subdir in $(WARP_CLEAN_SUBDIRS) ; \
		do \
			echo; \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make $@) || exit 1; \
		done
	rm -f *~ config.mk;
	rm -f lib/lib* mod/* obj/*

warp:
	mkdir -p obj
	mkdir -p mod
	ln -sf config/config.mk config.mk
	ln -sf Common_real.mk Common.mk

	@for subdir in $(WARP_SUBDIRS) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make) || exit 1; \
		done
	(cd lib && make)
	(cd src/f2py && make)

