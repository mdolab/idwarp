#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Gaetan Kenway                                   *
#      * Starting date: 05-23-2010                                      *
#      * Last modified: 12-23-2012                                      *
#      *                                                                *
#      ******************************************************************

SUBDIR_SRC    = src/modules	\
		src/warp	\
	        src/IO		\
	        src/utils 

WARP_SUBDIRS       = $(SUBDIR_SRC)
WARP_CLEAN_SUBDIRS = $(SUBDIR_SRC)

#      ******************************************************************
#      *                                                                *
#      * General targets.                                               *
#      *                                                                *
#      ******************************************************************

default:
	@echo "Usage: make <arch>"
	@echo "Supported architectures: LINUX_GFORTRAN_OPENMPI"
	@echo "                         LINUX_INTEL_OPENMPI"
	@echo "                         LINUX_INTEL_OPENMPI_SCINET"

all:	 default

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

#      ******************************************************************
#      *                                                                *
#      * The actual make. This is not a direct target, but is called    *
#      * from the architectures.                                        *
#      *                                                                *
#      ******************************************************************

warp:
	@for subdir in $(WARP_SUBDIRS) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make) || exit 1; \
		done
	(cd lib && make)

#      ******************************************************************
#      *                                                                *
#      * Platform specific targets.                                     *
#      *                                                                *
#      ******************************************************************

LINUX_GFORTRAN_OPENMPI:
	mkdir -p obj
	mkdir -p mod
	if [ ! -f "config/config.LINUX_GFORTRAN_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_GFORTRAN_OPENMPI.mk" ./config; fi
	ln -sf config/config.LINUX_GFORTRAN_OPENMPI.mk config.mk
	make warp
	(cd src/f2py && make)

LINUX_INTEL_OPENMPI:
	mkdir -p obj
	mkdir -p mod
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL_OPENMPI.mk config.mk
	make warp
	(cd src/f2py && make)

LINUX_INTEL_OPENMPI_SCINET:
	mkdir -p obj
	mkdir -p mod
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI_SCINET.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI_SCINET.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL_OPENMPI_SCINET.mk config.mk
	make warp
	(cd src/f2py && make)

