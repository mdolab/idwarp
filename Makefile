# Makefile for pymissionanalysis

SUBDIR_SRC    = src/modules\
		src/IO\
		src/geoCalcs\
		src/sort\
	        src/warp\
#		src/common\
		src/kriging\
		src/rbf\
		src/engine\
		src/weightAndBalance\
		src/missionAnalysis\
		src/derivativeRoutines\
		src/tapenade/output\

#      ******************************************************************
#      *                                                                *
#      * General targets.                                               *
#      *                                                                *
#      ******************************************************************

default:
	@echo "Usage: make <arch>"
	@echo "Supported architectures: LINUX_INTEL"
	@echo "                         LINUX_INTEL_OPENMPI"
	@echo "                         LINUX_INTEL_OPENMPI_SCINET"
	@echo "                         LINUX_GFORTRAN"
	@echo "                         LINUX_GFORTRAN_OPENMPI"
	@echo "                         LINUX_INTEL_SCINET"
	@echo "                         BASALT"

all:	 default

clean:
	ln -sf Common_real.mk Common.mk

	@echo " Making clean ... "

	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo; \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make $@) || exit 1; \
		done
	rm -f *~ config.mk;
	rm -f lib/lib* mod/*.mod obj/*

#      ******************************************************************
#      *                                                                *
#      * The actual make. This is not a direct target, but is called    *
#      * from the architectures.                                        *
#      *                                                                *
#      ******************************************************************

module:
	@for subdir in $(SUBDIR_SRC) ; \
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

LINUX_INTEL:
	mkdir -p obj
	if [ ! -f "config/config.LINUX_INTEL.mk" ]; then cp "config/defaults/config.LINUX_INTEL.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL.mk config.mk
	ln -sf Common_real.mk Common.mk
	make module
	(cd src/f2py && make)

LINUX_INTEL_OPENMPI:
	mkdir -p obj
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL_OPENMPI.mk config.mk
	ln -sf Common_real.mk Common.mk
	make module
	(cd src/f2py && make)

LINUX_GFORTRAN:
	mkdir -p obj
	if [ ! -f "config/config.LINUX_GFORTRAN.mk" ]; then cp "config/defaults/config.LINUX_GFORTRAN.mk" ./config; fi
	ln -sf config/config.LINUX_GFORTRAN.mk config.mk
	ln -sf Common_real.mk Common.mk
	make module
	(cd src/f2py && make)

LINUX_GFORTRAN_OPENMPI:
	mkdir -p obj
	if [ ! -f "config/config.LINUX_GFORTRAN_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_GFORTRAN_OPENMPI.mk" ./config; fi
	ln -sf config/config.LINUX_GFORTRAN_OPENMPI.mk config.mk
	ln -sf Common_real.mk Common.mk
	make module
	(cd src/f2py && make)

LINUX_INTEL_SCINET:
	mkdir -p obj
	if [ ! -f "config/config.LINUX_INTEL_SCINET.mk" ]; then cp "config/defaults/config.LINUX_INTEL_SCINET.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL_SCINET.mk config.mk
	ln -sf Common_real.mk Common.mk
	make module
	(cd src/f2py && make)

LINUX_INTEL_OPENMPI_SCINET:
	mkdir -p obj
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI_SCINET.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI_SCINET.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL_OPENMPI_SCINET.mk config.mk
	ln -sf Common_real.mk Common.mk
	make module
	(cd src/f2py && make)


TAPENADE:
	(cd src/missionAnalysis && make tapenade)
