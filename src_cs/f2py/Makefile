MAIN_DIR = ../..
COMMON_FILE = $(MAIN_DIR)/Common.mk
RULES_FILE  = $(MAIN_DIR)/Rules.mk
include ${COMMON_FILE}
include ${RULES_FILE}

LINKER_ALL_FLAGS = $(FF90_FLAGS) -L$(LIBDIR) -lwarpustruct_cs $(PETSC_LINKER_FLAGS) $(LINKER_FLAGS) $(CGNS_LINKER_FLAGS)

PYTHON_OBJECTS = fortranobject.o \
		 warpustruct_csmodule.o \
		 warpustruct_cs-f2pywrappers2.o 

default: all
all:

# Generate Python inlude directory
	 $(eval PYTHON_INCLUDES = $(shell $(PYTHON-CONFIG) --includes))
	 @echo "#------------------------------------------------------#"
	 @echo Python Inclue Flags $(PYTHON_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate Numpy inlude directory
	$(eval NUMPY_INCLUDES = $(shell $(PYTHON) -c 'import numpy; print numpy.get_include()'))
	 @echo "#------------------------------------------------------#"
	 @echo Numpy Include Directory: $(NUMPY_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate f2py root directory
	$(eval F2PY_ROOT = $(shell $(PYTHON) get_f2py.py))
	 @echo "#------------------------------------------------------#"
	 @echo f2py root directory: $(F2PY_ROOT)
	 @echo "#------------------------------------------------------#"

	f2py warpustruct.pyf.autogen
	@echo " "
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-I$(F2PY_ROOT)/src -c warpustruct_csmodule.c

	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) -c \
	$(F2PY_ROOT)/src/fortranobject.c -o fortranobject.o

        # Compiled f2py generated wrapper file
        # $(FF90) $(FF90_ALL_FLAGS) -I$(MAIN_DIR)/mod -c warpustruct-f2pywrappers.f
	$(FF90) $(FF90_ALL_FLAGS) -I$(MAIN_DIR)/mod -c warpustruct_cs-f2pywrappers2.f90 

	$(FF90) -shared $(PYTHON_OBJECTS) $(LINKER_ALL_FLAGS) -o warpustruct_cs.so
	$(PYTHON) importTest.py
	mv warpustruct_cs.so ../../python

# Clean up temporary files
	-rm *.o
	-rm *.c
	-rm .f2py_f2cmap
clean:
	@echo "        Making clean in f2py... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)
