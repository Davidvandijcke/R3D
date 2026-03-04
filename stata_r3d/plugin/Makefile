# Makefile for STATA R3D plugin
# This builds a shared library that can be loaded by STATA

# Compiler settings
CC = gcc
FC = gfortran
CFLAGS = -O3 -fPIC -Wall
FFLAGS = -O3 -fPIC -Wall
LDFLAGS =

# STATA plugin header location
# You may need to adjust this path to match your STATA installation
STATA_INCLUDE = /usr/local/stata/ado/base/

# Platform-specific settings
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
    PLUGIN_EXT = plugin
    LDFLAGS += -shared -llapack -lblas -lgfortran
    CFLAGS += -I$(STATA_INCLUDE) -DSYSTEM=LINUX
endif
ifeq ($(UNAME_S),Darwin)
    PLUGIN_EXT = plugin
    # Use dynamiclib like the working quadprog plugin
    GFORTRAN_STATIC := $(shell gfortran -print-file-name=libgfortran.a)
    QUADMATH_STATIC := $(shell gfortran -print-file-name=libquadmath.a)
    LDFLAGS += -dynamiclib -undefined dynamic_lookup -framework Accelerate
    LDFLAGS += $(GFORTRAN_STATIC) $(QUADMATH_STATIC)
    CFLAGS += -I$(STATA_INCLUDE) -DSYSTEM=APPLEMAC
endif
ifeq ($(OS),Windows_NT)
    PLUGIN_EXT = dll
    LDFLAGS += -shared -llapack -lblas -lgfortran
    CFLAGS += -I$(STATA_INCLUDE) -DSYSTEM=STWIN32 -fdeclspec
endif

# Targets
TARGET = r3d_plugin.$(PLUGIN_EXT)
F_SOURCES = r3d_fortran.f90
C_SOURCES = r3d_plugin.c
F_OBJECTS = $(F_SOURCES:.f90=.o)
C_OBJECTS = $(C_SOURCES:.c=.o)
OBJECTS = $(F_OBJECTS) $(C_OBJECTS)

# Build rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

install: $(TARGET)
	@echo "Copy $(TARGET) to your STATA ado-path directory"
	@echo "Example: cp $(TARGET) ~/ado/personal/"
	cp $(TARGET) ../ado/

.PHONY: all clean install