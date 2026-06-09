################################################################################
# Cminor - Chemical Kinetics Solver
################################################################################
# Copyright (C) 2025 Levin Rug, Willi Schimmel
# Contact: l.rug@lmu.de
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# See ./SRC/Cminor.f90 for the full copyright notice
# See ./LICENSE for complete license information
# SPDX-License-Identifier: GPL-3.0
################################################################################
#
# REQUIRED USER CONFIGURATION:
# 1. Set the following paths according to your system:
#    - NETCDF_DIR: Path to NetCDF installation (must contain 'include' and 'lib')
#    - BLAS_LIB_DIR: Path to BLAS libraries
#    - LAPACK_LIB_DIR: Path to LAPACK libraries
#
# 2. Compiler settings (if needed):
#    - FC: Change compiler if not using gfortran (default)
#    - FFLAGS_OPT: Adjust optimization flags if using different compiler
#    - FFLAGS_DBG: Adjust debug flags if using different compiler
#
# Example configurations for different systems are provided below.
################################################################################

#------------------------------------------------------------------------------
# System-specific paths - adjust these for your system
#------------------------------------------------------------------------------
# Provide the paths to the netcdf installation directory (includes subdirectories 
# `include` and `lib`) and paths to the blas and lapack `lib` (libraries) 

# --- Mac paths
# BLAS and LAPACK library directories
# On macOS, use the Accelerate framework for BLAS+LAPACK (see LIBS below);
# BLAS_LIB_DIR/LAPACK_LIB_DIR kept for compatibility with the Linux branch.
NETCDF_DIR ?= /opt/homebrew/opt/netcdf-fortran
NETCDF_C_DIR ?= /opt/homebrew/opt/netcdf
HDF5_DIR ?= /opt/homebrew/opt/hdf5
BLAS_LIB_DIR ?= /opt/homebrew/lib
LAPACK_LIB_DIR ?= /opt/homebrew/lib

# --- TROPOS dusti paths
# NETCDF_DIR=/opt/tools/packages/gcc-7.2.0/netcdf-4.4.1.1
# BLAS_LIB_DIR=/opt/tools/packages/gcc-7.2.0/lapack-3.8.0/lib64
# LAPACK_LIB_DIR=/opt/tools/packages/gcc-7.2.0/lapack-3.8.0/lib64


#------------------------------------------------------------------------------
# Compiler settings
#------------------------------------------------------------------------------
# Default gfortran; overrides conda/shell FC=f77. CI/mac: make Cminor FC="${FC}"
FC = gfortran# Fortran compiler
PYTHON ?= python
# Free-form Fortran flags to handle line length
FFLAGS_FREE = -ffree-form -ffixed-line-length-none -ffree-line-length-none

#------------------------------------------------------------------------------
# Directory structure
#------------------------------------------------------------------------------
SRC_DIR=SRC# Source code directory
LIB_DIR=LIB# Optimized build library directory
LIB_DBG_DIR=LIB_D# Debug build library directory
METHODS_DIR=METHODS# Methods directory

#------------------------------------------------------------------------------
# Compiler flags
#------------------------------------------------------------------------------
# Optimization flags
FFLAGS_OPT = -J$(LIB_DIR) -O3  # -J specifies where to put .mod files, -O3 for optimization

# Debug flags with extensive error checking
FFLAGS_DBG = -J$(LIB_DBG_DIR) -g -C -O0 -Warray-bounds -Wextra -fbacktrace \
             -ffpe-trap=zero -fimplicit-none -fcheck=all -Wall \
             -Wno-uninitialized -Wno-compare-reals

# Add include path to compiler flags
FFLAGS_OPT += -I$(METHODS_DIR) -I. -cpp
FFLAGS_DBG += -I$(METHODS_DIR) -I. -cpp

#------------------------------------------------------------------------------
# Include paths for header files and modules
# Linux: nf-config supplies Fortran module path (not always under include/)
#------------------------------------------------------------------------------
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
  NETCDF_FFLAGS := $(shell nf-config --fflags 2>/dev/null)
  INCLUDES_OPT = -I$(LIB_DIR) -I$(NETCDF_DIR)/include $(NETCDF_FFLAGS)
  INCLUDES_DBG = -I$(LIB_DBG_DIR) -I$(NETCDF_DIR)/include $(NETCDF_FFLAGS)
else
  INCLUDES_OPT = -I$(LIB_DIR) -I$(NETCDF_DIR)/include
  INCLUDES_DBG = -I$(LIB_DBG_DIR) -I$(NETCDF_DIR)/include
endif

#------------------------------------------------------------------------------
# External libraries required for linking
# Add the library paths to the runtime library search path using the -Wl,-rpath flag during linking
# macOS: use Accelerate for BLAS+LAPACK (no separate -lblas -llapack)
# Linux (CI/apt): OpenBLAS + LAPACK from distro packages
#------------------------------------------------------------------------------
ifeq ($(UNAME_S),Linux)
# Debian/Ubuntu: HDF5 libs live under hdf5/serial (libhdf5_serial, not libhdf5)
HDF5_LIBDIR := /usr/lib/x86_64-linux-gnu/hdf5/serial
LIBS = -Wl,-rpath,$(NETCDF_DIR)/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu \
       -Wl,-rpath,$(HDF5_LIBDIR) \
       -lcurl -L$(NETCDF_DIR)/lib -L/usr/lib/x86_64-linux-gnu -L$(HDF5_LIBDIR) \
       -lnetcdff -lnetcdf -lhdf5_serial_hl -lhdf5_serial -llapack -lopenblas
else
LIBS = -Wl,-rpath,$(NETCDF_DIR)/lib,-rpath,$(NETCDF_C_DIR)/lib,-rpath,$(HDF5_DIR)/lib \
       -lcurl -L$(NETCDF_DIR)/lib -L$(NETCDF_C_DIR)/lib -L$(HDF5_DIR)/lib \
       -lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -framework Accelerate
endif

#------------------------------------------------------------------------------
# Optional ISSA reduction (private submodule / local scaffold)
#------------------------------------------------------------------------------
ISSA_DIR = EXTERNALS/cminor-issa/SRC
ISSA_SRCS = $(ISSA_DIR)/issa_flux_io.f90 \
            $(ISSA_DIR)/issa_graph_tools.f90 \
            $(ISSA_DIR)/issa_reduce.f90
ISSA_PRESENT := $(wildcard $(ISSA_DIR)/issa_reduce.f90)
ifneq ($(ISSA_PRESENT),)
  FFLAGS_OPT += -DISSA
  FFLAGS_DBG += -DISSA
  ISSA_OBJS_OPT = $(LIB_DIR)/issa_flux_io.o $(LIB_DIR)/issa_graph_tools.o $(LIB_DIR)/issa_reduce.o
  ISSA_OBJS_DBG = $(LIB_DBG_DIR)/issa_flux_io.o $(LIB_DBG_DIR)/issa_graph_tools.o $(LIB_DBG_DIR)/issa_reduce.o
else
  ISSA_OBJS_OPT =
  ISSA_OBJS_DBG =
endif

#------------------------------------------------------------------------------
# Source files - all Fortran 90 files that make up the project
#------------------------------------------------------------------------------
SRCS = $(addprefix $(SRC_DIR)/, \
	UniRnk_Mod.f90 LexicalStringSort.f90 Kind_Mod.f90 Control_Mod.f90 \
	Reac_Mod.f90 Meteo_Mod.f90 InitRoutines_Mod.f90 NetCDF_Mod.f90 \
	Sparse_Mod.f90 String_Mod.f90 HashStr_Mod.f90 InputTool_Mod.f90 \
	ChemSys_Mod.f90 IO_Mod.f90 ChemKinInput_Mod.f90 fp_parameters.f90 \
	fparser.f90 Rates_Mod.f90 Rosenbrock_Mod.f90 Integration_Mod.f90 \
	Cminor.f90)

BASE_OBJS_OPT = $(patsubst $(SRC_DIR)/%.f90,$(LIB_DIR)/%.o,$(SRCS))
BASE_OBJS_DBG = $(patsubst $(SRC_DIR)/%.f90,$(LIB_DBG_DIR)/%.o,$(SRCS))

# Insert ISSA objects after Integration_Mod.o, before Cminor.o
OBJS_OPT = $(filter-out $(LIB_DIR)/Cminor.o,$(BASE_OBJS_OPT)) $(ISSA_OBJS_OPT) $(LIB_DIR)/Cminor.o
OBJS_DBG = $(filter-out $(LIB_DBG_DIR)/Cminor.o,$(BASE_OBJS_DBG)) $(ISSA_OBJS_DBG) $(LIB_DBG_DIR)/Cminor.o

#------------------------------------------------------------------------------
# Build targets
#------------------------------------------------------------------------------
# Default target: build both optimized and debug versions
all: Cminor Cminor_dbg

# Create necessary directories
$(LIB_DIR):
	mkdir -p $(LIB_DIR)

$(LIB_DBG_DIR):
	mkdir -p $(LIB_DBG_DIR)

# Optimized version
Cminor: $(OBJS_OPT)
	$(FC) $(FFLAGS_OPT) $(INCLUDES_OPT) -o $@ $^ $(LIBS)

# Debug version
Cminor_dbg: $(OBJS_DBG)
	$(FC) $(FFLAGS_DBG) $(INCLUDES_DBG) -o $@ $^ $(LIBS)

# Pattern rules for compilation
%.o: %.f90
	$(FC) $(FFLAGS_FREE) $(FFLAGS_OPT) $(INCLUDES_OPT) -c $< -o $@

$(LIB_DIR)/%.o: $(SRC_DIR)/%.f90 | $(LIB_DIR)
	$(FC) $(FFLAGS_FREE) $(FFLAGS_OPT) $(INCLUDES_OPT) -c $< -o $@

$(LIB_DBG_DIR)/%.o: $(SRC_DIR)/%.f90 | $(LIB_DBG_DIR)
	$(FC) $(FFLAGS_FREE) $(FFLAGS_DBG) $(INCLUDES_DBG) -c $< -o $@

ifneq ($(ISSA_PRESENT),)
$(LIB_DIR)/issa_flux_io.o: $(ISSA_DIR)/issa_flux_io.f90 | $(LIB_DIR)
	$(FC) $(FFLAGS_FREE) $(FFLAGS_OPT) $(INCLUDES_OPT) -c $< -o $@

$(LIB_DIR)/issa_graph_tools.o: $(ISSA_DIR)/issa_graph_tools.f90 $(LIB_DIR)/issa_flux_io.o | $(LIB_DIR)
	$(FC) $(FFLAGS_FREE) $(FFLAGS_OPT) $(INCLUDES_OPT) -c $< -o $@

$(LIB_DIR)/issa_reduce.o: $(ISSA_DIR)/issa_reduce.f90 $(LIB_DIR)/issa_graph_tools.o | $(LIB_DIR)
	$(FC) $(FFLAGS_FREE) $(FFLAGS_OPT) $(INCLUDES_OPT) -c $< -o $@

$(LIB_DBG_DIR)/issa_flux_io.o: $(ISSA_DIR)/issa_flux_io.f90 | $(LIB_DBG_DIR)
	$(FC) $(FFLAGS_FREE) $(FFLAGS_DBG) $(INCLUDES_DBG) -c $< -o $@

$(LIB_DBG_DIR)/issa_graph_tools.o: $(ISSA_DIR)/issa_graph_tools.f90 $(LIB_DBG_DIR)/issa_flux_io.o | $(LIB_DBG_DIR)
	$(FC) $(FFLAGS_FREE) $(FFLAGS_DBG) $(INCLUDES_DBG) -c $< -o $@

$(LIB_DBG_DIR)/issa_reduce.o: $(ISSA_DIR)/issa_reduce.f90 $(LIB_DBG_DIR)/issa_graph_tools.o | $(LIB_DBG_DIR)
	$(FC) $(FFLAGS_FREE) $(FFLAGS_DBG) $(INCLUDES_DBG) -c $< -o $@

$(LIB_DIR)/issa_reduce_main.o: EXTERNALS/cminor-issa/issa_reduce_main.f90 $(LIB_DIR)/issa_reduce.o | $(LIB_DIR)
	$(FC) $(FFLAGS_FREE) $(FFLAGS_OPT) $(INCLUDES_OPT) -c $< -o $@

ISSA_REDUCE_SKIP = $(LIB_DIR)/Cminor.o $(LIB_DIR)/Integration_Mod.o $(LIB_DIR)/Rosenbrock_Mod.o \
                   $(LIB_DIR)/NetCDF_Mod.o $(LIB_DIR)/Rates_Mod.o $(LIB_DIR)/ChemKinInput_Mod.o
ISSA_REDUCE_OBJS = $(filter-out $(ISSA_REDUCE_SKIP),$(BASE_OBJS_OPT)) $(ISSA_OBJS_OPT) $(LIB_DIR)/issa_reduce_main.o

issa_reduce: $(ISSA_REDUCE_OBJS)
	$(FC) $(FFLAGS_OPT) $(INCLUDES_OPT) -o $@ $^ $(LIBS)
endif

#------------------------------------------------------------------------------
# Test and cleanup targets
#------------------------------------------------------------------------------
test: Cminor
	./Cminor RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP.run
	./Cminor RUN/TESTRUN/RACM_ML/RACM_ML.run
	./Cminor RUN/TESTRUN/MCM/MCM.run
	./Cminor RUN/TESTRUN/kreidenweis2003_parcel/kreidenweis2003_parcel.run
	./Cminor RUN/TESTRUN/RACM+CAPRAM/RACM+C24.run
	./Cminor RUN/TESTRUN/MCM+CAPRAM/MCM+CAPRAM.run
	./Cminor RUN/TESTRUN/ERC_nHeptane/ERC_nHeptane.run
	./Cminor RUN/TESTRUN/LLNL_nHeptane/LLNL_nHeptane.run
	./Cminor RUN/TESTRUN/LLNL_MD/LLNL_MD.run

# Run full regression checks against NetCDF references
test-regression: Cminor
	$(PYTHON) PYTHONSCRIPTS/all_test_Cminor.py

# Run diagnosis checks, create plots of the test and reference NetCDF files against NetCDF references
test-diagnosis: test-regression
	mkdir -p RUN/TESTRUN/diagnostics
	$(PYTHON) PYTHONSCRIPTS/plot_test_diagnosis.py \
		--test RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP_test.nc \
		--ref RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP_reference.nc \
		--mode atm \
		--label SmallStratoKPP \
		--out PYTHONSCRIPTS/Figures/test_diagnosis/SmallStratoKPP_diagnostic.png
	$(PYTHON) PYTHONSCRIPTS/plot_test_diagnosis.py \
		--test RUN/TESTRUN/MCM/MCM_test.nc \
		--ref RUN/TESTRUN/MCM/MCM_reference.nc \
		--mode atm \
		--label MCM \
		--out PYTHONSCRIPTS/Figures/test_diagnosis/MCM_diagnostic.png
	$(PYTHON) PYTHONSCRIPTS/plot_test_diagnosis.py \
		--test RUN/TESTRUN/LLNL_MD/LLNL_MD_test.nc \
		--ref RUN/TESTRUN/LLNL_MD/LLNL_MD_reference.nc \
		--mode comb \
		--label LLNL_MD \
		--out PYTHONSCRIPTS/Figures/test_diagnosis/LLNL_MD_diagnostic.png

# Clean build artifacts - only within project directories
clean:
	rm -f $(LIB_DIR)/*.o $(LIB_DIR)/*.mod $(LIB_DIR)/*.a
	rm -f $(LIB_DBG_DIR)/*.o $(LIB_DBG_DIR)/*.mod $(LIB_DBG_DIR)/*.a
	rm -f Cminor Cminor_dbg issa_reduce

# Compare full vs ISSA-reduced NetCDF (requires both NC files from prior Cminor runs)
issa-verify:
	mkdir -p RUN/TESTRUN/diagnostics
	$(PYTHON) PYTHONSCRIPTS/compare_issa_verification.py \
		--full-run    RUN/TESTRUN/MCM+CAPRAM/MCM+CAPRAM_issa.run \
		--reduced-run RUN/TESTRUN/MCM+CAPRAM/MCM+CAPRAM_issa_reduced.run \
		--label       MCM32_CAPRAM40_issa \
		--out-dir     RUN/TESTRUN/diagnostics

# Declare phony targets (targets that don't create files)
.PHONY: all clean test test-regression test-diagnosis issa_reduce issa-verify
