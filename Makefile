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
NETCDF_DIR=/opt/homebrew
BLAS_LIB_DIR=/opt/homebrew/lib
LAPACK_LIB_DIR=/opt/homebrew/lib

# --- TROPOS dusti paths
# NETCDF_DIR=/opt/tools/packages/gcc-7.2.0/netcdf-4.4.1.1
# BLAS_LIB_DIR=/opt/tools/packages/gcc-7.2.0/lapack-3.8.0/lib64
# LAPACK_LIB_DIR=/opt/tools/packages/gcc-7.2.0/lapack-3.8.0/lib64


#------------------------------------------------------------------------------
# Compiler settings
#------------------------------------------------------------------------------
FC=gfortran# Fortran compiler
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
FFLAGS_OPT += -I$(METHODS_DIR) -I.
FFLAGS_DBG += -I$(METHODS_DIR) -I.

#------------------------------------------------------------------------------
# Include paths for header files and modules
#------------------------------------------------------------------------------
INCLUDES_OPT = -I$(LIB_DIR) -I$(NETCDF_DIR)/include
INCLUDES_DBG = -I$(LIB_DBG_DIR) -I$(NETCDF_DIR)/include

#------------------------------------------------------------------------------
# External libraries required for linking
# Add the library paths to the runtime library search path using the -Wl,-rpath flag during linking
#------------------------------------------------------------------------------
LIBS = -Wl,-rpath,$(NETCDF_DIR)/lib,-rpath,$(BLAS_LIB_DIR),-rpath,$(LAPACK_LIB_DIR) \
       -lcurl -L$(NETCDF_DIR)/lib -L$(BLAS_LIB_DIR) -L$(LAPACK_LIB_DIR) \
       -lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -llapack -lblas

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

#------------------------------------------------------------------------------
# Object files - automatically generated from source files
#------------------------------------------------------------------------------
OBJS_OPT = $(patsubst $(SRC_DIR)/%.f90,$(LIB_DIR)/%.o,$(SRCS))
OBJS_DBG = $(patsubst $(SRC_DIR)/%.f90,$(LIB_DBG_DIR)/%.o,$(SRCS))

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

#------------------------------------------------------------------------------
# Test and cleanup targets
#------------------------------------------------------------------------------
# Run test suite
test: Cminor
	./Cminor RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP.run
	./Cminor RUN/TESTRUN/RACM_ML/RACM_ML.run
	./Cminor RUN/TESTRUN/MCM/MCM.run
	./Cminor RUN/TESTRUN/kreidenweis2003_parcel/kreidenweis2003_parcel.run
	./Cminor RUN/TESTRUN/MCM+CAPRAM/MCM+CAPRAM.run
	./Cminor RUN/TESTRUN/ERC_nheptane/ERC_nheptane.run
	./Cminor RUN/TESTRUN/LLNL_nHeptane/LLNL_nHeptane.run
	./Cminor RUN/TESTRUN/LLNL_MD/LLNL_MD.run

# Clean build artifacts - only within project directories
clean:
	rm -f $(LIB_DIR)/*.o $(LIB_DIR)/*.mod $(LIB_DIR)/*.a
	rm -f $(LIB_DBG_DIR)/*.o $(LIB_DBG_DIR)/*.mod $(LIB_DBG_DIR)/*.a
	rm -f Cminor Cminor_dbg

# Declare phony targets (targets that don't create files)
.PHONY: all clean test
