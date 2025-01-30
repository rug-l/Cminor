# Cminor
We present the Chemical Mechanism Integrator (Cminor), a fully modularized modern Fortran software package for the computational simulation of skeletal and detailed chemical kinetic systems derived from atmospheric and combustion chemistry.
  Cminor aims for the efficient simulation of complex chemical mechanisms by using various mathematical techniques. These are tailored to systems of ordinary differential equations (ODEs), having the specific structure arising from chemical reaction systems. Additionally, a high-speed mechanism parser allows the user to interchange reactions or their parameters in an ASCII format text file and immediately start a new simulation without recompiling, enabling fast and numerous simulations.
  Cminor's solver technique is based on Rosenbrock methods. Different measures of local errors and an analytical Jacobian matrix approach are implemented, where efficiency is obtained by exploiting the sparsity structure of the Jacobian. 
  
Cminor can be run in one of three configurations: 
  - A box-model framework for either pure gas-phase mechanisms or a multi-modal aerosol distribution dissolved in mono-dispersed cloud droplets.
  - A rising adiabatic parcel, in which the activation of multi-modal aerosols is represented by solving the droplet condensation equation.
  - A constant volume environment, where thermodynamic properties are evaluated by polynomial functions of temperature according to the standards of the Chemkin thermodynamic data base.


## Installation

1. Configure the required paths in the Makefile:
   - `NETCDF_DIR`: Path to NetCDF installation (must contain 'include' and 'lib')
   - `BLAS_LIB_DIR`: Path to BLAS libraries
   - `LAPACK_LIB_DIR`: Path to LAPACK libraries

2. If not using gfortran, adjust the compiler settings in the Makefile:
   - `FC`: Compiler name
   - `FFLAGS_OPT`: Optimization flags
   - `FFLAGS_DBG`: Debug flags

3. Build the program:
   ```bash
   make Cminor     # for optimized version
   # or
   make Cminor_dbg # for debugging version
   # or
   make
   ```

4. Run the test suite:
   ```bash
   make test
   ```

5. Execute a specific mechanism:
   ```bash
   ./Cminor RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP.run
   ```