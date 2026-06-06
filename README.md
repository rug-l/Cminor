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

4. Install Python dependencies for regression tests:
   ```bash
   pip install -r requirements.txt
   ```

5. Run the test suites:
   ```bash
   make test             # smoke tests (legacy behavior)
   make test-regression  # NetCDF regression checks
   make test-diagnosis   # NetCDF regression + diagnostic PNGs
   ```

   Optional environment flags for regression runner:
   ```bash
   CMINOR_TEST_PROFILE=quick make test-regression      # faster subset
   CMINOR_FAIL_ON_WARNING=1 make test-regression       # fail on large relative error warnings
   ```

   Example: generate a diagnostic PNG for one mechanism:
   ```bash
   python PYTHONSCRIPTS/plot_test_diagnosis.py \
     --test RUN/TESTRUN/LLNL_nHeptane/LLNL_nHeptane_test.nc \
     --ref RUN/TESTRUN/LLNL_nHeptane/LLNL_nHeptane_reference.nc \
     --mode comb \
     --label LLNL_nHeptane \
     --out RUN/TESTRUN/diagnostics/LLNL_nHeptane_diagnostic.png
   ```

   Diagnostic figures are written to:
   ```bash
   RUN/TESTRUN/diagnostics/
   ```

6. Execute a specific mechanism:
   ```bash
   ./Cminor RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP.run
   ```

## Continuous integration

GitHub Actions (`.github/workflows/ci.yml`) builds `Cminor` on Ubuntu and macOS, then runs `make test-regression` with `CMINOR_TEST_PROFILE=quick`.

Override NetCDF paths without editing the Makefile (as CI does):

```bash
export NETCDF_DIR=/usr
export NETCDF_C_DIR=/usr
export HDF5_DIR=/usr
make Cminor
```

On macOS runners, Homebrew `gcc` supplies `gfortran-*`; CI sets `FC` accordingly. The project still uses the Makefile as the canonical build; [fpm](https://github.com/fortran-lang/fpm) is optional for local workflows (`brew install fpm`) but not required for CI.