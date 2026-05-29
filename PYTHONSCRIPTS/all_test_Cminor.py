
#
# Copyright (C) 2025 Levin Rug (E-Mail: l.rug@lmu.de)
# See ./SRC/Cminor.f90 for the copyright notice
# See ./LICENSE for license information
# SPDX-License-Identifier: GPL-3.0
#

import numpy as np
import os
import sys
import subprocess
from pathlib import Path
from fcns_all_test_Cminor import ncdfcompare, Tee

# routine to run exemplary atmospheric and combustion mechanisms with Cminor to see if they produce the results they should

# runs to skip, enter numbers between >=1 and <=9 and the respective simulations will be skipped
iSkipRuns = np.array([])

# test profile:
#   full  -> run all cases (default)
#   quick -> run reduced smoke set for faster CI cycles
test_profile = os.getenv("CMINOR_TEST_PROFILE", "full").lower()
if test_profile not in ("full", "quick"):
    print("ERROR: unknown CMINOR_TEST_PROFILE =", test_profile)
    print("       use one of: full, quick")
    sys.exit(2)

# optional strict mode: fail process on large relative errors
fail_on_warning = os.getenv("CMINOR_FAIL_ON_WARNING", "0") == "1"

CMINOR_DIR = Path(__file__).parent.parent.resolve()
# run files to be tested
RUN_Files = CMINOR_DIR / np.array([ \
                       "RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP.run"                          \
                     , "RUN/TESTRUN/MCM/MCM.run"                                                \
                     , "RUN/TESTRUN/RACM_ML/RACM_ML.run"                                        \
                     , "RUN/TESTRUN/kreidenweis2003_parcel/kreidenweis2003_parcel.run"          \
                     , "RUN/TESTRUN/RACM+CAPRAM/RACM+C24.run"                                   \
                     , "RUN/TESTRUN/MCM+CAPRAM/MCM+CAPRAM.run"                                  \
                     , "RUN/TESTRUN/ERC_nHeptane/ERC_nHeptane.run"                              \
                     , "RUN/TESTRUN/LLNL_nHeptane/LLNL_nHeptane.run"                            \
                     , "RUN/TESTRUN/LLNL_MD/LLNL_MD.run"                                        \
                     ])

# netcdf files to be assumed the truth
reference = CMINOR_DIR / np.array([ \
                       "RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP_reference.nc"                 \
                     , "RUN/TESTRUN/MCM/MCM_reference.nc"                                       \
                     , "RUN/TESTRUN/RACM_ML/RACM_ML_reference.nc"                               \
                     , "RUN/TESTRUN/kreidenweis2003_parcel/kreidenweis2003_parcel_reference.nc" \
                     , "RUN/TESTRUN/RACM+CAPRAM/RACM+C24_reference.nc"                          \
                     , "RUN/TESTRUN/MCM+CAPRAM/MCM32+CAPRAM_reference.nc"                       \
                     , "RUN/TESTRUN/ERC_nHeptane/ERC_nHeptane_reference.nc"                     \
                     , "RUN/TESTRUN/LLNL_nHeptane/LLNL_nHeptane_reference.nc"                   \
                     , "RUN/TESTRUN/LLNL_MD/LLNL_MD_reference.nc"                               \
                     ])


# netcdf files that are generated (have to be the same as in the run files RUN_Files)
test_ncdf = CMINOR_DIR / np.array([ \
                       "RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP_test.nc"                      \
                     , "RUN/TESTRUN/MCM/MCM_test.nc"                                            \
                     , "RUN/TESTRUN/RACM_ML/RACM_ML_test.nc"                                    \
                     , "RUN/TESTRUN/kreidenweis2003_parcel/kreidenweis2003_parcel_test.nc"      \
                     , "RUN/TESTRUN/RACM+CAPRAM/RACM+C24_test.nc"                               \
                     , "RUN/TESTRUN/MCM+CAPRAM/MCM32+CAPRAM_test.nc"                            \
                     , "RUN/TESTRUN/ERC_nHeptane/ERC_nHeptane_test.nc"                          \
                     , "RUN/TESTRUN/LLNL_nHeptane/LLNL_nHeptane_test.nc"                        \
                     , "RUN/TESTRUN/LLNL_MD/LLNL_MD_test.nc"                                    \
                     ])

# mode specifier needed for ncdfcompare routine
modes = np.array([ \
                    "atm"  \
                  , "atm"  \
                  , "atm"  \
                  , "atm"  \
                  , "atm"  \
                  , "atm"  \
                  , "comb" \
                  , "comb" \
                  , "comb" \
                 ])


# threshold to consider a relative error large
rel_threshold = 1E-2

CMINOR_BIN = CMINOR_DIR / "Cminor"
CMINOR_OUTPUT_LOG = CMINOR_DIR / "RUN" / "TESTRUN" / "output_during_tests.txt"
CMINOR_LOG_TAIL_LINES = 40


def _print_cminor_log_tail():
    if not CMINOR_OUTPUT_LOG.is_file():
        print("  (no Cminor log at", CMINOR_OUTPUT_LOG, ")")
        return
    lines = CMINOR_OUTPUT_LOG.read_text(errors="replace").splitlines()
    tail = lines[-CMINOR_LOG_TAIL_LINES:]
    print("  --- last lines of", CMINOR_OUTPUT_LOG.name, "---")
    for line in tail:
        print("  ", line)
    print("  --- end log tail ---")


# file where Cminors output is written
outfile = open(CMINOR_OUTPUT_LOG, "w")
logfile = open(CMINOR_DIR / "RUN" / "TESTRUN" / "log_during_tests.txt", "w")
original_stdout = sys.stdout
sys.stdout = Tee(sys.stdout, logfile)
# global counter to count large relative errors and show at the end
warn_counter = 0
fail_counter = 0

if test_profile == "quick":
    selected_runs = np.array([0, 1, 2, 3, 7])  # 3 atm + parcel + 1 combustion
else:
    selected_runs = np.arange(len(RUN_Files))

for iRun in selected_runs:
    runfile = RUN_Files[iRun]
    testfile = test_ncdf[iRun]
    referencefile = reference[iRun]
    if iRun in iSkipRuns:
        continue

    label = runfile.name.replace(".run", "")
    print("\n####### LABEL: " + label + " #######\n")

    if not referencefile.is_file():
        print("ERROR: reference NetCDF missing:", referencefile)
        fail_counter += 1
        continue

    if not CMINOR_BIN.is_file():
        print("ERROR: Cminor executable missing:", CMINOR_BIN)
        print("       build with: make Cminor")
        fail_counter += 1
        continue

    print("Running " + str(runfile) + "...", end="\r")
    subprocess.run(["rm", "-f", str(testfile)], cwd=CMINOR_DIR, check=False)
    run_arg = runfile.relative_to(CMINOR_DIR)
    run_status = subprocess.run(
        [str(CMINOR_BIN), str(run_arg)],
        cwd=CMINOR_DIR,
        stdin=subprocess.DEVNULL,
        stdout=outfile,
        stderr=outfile,
    )
    outfile.flush()
    print("Running " + str(runfile) + "... Finished.")
    if run_status.returncode != 0:
        print("ERROR: Cminor run failed with return code", run_status.returncode)
        _print_cminor_log_tail()
        fail_counter += 1
        continue

    if not testfile.is_file():
        print("ERROR: expected output NetCDF missing:", testfile)
        _print_cminor_log_tail()
        fail_counter += 1
        continue

    print("\nChecking " + str(runfile) + " results...", end="\r")
    reldata, absdata = ncdfcompare(str(testfile), str(referencefile), modes[iRun])
    print("Checking " + str(runfile) + " results... Finished.")

    print("\nStatistics:\n")
    print("  Largest relative error: "+"%2.3f"%(reldata[0]*100)+"% for "+reldata[1]+" at values "+"%.5e"%reldata[2][0]+" and "+"%.5e"%reldata[2][1])
    print("  Largest absolute error: "+"%.3e"%absdata[0]+" for "+absdata[1]+" at values "+"%.5e"%absdata[2][0]+" and "+"%.5e"%absdata[2][1])
    if reldata[0]<rel_threshold:
        print("\nSuccess!")
    else:
        print("\nWARNING !!! Large relative error !")
        warn_counter += 1

if fail_counter > 0:
    print("\n  Finished testing Cminor. There were failures in "
          + str(fail_counter) + " run" + ("s" if fail_counter != 1 else "") + ".")
    print("  See", CMINOR_OUTPUT_LOG.name, "under RUN/TESTRUN/ for Cminor output.")
elif warn_counter > 0:
    print("\n  Finished testing Cminor. There were large relative errors in "
          + str(warn_counter) + " run" + ("s" if warn_counter != 1 else "") + "!")
    print("")
    print("  Please note that the indicator for large errors might already trigger if tolerances have been adjusted or similar.")
    print("  It might still be insignificant and is just a neccessary condition for a true failure.")
    print("  Please further inspect by looking in the respective NetCDF files to compare the numerical solutions.")
else:
    print("\nFinished testing Cminor. Everything went fine.")

print("")
print("")
sys.stdout = original_stdout
outfile.close()
logfile.close()

if fail_counter > 0:
    sys.exit(1)

if warn_counter > 0 and fail_on_warning:
    sys.exit(1)
