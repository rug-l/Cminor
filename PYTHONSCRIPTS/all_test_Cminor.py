
#
# Copyright (C) 2025 Levin Rug (E-Mail: l.rug@lmu.de)
# See ./SRC/Cminor.f90 for the copyright notice
# See ./LICENSE for license information
# SPDX-License-Identifier: GPL-3.0
#

from fcns_all_test_Cminor import *

# routine to run exemplary atmospheric and combustion mechanisms with Cminor to see if they produce the results they should


# run files to be tested
RUN_Files = np.array([ \
                       "RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP.run"                          \
                     , "RUN/TESTRUN/MCM/MCM.run"                                                \
                     , "RUN/TESTRUN/RACM_ML/RACM_ML.run"                                        \
                     , "RUN/TESTRUN/kreidenweis2003_parcel/kreidenweis2003_parcel.run"          \
                     , "RUN/TESTRUN/RACM+CAPRAM/RACM+C24.run"                                   \
                     , "RUN/TESTRUN/MCM+CAPRAM/MCM+CAPRAM.run"                                  \
                     , "RUN/TESTRUN/ERC_nHeptane/ERC_nheptane.run"                              \
                     , "RUN/TESTRUN/LLNL_nHeptane/LLNL_nHeptane.run"                            \
                     , "RUN/TESTRUN/LLNL_MD/LLNL_MD.run"                                        \
                     ])


# netcdf files to be assumed the truth
reference = np.array([ \
                       "RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP_reference.nc"                 \
                     , "RUN/TESTRUN/MCM/MCM_reference.nc"                                       \
                     , "RUN/TESTRUN/RACM_ML/RACM_ML_reference.nc"                               \
                     , "RUN/TESTRUN/kreidenweis2003_parcel/kreidenweis2003_parcel_reference.nc" \
                     , "RUN/TESTRUN/RACM+CAPRAM/RACM+C24_reference.nc"                          \
                     , "RUN/TESTRUN/MCM+CAPRAM/MCM32+CAPRAM_reference.nc"                       \
                     , "RUN/TESTRUN/ERC_nHeptane/ERC_nheptane_reference.nc"                     \
                     , "RUN/TESTRUN/LLNL_nHeptane/LLNL_nHeptane_reference.nc"                   \
                     , "RUN/TESTRUN/LLNL_MD/LLNL_MD_reference.nc"                               \
                     ])

# netcdf files that are generated (have to be the same as in the run files RUN_Files)
test_ncdf = np.array([ \
                       "RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP_test.nc"                      \
                     , "RUN/TESTRUN/MCM/MCM_test.nc"                                            \
                     , "RUN/TESTRUN/RACM_ML/RACM_ML_test.nc"                                    \
                     , "RUN/TESTRUN/kreidenweis2003_parcel/kreidenweis2003_parcel_test.nc"      \
                     , "RUN/TESTRUN/RACM+CAPRAM/RACM+C24_test.nc"                               \
                     , "RUN/TESTRUN/MCM+CAPRAM/MCM32+CAPRAM_test.nc"                            \
                     , "RUN/TESTRUN/ERC_nHeptane/ERC_nheptane_test.nc"                          \
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

# file where Cminors output is written
outfile = open("RUN/TESTRUN/output_during_tests.txt", "w")
logfile = open("RUN/TESTRUN/log_during_tests.txt", "w")
original_stdout = sys.stdout
sys.stdout = Tee(sys.stdout, logfile)
# global counter to count large relative errors and show at the end
warn_counter = 0
for iRun, runfile in enumerate(RUN_Files):
    print("\n####### LABEL: "+runfile[runfile.rindex("/")+1:runfile.rindex(".run")]+" #######\n")

    print("Running "+runfile+"...", end="\r")
    subprocess.run(["rm", "-f", test_ncdf[iRun]])
    subprocess.run(["./Cminor", runfile], stdout=outfile)
    print("Running "+runfile+"... Finished.")

    print("\nChecking "+runfile+" results...", end="\r")
    reldata, absdata = ncdfcompare(test_ncdf[iRun], reference[iRun], modes[iRun])
    print("Checking "+runfile+" results... Finished.")

    print("\nStatistics:\n")
    print("  Largest relative error: "+"%2.3f"%(reldata[0]*100)+"% for "+reldata[1]+" at values "+"%.5e"%reldata[2][0]+" and "+"%.5e"%reldata[2][1])
    print("  Largest absolute error: "+"%.3e"%absdata[0]+" for "+absdata[1]+" at values "+"%.5e"%absdata[2][0]+" and "+"%.5e"%absdata[2][1])
    if reldata[0]<rel_threshold:
        print("\nSuccess!")
    else:
        print("\nWARNING !!! Large relative error !")
        warn_counter += 1

if warn_counter>0:
    print("\nFinished testing Cminor. There were large relative errors in "+str(warn_counter)+" run"+(warn_counter!=1)*"s"+"!")
else:
    print("\nFinished testing Cminor. Everything went fine.")
print("")
print("")
sys.stdout = original_stdout
outfile.close()
logfile.close()
