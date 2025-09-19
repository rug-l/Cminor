
#
# Copyright (C) 2025 Levin Rug (E-Mail: l.rug@lmu.de)
# See ./SRC/Cminor.f90 for the copyright notice
# See ./LICENSE for license information
# SPDX-License-Identifier: GPL-3.0
#

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import netCDF4 as nc
import sys

custompalette_colorblind=['#ebac23','#b80058','#008cf9','#00bbad', '#cc38b8', '#999999', '#636363', '#000000']
sns.set_palette(custompalette_colorblind)
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
plt.rc('font',**{'family':'serif','serif':['Times']})

print("")
print("  Starting NetCDF overview plots.\n")

mechanisms = []

mechanisms = [ *mechanisms, "ERC_nHeptane"]
mechanisms = [ *mechanisms, "LLNL_MD"]
mechanisms = [ *mechanisms, "LLNL_nHeptane"]
mechanisms = [ *mechanisms, "SmallStratoKPP"]
mechanisms = [ *mechanisms, "RACM+C24"]
mechanisms = [ *mechanisms, "MCM_CAPRAM"]
#mechanisms = [ *mechanisms, "MCM"]

Cminor_dir = "/home/l/L.Rug/CminorKPP/Cminor/"

sublabel=""
for mechanism in mechanisms:
    LABEL = mechanism
    if mechanism=="MCM":
        name = "MCMv3.2"
        fn=Cminor_dir+'RUN/TESTRUN/MCM/MCM_test.nc'
        atm = True
        plt_vars = np.array([["O1D"], ["NO", "NO2", "NO3"], ["NO"], ["NO2"], ["NO3"], ["OH"], ["O3"]], dtype=object)
        plt_vars = np.array([["O1D"], ["NO"], ["NO2"], ["O3"], ["H2O2"]], dtype=object)
        spcnames = np.array(["O($^1$D)", "NO$_x$", "NO", "NO$_2$", "NO$_3$", "OH", "O$_3$"])
        spcnames = np.array(["O($^1$D)", "NO", "NO$_2$", "O$_3$", "H$_2$O$_2$"])

    elif mechanism=="SmallStratoKPP":
        name = "Chapman-like"
        sublabel = "(a)"
        fn=Cminor_dir+'RUN/TESTRUN/SmallStratoKPP/SmallStratoKPP_test.nc'
        atm = True
        plt_vars = np.array([["O1D"], ["O"], ["O3"], ["NO"], ["NO2"]], dtype=object)
        spcnames = np.array(["O($^1$D)", "O", "O$_3$", "NO", "NO$_2$"])

    elif mechanism=="RACM+C24":
        name = "RACM+CAPRAMv2.4"
        sublabel = "(b)"
        fn=Cminor_dir+'RUN/TESTRUN/RACM+CAPRAM/RACM+C24_test.nc'
        atm = True
        plt_vars = np.array([["aNO_1_m3", "aNO2_1_m3", "aNO3_1_m3"], ["NO", "NO2", "NO3"], ["O3"], ["aO3_1_m3"], ["HO"], ["aHO_1_m3"]], dtype=object)
        spcnames = np.array(["aNO$_x$", "NO$_x$", "O$_3$", "aO$_3$", "OH", "aOH"])

    elif mechanism=="MCM_CAPRAM":
        name = "MCMv3.2+CAPRAMv4.0"+r"$\alpha$"
        sublabel = "(c)"
        fn=Cminor_dir+'RUN/TESTRUN/MCM+CAPRAM/MCM32+CAPRAM_test.nc'
        atm = True
        plt_vars = np.array([["aNO_1_m3", "aNO2_1_m3", "aNO3_1_m3"], ["NO", "NO2", "NO3"], ["O3"], ["aO3_1_m3"], ["OH"], ["aHO_1_m3"]], dtype=object)
        spcnames = np.array(["aNO$_x$", "NO$_x$", "O$_3$", "aO$_3$", "OH", "aOH"])

    elif mechanism=="ERC_nHeptane":
        name = "ERC $n$-heptane"
        sublabel = "(a)"
        fn=Cminor_dir+'RUN/TESTRUN/ERC_nHeptane/ERC_nHeptane_test.nc'
        atm = False
        plt_vars = np.array([["oh"], ["co2"], ["co"], ["h2o"], ["ho2"], ["nc7h16"]], dtype=object)
        spcnames = np.array(["OH", "CO$_2$", "CO", "H$_2$O", "HO$_2$", "C$_7$H$_{16}$"])

    elif mechanism=="LLNL_nHeptane":
        name = "LLNL $n$-heptane"
        sublabel = "(b)"
        fn=Cminor_dir+'RUN/TESTRUN/LLNL_nHeptane/LLNL_nHeptane_test.nc'
        atm = False
        plt_vars = np.array([["oh"], ["co2"], ["co"], ["h2o"], ["ho2"], ["nc7h16"]], dtype=object)
        spcnames = np.array(["OH", "CO$_2$", "CO", "H$_2$O", "HO$_2$", "C$_7$H$_{16}$"])

    elif mechanism=="LLNL_MD":
        name = "LLNL methyl-decanoate"
        sublabel = "(c)"
        fn=Cminor_dir+'RUN/TESTRUN/LLNL_MD/LLNL_MD_test.nc'
        atm = False
        plt_vars = np.array([["oh"], ["co2"], ["co"], ["h2o"], ["md"]], dtype=object)
        spcnames = np.array(["OH", "CO$_2$", "CO", "H$_2$O", "C$_{11}$H$_{22}$O$_2$"])

    ds   = nc.Dataset(fn)
    keys = ds.variables.keys()

    print('  Done reading dataset "'+fn+'".\n')

    nD = 0

    if "LWC_Level" in keys:
        LWC = True
    else:
        LWC = False

    for i in range(1, 500):
        if "pH_Value_"+str(i) in keys:
            nD = i
        else:
            break

    print("  Found nD = "+str(nD)+".\n")


    # possible conversion factor
    #c = 1 / 6.02295e17 * 0.0820574 * 298.15 * 0.001     # molec/cm3 -> partial pressure
    c = 1

    if atm:
        time_values = np.array(ds["time"]) # [hours]
    else:
        time_values = np.array(ds["time"])*1000 # [ms]

    n_t = time_values.size
    n_v = plt_vars.size


    values = np.zeros((n_v+1, n_t))
    ranges = np.zeros((n_v+1, 3))

    ##############################################################################
    # read values and plot variables
    ##############################################################################
    #plt.figure(figsize = (10, 2.5), dpi=150)
    if atm:
        fig, ax = plt.subplots(figsize = (10, 3.5), dpi=300)
    else:
        fig, ax = plt.subplots(figsize = (10, 3.5), dpi=300)

    for iVar in range(n_v):
        for ncdf_var in plt_vars[iVar]:
            values[iVar, :] += np.array(ds[ncdf_var])
        
        # scale [0, max] -> [0,1]
        if atm:
            ranges[iVar, :] = [np.min(values[iVar, :]), np.max(values[iVar, :]), np.max(values[iVar, :])-np.min(values[iVar, :])]
            if ncdf_var[0]=="a":
                if ncdf_var[-2:]=="m3":
                    plt.plot(time_values, values[iVar, :]/ranges[iVar, 1], label=spcnames[iVar]+" ["+"%.2e" % (6.022e+17*ranges[iVar, 1])+" molec cm$^{-3}$]")
                else:
                    plt.plot(time_values, values[iVar, :]/ranges[iVar, 1], label=spcnames[iVar]+" ["+"%.2e" % (ranges[iVar, 1])+" mol l$^{-1}$]")
            else:
                plt.plot(time_values, values[iVar, :]/ranges[iVar, 1], label=spcnames[iVar]+" ["+"%.2e" % (ranges[iVar, 1])+" molec cm$^{-3}$]")
        else:
            ranges[iVar, :] = [np.min(values[iVar, 1:]), np.max(values[iVar, 1:]), np.max(values[iVar, 1:])-np.min(values[iVar, 1:])]
            #plt.plot(time_values[1:], values[iVar, 1:]/ranges[iVar, 1], label=spcnames[iVar]+" ["+"%.2e" % (ranges[iVar, 1])+" mol mol$^{-1}$]")
            plt.plot(time_values[1:], values[iVar, 1:], label=spcnames[iVar])


    if atm:
        last_var  = "Zenith"
        last_name = "Solar altitude angle"
    else:
        last_var  = "Temperature"
        last_name = last_var

    ##############################################################################
    # read values and plot additional, scenario-specific variables and plot settings
    ##############################################################################
    iVar = n_v
    ncdf_var = last_var
    values[iVar, :] += np.array(ds[ncdf_var])
    if atm:
        values[iVar, :] = -(values[iVar, :]-np.pi/2) * (180/np.pi)
        values[iVar, :][values[iVar, :]<0] = 0.0

        ranges[iVar, :] = [np.min(values[iVar, :]), np.max(values[iVar, :]), np.max(values[iVar, :])-np.min(values[iVar, :])]

        plt.plot(time_values, values[iVar, :]/ranges[iVar, 1], label=last_name+" ["+"%.0f" % (ranges[iVar, 1])+"°]", alpha=0.5, color='gold', linestyle='--')
        if LWC:
            LWC_val = np.array(ds["LWC_Level"])
            LWC_min = np.min(LWC_val)
            LWC_max = np.max(LWC_val)
            plt.fill_between(time_values, LWC_val/LWC_max, label="LWC"+" ["+"%.1e" % LWC_max+" l m$^{-3}$]", alpha=0.3, color='deepskyblue')

        #ax.set_ylabel("Concentration [molec/cm$^3$]")
        ax.set_xlabel("Time [h]")
        ax.set_ylabel("Concentration")
        ax.legend(loc='center left', bbox_to_anchor=(1.01, 0.5))
        plt.subplots_adjust(right=0.65)
    else:
        plt.yscale('log')
        ranges[iVar, :] = [np.min(values[iVar, :]), np.max(values[iVar, :]), np.max(values[iVar, :])-np.min(values[iVar, :])]
        
        ax2 = plt.twinx()
        ax2.plot(time_values[1:], values[iVar, 1:], 'k--')
        ax.plot([0.0], [0.0], 'k--', label="Temperature")

        ax.set_xlabel("Time [ms]")
        ax.set_ylabel("Mole Fraction [mol/mol]")
        ax2.set_ylabel("Temperature [K]")

        ax.legend(loc='center left', bbox_to_anchor=(1.07, 0.5))
        plt.subplots_adjust(right=0.85)


    ##############################################################################
    # general plot settings
    ##############################################################################

    if sublabel!="":
        plt.text(-0.08*(plt.xlim()[1] - plt.xlim()[0]), plt.ylim()[1], sublabel, fontsize=15)

    plt.title(name)
    plt.xlim(time_values[0], time_values[-1])
    plt.tight_layout()
    ax.grid(which='major', axis='both')

    plt.savefig("PYTHONSCRIPTS/Figures/"+LABEL+"_trajectories.pdf")
    #plt.show()

print("  Finished gas distribution plots.\n")
