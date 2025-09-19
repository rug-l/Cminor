
#
# Copyright (C) 2025 Levin Rug (E-Mail: l.rug@lmu.de)
# See ./SRC/Cminor.f90 for the copyright notice
# See ./LICENSE for license information
# SPDX-License-Identifier: GPL-3.0
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import netCDF4 as nc
from jaruga_data import *

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
plt.rc('font',**{'family':'serif','serif':['Times']})
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 13}

plt.rc('font', **font)

initial_color = 'black'
final_color   = 'olive'
row2_color    = 'black'
jaruga_color  = "orange"

fontsize_t=17
fontsize_l=13


print("")

fn='RUN/TESTRUN/kreidenweis2003_parcel/kreidenweis2003_parcel_reference.nc'

ds   = nc.Dataset(fn)
keys = ds.variables.keys()

print('  Done reading dataset "'+fn+'".\n')

nD = 0

for i in range(1, 5000):
    if "wetRadius_"+str(i) in keys:
        nD = i
    else:
        break


if nD!=50:
    print("  Found "+str(nD)+" droplet classes. The script only works for the exact conditions of the Kreidenweis et al. 2003 scenario.")

time_values       = np.array(ds["time"])
LWC_values        = np.zeros(time_values.size)
SO2_values        = np.zeros(time_values.size)
diss_values       = np.zeros((nD, time_values.size))
NH3_henry_values  = np.zeros((nD, time_values.size))
meanpH_values     = np.zeros(time_values.size)
SVI_values        = np.zeros((nD, 2))
NH4p_values       = np.zeros((nD, 2))
aerosol_diameters = np.zeros((nD, 2)) 

tstart = np.argmin(abs(time_values-196/3600))


n_activated = 14

dropletclass_numbers = np.array([1428754.4401235012      ,\
   2029517.2410624390      ,\
   2823320.1669173865      ,\
   3846447.0611924268      ,\
   5132059.0305209681      ,\
   6705880.1971429028      ,\
   8581282.7109945714      ,\
   10754269.371638030      ,\
   13199026.441113003      ,\
   15864821.209781893      ,\
   18675006.700654089      ,\
   21528742.176978022      ,\
   24305740.578107223      ,\
   26873941.949977785      ,\
   29099545.743447125      ,\
   30858397.230044454      ,\
   32047403.164526761      ,\
   32594524.271948516      ,\
   32465999.249659121      ,\
   31669793.120090187      ,\
   30254780.208793700      ,\
   28305778.272943079      ,\
   25935135.334596634      ,\
   23272029.239278316      ,\
   20450894.238888204      ,\
   17600404.982386172      ,\
   14834240.765341640      ,\
   12244478.266792417      ,\
   9898002.0238953233      ,\
   7835867.3821505308      ,\
   6075176.4142080545      ,\
   4612782.5224030018      ,\
   3430041.1124258041      ,\
   2497859.0140571594      ,\
   1781429.9865595102      ,\
   1244233.3100526333      ,\
   851073.57543301582      ,\
   570117.79124033451      ,\
   374019.52382957935      ,\
   240301.30533993244      ,\
   151199.42876362801      ,\
   93170.048003435135      ,\
   56225.667601585388      ,\
   33229.594268918037      ,\
   19233.023209214211      ,\
   10901.899003267288      ,\
   6051.8592127561569      ,\
   3290.0880628824234      ,\
   1751.6940444707870      ,\
   913.35789859294891     
], dtype=int)


mol2part = 6.0221408e+17

LWC_values = np.array(ds["LWC_Level"]) * 1000.0 / np.array(ds["rho_parcel"]) # g/kg
SO2_values = np.array(ds["SO2"]) / mol2part

[rho_0, rho_cloudbase, rho_end] = np.array(ds["rho_parcel"])[[0, tstart, -1]]

sum_drop_vols = np.zeros(time_values.size)

for i in range(nD):

    # add aqueous SO2 concentration as well (like in kreidenweis/jaruga)
    SO2_values += np.array(ds["aSO2_"+str(i+1)+"_m3"]) + np.array(ds["HSO3m_"+str(i+1)+"_m3"]) + np.array(ds["SO3mm_"+str(i+1)+"_m3"])
    SVI_values[i, :]  = np.array(ds["HSO4m_"+str(i+1)+"_m3"][[0, -1]]) + np.array(ds["SO4mm_"+str(i+1)+"_m3"][[0, -1]])
    NH4p_values[i, :] = np.array(ds["NH4p_"+str(i+1)+"_m3"][[0, -1]])           # mol/m3

    SVI_values[i, :]  = SVI_values[i, :]  * [1.0, rho_0/rho_end] * 97           # g in parcel
    NH4p_values[i, :] = NH4p_values[i, :] * [1.0, rho_0/rho_end] * 18           # g in parcel

    meanpH_values += 10**(-np.array(ds["pH_Value_"+str(i+1)])) * dropletclass_numbers[i] * np.array(ds["rho_parcel"]) / rho_0 \
                     * (4/3)*np.pi*np.array(ds["wetRadius_"+str(i+1)])**3 / np.array(ds["LWC_Level"]) * 1000

# turn SO2 to ppb
SO2_values *= 1000000000 * 8.31446261815324 * np.array(ds["Temperature"]) / np.array(ds["pressure"]) # ppb

aerosol_diameters[:, 0] = 2 * 1e6 * (3*((SVI_values[:, 0] + NH4p_values[:, 0]) / dropletclass_numbers / 1800000.0) / 4 / np.pi)**(1/3) # µm
aerosol_diameters[:, 1] = 2 * 1e6 * (3*((SVI_values[:, 1] + NH4p_values[:, 1]) / dropletclass_numbers / 1800000.0) / 4 / np.pi)**(1/3) # µm


print("initial aerosol mass: ")
print("  SO4mm = ", 10**6 * sum(SVI_values[:,0]), "µg/m3")
print("  NH4p  = ", 10**6 * sum(NH4p_values[:,0]), "µg/m3")


f, axes = plt.subplot_mosaic('FAAAAG;CCDDEE', figsize=(10,8), dpi=500)
axes['F'].axis('off')
axes['G'].axis('off')

# jaruga plots
axes['A'].scatter(jaruga_final_x, jaruga_final_y, color=jaruga_color, marker='.', label="Jaruga et al. 2018", zorder=100000)
axes['C'].scatter(jaruga_LWC_x, jaruga_LWC_y, color=jaruga_color, marker='.', label="Jaruga et al. 2018", zorder=100000)
axes['D'].scatter(jaruga_SO2_x, jaruga_SO2_y, color=jaruga_color, marker='.', label="Jaruga et al. 2018", zorder=100000)
axes['E'].scatter(jaruga_pH_x , jaruga_pH_y, color=jaruga_color, marker='.', label="Jaruga et al. 2018", zorder=100000)

axes['A'].plot(aerosol_diameters[1:,0]           , SVI_values[1:,0]*10**6/abs(np.log10(aerosol_diameters[:-1,0])-np.log10(aerosol_diameters[1:,0]))               , color=initial_color , marker='s', markersize=3, label="initial")
axes['A'].plot(aerosol_diameters[1:n_activated,1], SVI_values[1:n_activated,1]*10**6/abs(np.log10(aerosol_diameters[:n_activated-1,1])-np.log10(aerosol_diameters[1:n_activated,1])), color=final_color , marker='s', markersize=3, label="final")
axes['A'].plot(aerosol_diameters[n_activated:,1] , SVI_values[n_activated:,1]*10**6/abs(np.log10(aerosol_diameters[n_activated-1:-1,1])-np.log10(aerosol_diameters[n_activated:,1])), color=final_color , marker='s', markersize=3)
axes['A'].text(0.017, 1.5, "N = "+ "%d" % (sum(dropletclass_numbers[n_activated:])*rho_cloudbase/rho_0/10**6)+" cm"+r"$^{-3}$", fontsize=fontsize_t)

axes['A'].set_xscale('log')
axes['A'].set_yscale('log')
axes['A'].set_xlabel("particle diameter (µm)", fontsize=fontsize_t)
axes['A'].set_ylabel("dS(VI)/d(logD)\n(µg/m"+r"$^3$"+" per log"+r"$_{10}$"+" size interval)", fontsize=fontsize_t)
axes['A'].grid()
axes['A'].legend(fontsize=fontsize_l)

axes['C'].plot(LWC_values[tstart:], time_values[tstart:]*3600-196, color=row2_color)
axes['C'].grid()
axes['D'].plot(SO2_values[tstart:], time_values[tstart:]*3600-196, color=row2_color)
axes['D'].grid()
axes['E'].plot(-np.log10(meanpH_values[tstart:]), (time_values[tstart:]*3600-196)/2000, color=row2_color)
axes['E'].grid()

axes['C'].set_xlabel("LWC (g/kg)", fontsize=fontsize_t)
axes['C'].set_ylabel("time above cloud base (s)", fontsize=fontsize_t)
axes['C'].yaxis.set_ticks(np.arange(0.0, 2600, 200))
axes['C'].set_ylim([0.0, 2400.0])
axes['C'].set_xlim([0.0, 2.5])
axes['C'].xaxis.set_ticks(np.arange(0.0, 3.0, 0.5))

axes['D'].set_xlabel(r"SO$_2$ (ppb)", fontsize=fontsize_t)
axes['D'].get_yaxis().set_ticklabels([])
axes['D'].set_ylim([0.0, 2400.0])
axes['D'].set_xlim([0.0, 0.2])
axes['D'].yaxis.set_ticks(np.arange(0.0, 2600, 200))
axes['D'].xaxis.set_ticks(np.arange(0.0, 0.205, 0.05))

axes['E'].set_xlabel("(volume-)mean pH", fontsize=fontsize_t)
axes['E'].yaxis.set_label_position("right")
axes['E'].yaxis.tick_right()
axes['E'].set_ylabel("height above cloud base (km)", fontsize=fontsize_t)
axes['E'].yaxis.set_ticks(np.arange(0.0, 1.3, 0.1))
axes['E'].set_ylim([0.0, 1.2])
axes['E'].set_xlim([3.8, 5.0])
axes['E'].xaxis.set_ticks(np.arange(3.6, 5.2, 0.2))
axes['E'].set_xticklabels(["", "3.8", "4.0", "4.2", "4.4", "4.6", "4.8", "5.0"])

axes['A'].set_xlim([0.01,1.0])
axes['A'].xaxis.set_ticks([0.01, 0.1, 1.0])
axes['A'].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes['A'].set_xticklabels(["0.01", "0.1", "1"])
axes['A'].set_ylim([0.01,60.0])
axes['A'].yaxis.set_ticks([0.01, 0.1, 1.0, 10])
axes['A'].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes['A'].set_yticklabels(["0.01", "0.1", "1", "10"])

for ax in axes.values():
    ax.tick_params(labelsize=fontsize_l)

plt.tight_layout()
plt.savefig("/Users/rug/Cminor/PYTHONSCRIPTS/Figures/kreidenweis2003_figure1.png")
#plt.show()

print("  Finished Kreidenweis2003 figure.\n")
