
#
# Copyright (C) 2025 Levin Rug (E-Mail: l.rug@lmu.de)
# See ./SRC/Cminor.f90 for the copyright notice
# See ./LICENSE for license information
# SPDX-License-Identifier: GPL-3.0
#

import netCDF4 as nc
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import fileinput

class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush() # If you want the output to be visible immediately
    def flush(self) :
        for f in self.files:
            f.flush()

def ncdfcompare(fn1, fn2, mode, rel_tol = 1E-3, abs_tol=0.0, eps=0.0):

    # rel_tol is the relative tolerance every value has to fulfill
    # abs_tol is the absolute tolerance every value has to fulfill
    # eps is a small value, below which values are not compared

    # make mode specific default values
    if abs_tol == 0.0:
        if mode=="atm":
            abs_tol = 1000
        if mode=="comb":
            abs_tol = 0.0001

    # make mode specific default values
    if eps == 0.0:
        if mode=="atm":
            eps = 1000
        if mode=="comb":
            eps = 1E-8

    global_maxrel      = -1.
    global_maxabs      = -1.
    global_maxrelkey   = "None"
    global_maxabskey   = "None"
    global_maxrelvals  = np.array([-1., -1.])
    global_maxabsvals  = np.array([-1., -1.])

    # read data sets out of netcdf files
    ds1=nc.Dataset(fn1)
    ds2=nc.Dataset(fn2)

    # collect all keys (variable names)
    keys1 = ds1.variables.keys()
    keys2 = ds2.variables.keys()
    #if keys1 != keys2:
    #    print("")
    #    print(" ERROR :: Different keys in ncdfcompare :")
    #    print("  ", (set(keys1)-set(keys2)) | (set(keys2)-set(keys1)))
    #    sys.exit()

    #keys = set(ds1.variables.keys())

    if keys1 != keys2:
        print("")
        print(" WARNING :: Different keys in ncdfcompare :")
        print("  ", (set(keys1)-set(keys2)) | (set(keys2)-set(keys1)))

    keys = set(keys1) & set(keys2)

    # delete keys that may vary
    keys = list(keys - {"trajectory", "Step_Size"})

    if "LWC_Level" in keys:
        arr1 = np.array(ds1["LWC_Level"]) + 1E-100
        arr2 = np.array(ds2["LWC_Level"]) + 1E-100
        #errorLWC = sum(abs(arr1-arr2))
        errorLWC = np.max(abs(arr1-arr2))
        if errorLWC>1E-5:
            print("ERROR in NetCDFcompare :: LWC levels are not equal! (Error: "+"%.3e"%(errorLWC)+") Results cannot be expected to be the same. Aborting.")
            return (np.nan, np.nan, np.nan), (np.nan, np.nan, np.nan)


    # loop through all keys and check maximum absolute and relative differences
    for key in keys:

        arr1 = np.array(ds1[key]) + 1E-100
        arr2 = np.array(ds2[key]) + 1E-100

        absdiff = abs(arr1 - arr2)
        reldiff = np.zeros_like(absdiff)
        for i in range(reldiff.size):
            reldiff[i] = abs((arr1[i]-arr2[i])/min(arr1[i],arr2[i]))

        max_rel_error = max(reldiff)
        max_abs_error = max(absdiff)

        iMaxRel = np.argmax(reldiff)
        iMaxAbs = np.argmax(absdiff)

        val1rel = arr1[iMaxRel]
        val2rel = arr2[iMaxRel]
        val1abs = arr1[iMaxAbs]
        val2abs = arr2[iMaxAbs]

        # update global values
        if max(val1rel, val2rel)>eps and max_rel_error>global_maxrel:
            global_maxrel     = max_rel_error
            global_maxrelkey  = key
            global_maxrelvals = np.array([val1rel, val2rel])
        if max(val1abs, val2abs)>eps and max_abs_error>global_maxabs:
            global_maxabs     = max_abs_error
            global_maxabskey  = key
            global_maxabsvals = np.array([val1abs, val2abs])

        #if max(val1rel, val2rel)>eps and reldiff[iMaxRel]>rel_tol:
        #    print("  WARNING :: Large relative difference: ", max(abs(reldiff)), " at values ", val1rel, val2rel)

        #if max(val1abs, val2abs)>eps and absdiff[iMaxAbs]>abs_tol:
        #    print("  WARNING :: Large absolute difference: ", max(abs(absdiff)), " at values ", val1abs, val2abs)

    return (global_maxrel, global_maxrelkey, global_maxrelvals), (global_maxabs, global_maxabskey, global_maxabsvals)
