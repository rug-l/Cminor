
#
# Copyright (C) 2025 Levin Rug (E-Mail: l.rug@lmu.de)
# See ./SRC/Cminor.f90 for the copyright notice
# See ./LICENSE for license information
# SPDX-License-Identifier: GPL-3.0
#

import numpy as np
import matplotlib.pyplot as plt
import sys

LABEL = "kreidenweis2003_parcel"
LABEL = "LLNL_MD"
LABEL = "LLNL_nheptane"
LABEL = "ERC_nheptane"
LABEL = "SmallStratoKPP"
LABEL = "MCM32+CAPRAM40"
LABEL = "RACM+C24"
matrixnames = np.array(["Miter", "LU_Miter"])

plt.figure(figsize=(8,4), dpi=150)
for iMatrix, matrixname in enumerate(matrixnames):

    filename = "MATRICES/"+matrixname+"_"+LABEL+".SparseMat"

    in_matrix_section = False
    dims = np.array([])
    nnz = 0

    with open(filename) as file:
        for line in file:
            if not in_matrix_section and dims.size==0 and "Dimension" in line:
                pos1 = line.index("Dimension")
                pos2 = line.index("x")
                dims = np.array([ int(line[pos1+10:pos2]) , int(line[pos2+1:]) ])
                M = np.zeros(dims)

            if not in_matrix_section and "Nonzeros" in line:
                pos1 = line.index("Nonzeros")
                nnz = int(line[pos1+9:])

            if (line.rstrip()[-6:]=="MATRIX"):
                in_matrix_section = True
                continue
        
            if in_matrix_section and line.rstrip()=="":
                in_matrix_section = False
        
            if in_matrix_section:
                content = line.strip()
                i = content[:content.index(" ")]
                content = content[content.index(" "):]
                content = content.strip()
                j = content[:content.index(" ")]
                content = content[content.index(" "):]
                content = content.strip()
                val = content

                i = int(i)-1
                j = int(j)-1
                val = float(val)+1
                M[i,j] = val

    plt.subplot(1, matrixnames.size, iMatrix+1)
    plt.spy(M, markersize=.6)
    plt.xticks([])
    plt.yticks([])
    plt.xlabel("Non-zeros: "+str(nnz))

#plt.subplots_adjust(top=1.4, bottom=-.4)
plt.tight_layout()
plt.savefig("PYTHONSCRIPTS/Figures/"+LABEL+"_"+matrixname+"_sparsity"+".png")
#plt.show()
