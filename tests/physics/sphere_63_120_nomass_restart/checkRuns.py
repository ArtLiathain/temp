#!/usr/bin/env python
"""
  This file is part of the FFEA simulation package

  Copyright (c) by the Theory and Development FFEA teams,
  as they appear in the README.md file.

  FFEA is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  FFEA is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.

  To help us fund FFEA development, we humbly ask that you cite
  the research papers on the package.
"""


import sys



# This test only works for single threaded runs. On
# multi-threaded runs, threads consume a different amount
# of random numbers every time, and so this direct comparison
# will fail.
#
# Still, it is important to use checkpoints, in order to have
# uncorrelated streams, as  otherwise artifacts could rise.
#
# sphere_63_120_nomass_6-10steps.cpt must be the same as
# sphere_63_120_nomass_10steps.cpt. Exactly.



CPT = ["sphere_63_120_nomass_6-10steps.fcp",
       "sphere_63_120_nomass_10steps.fcp"]
STA = []
for filename in CPT:
    with open(filename, 'r', encoding="utf-8") as file:
        STA.append(file.readlines())

for item in range(2):
    if STA[0][item] != STA[1][item]:
        print("checkpoint files: {CPT[0]} and {CPT[1]} "
              "should be exactly the same.")
        sys.exit(1)


# measurements can be SLIGHTLY different, because the trajectory
# positions that will be used for restarts were not saved with
# infinite accuracy.

MSM = ["sphere_63_120_nomass_I_measurement.out",
       "sphere_63_120_nomass_II_measurement.out"]
STA = []
for filename in MSM:
    tmpT = []
    with open(filename, 'r', encoding="utf-8") as file:
        while file.readline() != "Measurements:\n":
            continue
        file.readline() # skip the header
        for line in file:
            if line.count("RESTART"):
                continue
            tmpB = []
            for m in line.split():
                tmpB.append(float(m))
        tmpT.append(tmpB)
    STA.append(tmpT)


D_MAX = 0.0e0
D_MIN = 100.0e0
for i in range(len(STA[0])):
    for j in range(len(STA[0][i])):
        try:
            DATA = abs((STA[1][i][j] - STA[0][i][j])/STA[0][i][j])
        except ZeroDivisionError:
            DATA = 0
        if DATA < D_MIN:
            Dmin = DATA
        if DATA > D_MAX:
            Dmax = DATA

print("DMax: {D_MAX}")
print("Dmin: {D_MIN}")

if Dmax > 5e-4:
# the max error I got using Intel and GCC and different compiler
# flags was < 2e-4.
    sys.exit(1)
else:
    sys.exit(0)
