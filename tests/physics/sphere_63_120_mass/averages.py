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
import numpy as np


def read_and_average(file_in, first, last, heading):
    """
    This function reads the file in and returns the average. As understood when converting to
    Python 3 in June 2022.
    """
    data = []
    averages = {}
    fields = heading.keys()
    with open(file_in, 'r', encoding="utf-8") as file:
        while file.readline() != "Measurements:\n":
            continue

        for item in file.readline().split():
            data.append([item.strip()])
        for line in file:
            count = 0
            if line.count("RESTART"):
                continue
            for item in line.split():
                data[count].append(float(item))
                count += 1

    if first == 0:
        first = 1
    for item in data:
        for field in fields:
            if field == item[0]:
                print("{item[0]}  "
                      "{abs(np.mean(item[first:last])/heading[item[0]] - 1)}  "
                      "{np.std(item[first:last])}  "
                      "{len(item[first:last])}")
                #print("{}  {}  {}  {}".format(item[0],
                #                              abs(np.mean(item[first:last])/heading[item[0]] - 1),
                #                              np.std(item[first:last]),
                #                              len(item[first:last])))

                #print("next: {}".format(abs(np.mean(item[first:last])/heading[item[0]] - 1)))
                #print("next: {}".format(np.std(item[first:last])))
                #print("length: {}".format(len(item[first:last])))
                averages[item[0]] = [abs(np.mean(item[first:last])/heading[item[0]] - 1),
                                     np.std(item[first:last])]
    return averages

NODES = 63
KbT = 4.11e-21
energy = {"KineticEnergy": (KbT*(3*NODES)/2), "StrainEnergy": (KbT*(3*NODES - 6)/2)}
print("K_th: {energy['KineticEnergy']}")
print("E_th: {energy['StrainEnergy']}")
Tol = {"KineticEnergy": 0.04, "StrainEnergy":0.03}
INITIAL = 40
END = -1

filein = ["sphere_63_120_mass_measurement.out"]


ERROR = 0
for f in filein:
    print(f)
    averages_list = read_and_average(f, INITIAL, END, energy)
    for s in Tol:
        print(averages_list[s])
        if averages_list[s][0] < Tol[s]:
            #print s, ": correct ", A[s][0], " < ", Tol[s]
            print(f"{s}: correct {averages_list[s][0]} < {Tol[s]}\n")
        else:
            #print s, ": failed ", A[s][0], " > ", Tol[s]
            print(f"{s}: failed {averages_list[s][0]} > {Tol[s]}")
            ERROR = 1

    print("\n\n")


sys.exit(ERROR)
