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


def read_and_average(file_in, first, last, eq_energy, fields):
    """
    This function reads the file in and returns the average. As understood when converting to
    Python 3 in August 2022.
    """
    data = []
    averages = {}
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

                print("{item[0]} "
                      "{abs(np.mean(item[first:last])/eq_energy - 1)} "
                      "{np.std(item[first:last])} "
                      "{len(item[first:last])}")
                averages[item[0]] = [abs(np.mean(item[first:last])/eq_energy - 1),
                             np.std(item[first:last])]
    return averages

NODES = 63
KbT = 4.11e-21
ENERGY = KbT*(3*NODES - 6)/2
print("Equipartition th: {ENERGY}")
Tol = {"StrainEnergy":0.03}
INITIAL = 40
END = -1

filein = ["sphere_63_120_nomass_measurement.out"]


ERROR = 0
for f in filein:
    print(f)
    averages_list = read_and_average(f, INITIAL, END, ENERGY, Tol.keys())
    for s in Tol.keys():
        print(averages_list[s])
        if averages_list[s][0] < Tol[s]:
            print("{s}: correct {averages_list[s][0]} < {Tol[s]}")
        else:
            print("{s}: failed {averages_list[s][0]} > {Tol[s]}")
            ERROR = 1

    print("\n\n")


sys.exit(ERROR)
