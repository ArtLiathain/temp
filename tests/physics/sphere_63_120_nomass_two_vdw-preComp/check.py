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


def read_and_store(file_in):
    """
    This function reads the file in and stores the data. As understood when converting to
    Python 3 in August 2022.
    """
    data = []
    changed_data = {}
    with open(file_in, 'r', encoding="utf-8") as file:
        while file.readline() != "Measurements:\n":
            continue

        blob = ""
        readblob = False
        for item in file.readline().split():
            if item == "|":
                readblob = True
                continue
            if readblob:
                blob = item
                readblob = False
                print(item)
                continue
            data.append([blob+item])
        for line in file:
            cnt = 0
            if line.count("RESTART"):
                continue
            for item in line.split():
                data[cnt].append(float(item))
                cnt += 1

    for item in range(len(data)):
        changed_data[data[item][0]] = item
        data[item].pop(0)

    return data, changed_data


## compare the PreCompEnergy of the first step to a value that has been checked:
ERROR = 0
KNOWN_VALUE = 2.79280200000000000e-18
#KNOWN_VALUE = 1.0075630000000001e-18
#KNOWN_VALUE = 6.87781100000000001e-17

## ## ## STRAIGHT ## ## ##
filein = ["sphere_63_120_two_measurement.out"]
H, dH = read_and_store(filein[0])

if abs (H[dH["PreCompEnergy"]][-1] - KNOWN_VALUE) > 1e-32 :
    print("PreComputed energy should be {KNOWN_VALUE}, "
          "but was found to be {(H[dH['PreCompEnergy']][-1]))}")
    ERROR = 1


## ## ## PBC ## ## ##
filein = ["sphere_63_120_two-PBC_measurement.out"]
H, dH = read_and_store(filein[0])
if abs (H[dH["PreCompEnergy"]][-1] - KNOWN_VALUE) > 1e-32 :
    print("PBC PreComputed energy should be {KNOWN_VALUE}, "
          "but was found to be {(H[dH['PreCompEnergy']][-1])}")
    ERROR = 1



sys.exit(ERROR)
