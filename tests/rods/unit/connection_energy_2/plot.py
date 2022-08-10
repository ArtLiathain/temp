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
import subprocess
import matplotlib.pyplot as plt


def read(log_file):
    """
    Reads the log file and puts the information into a dictionary.
    """
    results_dict = {0:{}, 1: {}, 2: {}, 3:{} }
    #for key, value in results_dict.iteritems():
    for key, value in results_dict.items():
        results_dict[key] = { "+x": {}, "-x":{}, "+y":{}, "-y":{}, "+z":{}, "-z":{} }

    with open(log_file) as log:
        for line in log.readlines():
            if "ENERGYPLOT" in line:
                line_list = line.split(" ")
                displacement = float(line_list[2])
                node_index = int(line_list[4])
                results_dict[node_index]["+x"][displacement] = float(line_list[6])
                results_dict[node_index]["+y"][displacement] = float(line_list[7])
                results_dict[node_index]["+z"][displacement] = float(line_list[8])
                results_dict[node_index]["-x"][displacement] = float(line_list[9])
                results_dict[node_index]["-y"][displacement] = float(line_list[10])
                results_dict[node_index]["-z"][displacement] = float(line_list[11])
            if "EXPLODING" in line:
                line_list = line.split(" ")

    return results_dict

def plot_one(plusdict, minusdict, node_index, axis, filename):
    """
    Creates a single plot.
    """
    x = []
    y = []
    #minuskeys = minusdict.keys()
    minuskeys = list(minusdict.keys())
    minuskeys.reverse()
    minusvalues = list(minusdict.values())
    minusvalues.reverse()
    for key in minuskeys:
        x.append(key*-1)
    for key in plusdict.keys():
        x.append(key)
    for value in minusvalues:
        y.append(value)
    for value in plusdict.values():
        y.append(value)

    plt.plot(x, y)
    plt.xlabel("Displacement")
    plt.ylabel("Energy")
    plt.title("Node "+str(node_index)+" energies in the "+axis+" axis")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.clf()
    plt.close()
    return x, y

def plot_all(results_dict):
    """
    Call plot_one to create all the plots.
    """
    for index in [0,1,2,3]:
        for axis in ["x", "y", "z"]:
            plot_one(results_dict[index]["+"+axis],
                     results_dict[index]["-"+axis],
                     index,
                     axis,
                     "node_"+str(index)+"_"+axis+"_axis.pdf")
    return

def auto():
    """
    A helper function.
    """
    with open("log.txt", "wb") as fileout:
        subprocess.call(["../../../../src/ffea", "connection_energy_2.ffeatest"], stdout=fileout)
        plot_all(read("log.txt"))
        sys.exit(0)

if __name__ == "__main__":
    auto()
