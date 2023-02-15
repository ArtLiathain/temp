#!/usr/bin/env python
"""
This script tests the diffusion of a sphere with mass.
As understood when converting to Python 3 in August 2022.

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
from ffeatools import ffea_script
from ffeatools import ffea_material

# Load trajectory
START = 5000
# START = 0
END = 10000
script = ffea_script.ffea_script(sys.argv[1])

trajectory = script.load_trajectory(start=START, num_frames = END-START)
END = trajectory.num_frames + START
material = script.load_material(0)
topology = script.load_topology(0)

# Analyse trajctory in sets of 1ps and test against theoretical diffusion
r2 = [None for i in range(START, END-1)]

for f in trajectory.blob[0][0].frame:
    f.calc_centroid()

for index in range((END - START) - 1):
    r2[index] = (trajectory.blob[0][0].frame[index + 1].get_centroid() -
             trajectory.blob[0][0].frame[index].get_centroid())
    r2[index] = np.dot(r2[index], r2[index])

r2mean = np.mean(r2, axis=0)
r2stdev = np.std(r2, axis=0)
r2err = r2stdev / np.sqrt(len(r2))

mass = trajectory.blob[0][0].frame[0].calc_mass(topology, material)
r2theory = 3 * (script.params.kT / mass) * ((script.params.dt * script.params.check)**2)
#print "Calculated Diffusion: <r^2> = %6.2e +/- %6.2e" % (r2mean, r2err)
#print "Theoretical Diffusion: <r^2> = %6.2e" % (r2theory)
#print "Percent Error: d<r^2>/<r^2> = %6.2f%%" % (100 * np.fabs(r2mean - r2theory) / r2theory)
ERROR = 0
TOLERANCE = 0.1
if np.fabs(r2mean - r2theory) / r2theory <= TOLERANCE:
    print("Inertial Diffusion: correct \n"
          "\t {np.fabs(r2mean - r2theory) / r2theory}  < {TOLERANCE}")
else:
    print("Inertial Diffusion: failed \n"
          "\t {np.fabs(r2mean - r2theory) / r2theory}  > {TOLERANCE}")
    ERROR = 1
sys.exit(ERROR)
