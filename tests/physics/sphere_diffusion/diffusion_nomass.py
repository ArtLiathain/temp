#!/usr/bin/env python
"""
This script tests the diffusion of a sphere with no mass.
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

# This set of import exceptions looks odd to me. I do not understand
# what the aim of them is but I have left them here as it does no
# harm. It may be a check that the right virtual python environment
# is being used.
print("sys.argv[1]:- ")
print(sys.argv[1])
FFEATOOLS_FOUND = False
try:
    import ffeatools # python package
    #FFEATOOLS_FOUND = True
except ImportError:
    print("diffussion_nomass.py: Failure to import ffeatools module")
    sys.exit(1) # failure to import ffeatools
try:
    #import ffea_script = ffeatools.ffea_script
    from ffeatools import ffea_script as ffea_script
except ImportError:
    print("diffussion_nomass.py: Failure to import ffea_script module")
    sys.exit(1) # failure to import ffea_script
try:
    #ffea_material = ffeatools.ffea_material
    from ffeatools import ffea_material as ffea_material
except ImportError:
    print("diffussion_nomass.py: Failure to import ffea_material module")
    sys.exit(1) # failure to import ffea_material
FFEATOOLS_FOUND = True
# from ffeatools import ffea_script as ffea_script

# Load trajectory
START = 5000
# START = 0
END = 10000

script = ffea_script.ffea_script(sys.argv[1])
try:
    script = ffea_script.ffea_script(sys.argv[1])
except Exception:
    # This is a very general exception. I think ffea_scripts
    # handles this exception and it probably will not come back
    # this far.
    sys.exit("diffussion_nomass.py: Failure to process ffea script!")


traj = script.load_trajectory(start=START, num_frames=END-START)
END = traj.num_frames + START
stokes = script.load_stokes(0)

# Analyse trajctory in sets of 1ps and test against theoretical diffusion
r2 = [None for i in range(START, END-1)]

for f in traj.blob[0][0].frame:
    f.calc_centroid()

for i in range((END - START)-1):
    r2[i] = (traj.blob[0][0].frame[i + 1].get_centroid() -
             traj.blob[0][0].frame[i].get_centroid())
    r2[i] = np.dot(r2[i], r2[i])

r2mean = np.mean(r2, axis=0)
r2stdev = np.std(r2, axis=0)
r2err = r2stdev / np.sqrt(len(r2))

drag = 6 * np.pi * script.params.stokes_visc * sum(stokes.radius) * script.blob[0].scale
r2theory = 6 * (script.params.kT / drag) * script.params.dt * script.params.check

#print "Calculated Diffusion: <r^2> = %6.2e +/- %6.2e" % (r2mean, r2err)
#print "Theoretical Diffusion: <r^2> = %6.2e" % (r2theory)
#print "Percent Error: d<r^2>/<r^2> = %6.2f%%" % (100 * np.fabs(r2mean - r2theory) / r2theory)
ERROR = 0
TOLERANCE = 0.05
if np.fabs(r2mean - r2theory) / r2theory <= TOLERANCE:
    print("Non-Inertial Diffusion: correct \n"
          "\t{np.fabs(r2mean - r2theory) / r2theory }  < { TOLERANCE}")
else:
    print("Non-Inertial Diffusion: failed \n"
    "\t{np.fabs(r2mean - r2theory) / r2theory}  > {TOLERANCE}")
    ERROR = 1
sys.exit(ERROR)
