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
from math import sqrt, pi
import numpy as np
from ffeatools import ffea_script


# INPUT STUFF
FILEIN = "cyl_160_fine-EI.ffea"
E_t = 6e8 ## the Young's modulus we put in
F = 1e-12
L = 160e-9
r = 10e-9
I_t = pi * r**4 / 4.


ENDNODES = [60, 61, 63, 66, 67, 69, 72, 73, 75, 78,
            79, 81, 84, 85, 87, 90, 91, 93, 96, 97,
            99, 102, 103, 105, 108, 109, 111, 114, 116, 156,
            160, 164, 168, 172, 176, 180, 184, 188, 4310, 4311,
            4312, 4313, 4314, 4315, 4316, 4317, 4318, 4319, 4320, 4321,
            4322, 4323, 4324, 4325, 4326, 4327, 4328, 4329, 4330, 4331,
            4332, 4333, 4334, 4335, 4336, 4337, 4338, 4339, 4340, 4341,
            4342, 4343, 4344, 4345, 4346, 4347, 4348, 4349, 4350, 4351,
            4352, 4353, 4354, 4355, 4356, 4357, 4358, 4359, 4360, 4361,
            4362, 4363, 4364, 4365, 4366, 4367, 4368, 4369, 4370, 4371,
            4372, 4373, 4374, 4375, 4376, 4377, 4378, 4379, 4380, 4381,
            4382, 4383, 4384, 4385, 4386, 4387, 4388, 4389, 4390, 4391,
            4392, 4393, 4394, 4395, 4396, 4397, 4398, 4399, 4400, 4401]

script = ffea_script.ffea_script(FILEIN)


trj = script.load_trajectory()

## ## SHOW CONVERGENCE ## ##
CM0 = np.array([0.,0.,0.])
for f in range(trj.num_frames):
    CM = np.array([0.0,0.0,0.0])
    for n in ENDNODES:
        CM += trj.blob[0][0].frame[f].pos[n]
    CM /= len(ENDNODES)
    if f == 0:
        CM0 = CM
    print("{CM} {sqrt( (CM[0] - CM0[0])**2 + (CM[1] - CM0[1])**2 + (CM[2] - CM0[2])**2}")



## ## ## GET dL ## ##
CM = []
for f in [0, trj.num_frames-1]:
    cm = np.array([0.0,0.0,0.0])
    for n in ENDNODES:
        cm += trj.blob[0][0].frame[f].pos[n]
    cm /= len(ENDNODES)
    CM.append(cm)


dL = sqrt( (CM[1][0] -
            CM[0][0])**2 +
            (CM[1][1] - CM[0][1])**2 +
            (CM[1][2] - CM[0][2])**2 )

## ## ## Calculate the effective flexural rigidity:
EI_c = F * L**3 / (3 * dL)
EI_t = E_t * I_t

print("EI_c: {EI_c}")
print("EI_t: {EI_t}")
print("EI_c/EI_t: {EI_c/EI_t}")
print("abs(EI_c/EI_t -1): {abs(EI_c/EI_t -1)}")

ERROR = 0
if abs (EI_c/EI_t - 1) > 0.1:
    ERROR = 1

sys.exit(ERROR)
