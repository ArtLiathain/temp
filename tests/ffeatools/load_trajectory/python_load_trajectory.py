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


# test for importing ffea_trjectory
try:
    import ffeatools  # python package
except ImportError:
    try:
        print("python_load_trajectory: Failure to import ffeatools")
        from ffeatools import ffea_trajectory
    except ImportError:
        print("python_load_trajectory: Failure to import ffea_trajectory")
        sys.exit(1) # failure to import


# test for using ffea_trajectory
try:
    from ffeatools import ffea_trajectory
    test_load_trajectory = ffea_trajectory.ffea_trajectory("unit_test_traj.ftj")
    sys.exit(0)
except OSError:
    print("python_load_trajectory: ")
    print("Could not find trajectroy file in test directory.")
    sys.exit(1)
except SystemExit as error:
    #print(f"Unexpected {error=}, {type(error)=}")
    sys.exit(0)
except BaseException as error:
    print("python_load_trajectory: Unexpected error")
    print(f"Unexpected {error=}, {type(error)=}")
    sys.exit(1)
