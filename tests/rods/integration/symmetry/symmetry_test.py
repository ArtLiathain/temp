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
from ffeatools import wrap
from ffeatools import ffea_script_jl
#from ffeatools import ffea_rod_jl as FFEA_rod
#import ffeatools.rod.anal_rod as anal_rod
from ffeatools.rod import anal_rod


#try:
#    import wrap
#    import FFEA_script
#    import FFEA_rod
#except ImportError:
#    from ffeatools import wrap
#    from ffeatools import FFEA_script
#    from ffeatools import FFEA_rod

def main():
    """
    Main function for the rods symmetry test.
    """

    try:
        #wrap.wrap_process("ffeadev", ["symmetry_test_stretch_only.ffea"])
        wrap.wrap_process("../../../../src/ffea", ["symmetry_test_bend_only.ffea"])
        wrap.wrap_process("../../../../src/ffea", ["symmetry_test_stretch_only.ffea"])
    except OSError:
        print("symmetry_test.py: something went wrong.")
        sys.exit(1)

    #bend_script = FFEA_script.FFEA_script("symmetry_test_bend_only.ffea")
    bend_script = ffea_script_jl.ffea_script_jl("symmetry_test_bend_only.ffea")
    bend_rod = bend_script.rod[0]
    bend_rod.set_avg_energies()
    #bend_analysis = FFEA_rod.anal_rod(bend_rod)
    bend_analysis = anal_rod.analysis(bend_rod)
    bend_test_result = bend_analysis.do_bend_symmetry_test()

    #stretch_script = FFEA_script.FFEA_script("symmetry_test_stretch_only.ffea")
    stretch_script = ffea_script_jl.ffea_script_jl("symmetry_test_stretch_only.ffea")
    stretch_rod = stretch_script.rod[0]
    stretch_rod.set_avg_energies()
    #stretch_analysis = FFEA_rod.anal_rod(stretch_rod)
    stretch_analysis = anal_rod.analysis(stretch_rod)
    stretch_test_result = stretch_analysis.do_stretch_symmetry_test()

    if bend_test_result is False or stretch_test_result is False:
        sys.exit(1)

    sys.exit(0)



if __name__ == "__main__":
    main()
