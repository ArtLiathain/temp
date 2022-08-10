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

# allow old-style pythonpath or new module imports
import sys
#from ffeatools import wrap
from ffeatools import ffea_script_jl as ffea_script_jl
from ffeatools import ffea_rod_jl as FFEA_rod
#from ffeatools import ffea_rod_jl.cc_extractor as ndc_extractor
#import ndc_extractor = FFEA_rod.cc_extractor


#try:
#    import wrap
#    import FFEA_script
#    import FFEA_rod
#    ndc_extractor = FFEA_rod.cc_extractor
#except ImportError:
#    from ffeatools import wrap
#    from ffeatools import FFEA_script
#    from ffeatools import FFEA_rod
#    import ffeatools.rod.cc_extractor as ndc_extractor


import numpy as np
import matplotlib.pyplot as plt

#import sys
#sys.path.insert(0, '/home/rob/pCloudSync/scratchpad/ndc_full_reader/') # not my fault


def plot_energy_histogram(energy,
                          energy_name,
                          bins=100,
                          bin_range=None,
                          scatter_plot_mode=False,
                          do_plot=True,
                          line=None,
                          col='c'):
    """
    Uses matplotlib to create energy histograms.
    """
    flat_energy = energy[1:].flatten()
    if bin_range is None:
        bin_range = (0,np.max(flat_energy)*0.5)
    if do_plot:
        values, bins, patches = plt.hist(flat_energy,
                                         bins,
                                         range=bin_range,
                                         density=True,
                                         color=col,
                                         edgecolor='k')
        #plt.title(energy_name)
        plt.ylabel("Probability")
        plt.xlabel("Energy (J)")
        if line:
            plt.axvline(line, color='k', linestyle='dashed', linewidth=1)
        plt.savefig(energy_name+".pdf")
        plt.close()
    else:
        values, bins = np.histogram(flat_energy, bins=bins, range=bin_range)
    if scatter_plot_mode:
        bin_width = (bins[1] - bins[0])
        bins = bins+(bin_width/2.)
        bins = bins[:-1]
    return bins, values

def get_x(stretch_analysis, twist_analysis, bend_analysis):
    """
    Gets the x values as understood when converting 
    to Python 3 in August 2022
    """
    delta_omega, L_i = FFEA_rod.ndc_extractor.get_delta_omega(bend_analysis, fast=False)
    delta_omega = delta_omega.reshape( [np.shape(delta_omega)[0]*np.shape(delta_omega)[1], 2] )
    delta_omega_1 = delta_omega[:,0]
    delta_omega_2 = delta_omega[:,1]

    p_i = stretch_analysis.rod.get_p_i(stretch_analysis.rod.current_r)
    p_i_equil = stretch_analysis.rod.get_p_i(stretch_analysis.rod.equil_r)
    delta_p = np.linalg.norm(p_i, axis=2) - np.linalg.norm(p_i_equil, axis=2)

    twist_analysis.get_twist_amount(set_twist_amount=True)
    delta_theta = twist_analysis.twist_amount
    return delta_p.flatten(), delta_omega_1, delta_omega_2, delta_theta.flatten()

def plot_x_histograms(analyses, bin_no=100, bin_range=None):
    """
    Function to use matplotlib to create histograms as understood
    when converting to Python 3 in August 2022
    """

    def hist_plot(x, plotname, analytical_formula, k, xlabel="Thing", ylabel="Probability (normalized)", num_bins=None, bin_range=None, col='c'):
        values, bins, patches = plt.hist(x, bins=num_bins, range=bin_range, density=True, color=col)
        P = analytical_formula(300, bins, k)
        P_norm = normalize_P(P, bins)
        #plt.title(energy_name)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.plot(bins, P_norm, color="k")
        plt.savefig(plotname+".pdf")
        plt.close()

    delta_x_dof = get_x(analyses[0], analyses[1], analyses[2])

    L_i_str = np.average(np.linalg.norm(analyses[0].rod.get_p_i(analyses[0].rod.equil_r), axis=2))
    L_i_bend = np.average(np.linalg.norm(analyses[2].rod.get_p_i(analyses[2].rod.equil_r), axis=2))
    L_i_twist = np.average(np.linalg.norm(analyses[1].rod.get_p_i(analyses[1].rod.equil_r), axis=2))

    hist_plot(delta_x_dof[0],
              "stretch_x",
              get_P,
              analyses[0].rod.material_params[0][0][0]/L_i_str,
              "$\Delta P$",
              col="b",
              num_bins=bin_no)
    hist_plot(delta_x_dof[1],
              "bend_omega_1_x",
              get_P,
              analyses[2].rod.B_matrix[0][0][0]/L_i_bend,
              "$\Delta \omega_1$",
              col="r",
              num_bins=bin_no)
    hist_plot(delta_x_dof[2],
              "bend_omega_2_x",
              get_P,
              analyses[2].rod.B_matrix[0][0][0]/L_i_bend,
              "$\Delta \omega_2$",
              col="r",
              num_bins=bin_no)
    hist_plot(delta_x_dof[3],
              "twist_x",
              get_P,
              analyses[1].rod.material_params[0][0][1]/L_i_twist,
              r'$\Delta \theta $',
              col="g",
              num_bins=bin_no)

def get_P(T, bins, k): # bins are P
    P = []
    for curr_bin in bins:
        P.append(         np.exp(-0.5 * k/(T*1.38064852e-23) * curr_bin**2)       )
    return np.array(P)

def get_E_1deg(T, bins, k): # bins are E
    E = []
    for curr_bin in bins:
        E.append( 2*np.exp( -curr_bin/(T*1.38064852e-23) )/np.sqrt(2*curr_bin*k) )
    return np.array(E)

def get_E_2deg(T, bins, k):
    E = []
    for curr_bin in bins:
        E.append( ((2*np.pi)/k)*np.exp( -curr_bin/(T*1.38064852e-23) ) )
    return np.array(E)

def normalize_P(P, bins):
    area = 0
    for bin_interval in range(len(bins)-1):
        area += ((P[bin_interval]+P[bin_interval+1])/2. *
                 abs(bins[bin_interval+1]-bins[bin_interval]))
    scale_factor = 1./area
    return np.array(P)*scale_factor


def main():

    #try:
    #    wrap.wrap_process("../../../../src/ffea", ["realistic.ffea"])
    #except OSError:
    #    raise

    script = FFEA_script.FFEA_script("2/realistic.ffea")

    temp = 300 #kelvin
    analytical_kbT = temp*1.38064852 * 10**-23

    bendy_rod = script.rod[0]
    stretchy_rod = script.rod[1]
    twisty_rod = script.rod[2]

    #Bendy rod
    bendy_analysis = FFEA_rod.anal_rod(bendy_rod)
    bendy_analysis.get_equipartition()
    bendy_rod_avg_energy = np.average(bendy_analysis.bending_energy)
    print("Bendy energy = "+str(bendy_rod_avg_energy))
    # 2 axes of bending = more D.O.F., see paper!
    print("Equipartition energy = "+str(1*analytical_kbT))
    bend_equal = FFEA_rod.rod_math.approximately_equal(1*analytical_kbT, bendy_rod_avg_energy, 0.08)
    print("Bendy status: "+str(bend_equal))
    print("---")

    #Stretchy rod
    stretchy_analysis = FFEA_rod.anal_rod(stretchy_rod)
    stretchy_analysis.get_equipartition()
    stretchy_rod_avg_energy = stretchy_analysis.stretch_energy
    stretchy_time_avg_energy = np.average(stretchy_rod_avg_energy)
    print("Stretchy energy = "+str(stretchy_time_avg_energy))
    print("Equipartition energy = "+str(0.5*analytical_kbT))
    stretch_equal = FFEA_rod.rod_math.approximately_equal(0.5*analytical_kbT,
                                                          stretchy_time_avg_energy,
                                                          0.08)
    print("Stretchy status: "+str(stretch_equal))
    print("---")

    #Twisty rod
    twisty_analysis = FFEA_rod.anal_rod(twisty_rod)
    twisty_analysis.get_equipartition()
    twisty_rod_avg_energy = np.average(twisty_analysis.twist_energy[:100])
    print("Tiwsty energy = "+str(twisty_rod_avg_energy))
    print("Equipartition energy = "+str(0.5*analytical_kbT))
    twist_equal = FFEA_rod.rod_math.approximately_equal(0.5*analytical_kbT,
                                                        twisty_rod_avg_energy,
                                                        0.08)
    print("Twisty status: "+str(twist_equal))
    print("---")

    #plot_histogram(bendy_analysis.bending_energy,
    #                "bending_energy",
    #                bins=50,
    #                line=analytical_kbT,
    #                col='b')
    #plot_histogram(twisty_analysis.twist_energy,
    #                "twsting_energy",
    #                bins=50,
    #                line=analytical_kbT/2.0,
    #                col='g')
    #plot_histogram(stretchy_analysis.stretch_energy,
    #                "stretch_energy",
    #                bins=50, line=analytical_kbT/2.0,
    #                col='r')

    plot_x_histograms([stretchy_analysis, twisty_analysis, bendy_analysis],
                      bin_no=100,
                      bin_range=None)

    if twist_equal and bend_equal and stretch_equal:
        sys.exit(0)
    else:
        sys.exit(1)

#if __name__ == "__main__":
#    main()
