# -*- coding: utf-8 -*-
#
#  This file is part of the FFEA simulation package
#
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file.
#
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
#
#  To help us fund FFEA development, we humbly ask that you cite
#  the research papers on the package.
#

"""
        FFEA_rod.py
        Author: Rob Welch, University of Leeds
        Email: py12rw@leeds.ac.uk
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from mpl_toolkits.mplot3d import Axes3D  # do not remove
try:
    import rod.ndc_extractor as cc_extractor  # for old-style module imports
except:
    pass

global rod_creator_version
rod_creator_version = 1.3


"""
 _____     ______
|_   _|   / / __ \
  | |    / / |  | |
  | |   / /| |  | |
 _| |_ / / | |__| |
|_____/_/   \____/

"""


class FFEA_rod:
    """
    Class for an FFEA rod object containing the raw data for an FFEA rod object.
    The raw data is contained in attributes containing 3-D numpy arrays. One
    dimension for time, another for the index, and the third for the positions.

    Attributes:
        Header:
        filename: filename of the rod that's been loaded.
        end_of_header: line number where the header ends.
        rod_id: unique identifier for the rod. not used outside of the simulation
        itself.
        num_elements: number of nodes that make up the rod.
        length: length of the 1-D array that containing the rod's positional data.
        For example, a 10-element long 3-D rod has a num_elements of 30.
        num_rods: total number of rods in the simulation.

        Contents:
        current_r - current position of the nodes in the rod.
        current_m - current vectors for the material frames.
        equil_r - equilibrium configuration of the rod
        equil_m - equilibrium material frame of the rod
        perturbed x, y, z and twist energy, positive and negative - results
        from each stage of the numerical derivatives taken during the simulation.
        material_params - each node has the structure [stretch, twist, blank]
        B_matrix - the parameter for the bending energy. This one is a 4-element
        array specifying the contents of the 2x2 B matrix.

        Dragons:
        unperturbed_energy_type - an estimate of the unperturbed energy given
        by type (e.g. stretch, bend, twist)
        unperturbed_energy_dof - an estimate of the unperturbed energy, by
        degree of freedom (x, y, z, twist)
    """

    def __init__(self, filename=None, rod_no=0, num_rods=1, num_elements=0):
        """
        Initialize the rod object and load the contents of the trajectory.

        Params:
            filename - the path to the .rodtraj file to be loaded.
            rod_no - give the rod a unique ID, if you like

        Returns:
            nothing. But it populates every attribute in this object, save for
            dragons.
        """
        self.rod_no = rod_no

        if filename:
            self.filename = filename
            self.end_of_header = 0
            rod_file = open(filename, "r")
            for line in rod_file:
                self.end_of_header += 1
                if line.split(',')[0] == 'HEADER':
                    self.rod_id = int(line.split(',')[2])
                if line.split(',')[0] == 'num_elements':
                    self.num_elements = int(line.split(',')[1])
                if line.split(',')[0] == 'length':
                    self.length = int(line.split(',')[1])
                if line.split(',')[0] == 'num_rods':
                    self.num_rods = int(line.split(',')[1])
                if line == "---END HEADER---\n":
                    break
            rod_file.close()
            self.get_trajectory_length()
            self.load_trajectory()

        else:
            self.num_frames = 1
            self.num_elements = num_elements
            self.length = self.num_elements * 3
            self.num_rods = num_rods
            self.equil_r = np.zeros([self.num_frames, self.num_elements, 3])
            self.equil_m = np.zeros([self.num_frames, self.num_elements, 3])
            self.current_r = np.zeros([self.num_frames, self.num_elements, 3])
            self.current_m = np.zeros([self.num_frames, self.num_elements, 3])
            self.internal_perturbed_x_energy_positive = np.zeros(
                [self.num_frames, self.num_elements, 3])
            self.internal_perturbed_y_energy_positive = np.zeros(
                [self.num_frames, self.num_elements, 3])
            self.internal_perturbed_z_energy_positive = np.zeros(
                [self.num_frames, self.num_elements, 3])
            self.twisted_energy_positive = np.zeros(
                [self.num_frames, self.num_elements, 3])
            self.internal_perturbed_x_energy_negative = np.zeros(
                [self.num_frames, self.num_elements, 3])
            self.internal_perturbed_y_energy_negative = np.zeros(
                [self.num_frames, self.num_elements, 3])
            self.internal_perturbed_z_energy_negative = np.zeros(
                [self.num_frames, self.num_elements, 3])
            self.twisted_energy_negative = np.zeros(
                [self.num_frames, self.num_elements, 3])
            self.material_params = np.zeros(
                [self.num_frames, self.num_elements, 3])
            self.B_matrix = np.zeros([self.num_frames, self.num_elements, 4])
            self.steric_energy = np.zeros(
                [self.num_frames, self.num_elements, 3])
            self.steric_force = np.zeros(
                [self.num_frames, self.num_elements, 3])
            self.num_neighbours = np.zeros(
                [self.num_frames, self.num_elements])

        return

    def get_num_dimensions(self, row):
        """
        Get the number of dimensions for the row. For example, equil_r is stored
        in row 1, and it's positions, so it's 3D. This will return 4 for
        B_matrix because the B_matrix of each node is 4 elements long.
        """
        rod_file = open(self.filename, "r")
        for lines_skipped in range(self.end_of_header):
            rod_file.readline()
        for lines_skipped in range(row):
            rod_file.readline()
        line = rod_file.readline()
        rod_file.close()
        # What if this isn't an integer? Potential bug here
        return int(len(line.split(",")) / self.num_elements)

    def get_trajectory_length(self):
        """
        Get the number of frames in the trajectory.
        Input: nothing, but the file header info must already be loaded (see
        __init__).
        Output: sets the self.num_frames attribute.
        """
        self.num_frames = 0
        rod_file = open(self.filename, "r")
        for line in rod_file:
            if line.split(" ")[0] == "FRAME" and line.split(" ")[2] == "ROD":
                self.num_frames += 1
        rod_file.close()

    def load_trajectory(self):
        """
        Loads the trajectory. Fast, but also messy. The basic idea is that it allocates
        a load of very large numpy arrays, then iterates over every frame, calling
        rod_file.readline() for each line, loading that frame for that line into
        the array, and then does .readline() again. It does that 14 times for the
        14 rows in each frame. No, it's not pretty, you're right. It's not.
        At least it's not trying to be clever.

        Inputs: none, it already has the filename set in self.filename and all
        the meta info.
        Outputs: none, but it populates all the 'contents' arrays.
        """

        self.equil_r = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(1)])
        self.equil_m = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(2)])
        self.current_r = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(3)])
        self.current_m = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(4)])
        self.internal_perturbed_x_energy_positive = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(5)])
        self.internal_perturbed_y_energy_positive = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(6)])
        self.internal_perturbed_z_energy_positive = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(7)])
        self.twisted_energy_positive = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(8)])
        self.internal_perturbed_x_energy_negative = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(9)])
        self.internal_perturbed_y_energy_negative = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(10)])
        self.internal_perturbed_z_energy_negative = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(11)])
        self.twisted_energy_negative = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(12)])
        self.material_params = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(13)])
        self.B_matrix = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(14)])
        self.steric_energy = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(15)])
        self.steric_force = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(16)])
        self.num_neighbours = np.empty([self.num_frames, self.num_elements])

        # look, this is not pretty but it is really fast
        # Hard-coded some shapes at the end because I wasn't following the previous convention

        rod_file = open(self.filename, "r")
        frame_no = 0
        while True:
            line = rod_file.readline()
            if not line:
                break
            if line.split(" ")[0] == "FRAME":
                try:
                    row = 1
                    line = rod_file.readline()
                    self.equil_r[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.equil_r[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.equil_m[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.equil_m[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.current_r[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.current_r[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.current_m[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.current_m[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.internal_perturbed_x_energy_positive[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.internal_perturbed_x_energy_positive[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.internal_perturbed_y_energy_positive[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.internal_perturbed_y_energy_positive[frame_no].shape)
                    line = rod_file.readline()
                    row += 1
                    self.internal_perturbed_z_energy_positive[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.internal_perturbed_z_energy_positive[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.twisted_energy_positive[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.twisted_energy_positive[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.internal_perturbed_x_energy_negative[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.internal_perturbed_x_energy_negative[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.internal_perturbed_y_energy_negative[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.internal_perturbed_y_energy_negative[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.internal_perturbed_z_energy_negative[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.internal_perturbed_z_energy_negative[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.twisted_energy_negative[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.twisted_energy_negative[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.material_params[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.material_params[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.B_matrix[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.B_matrix[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.steric_energy[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.steric_energy[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.steric_force[frame_no] = np.fromstring(
                        line, sep=",").reshape(self.steric_force[frame_no].shape)
                    row += 1
                    line = rod_file.readline()
                    self.num_neighbours[frame_no] = np.fromstring(
                        line, sep=",", dtype=int)
                    frame_no += 1
                except ValueError as e:
                    raise ValueError(str(e) + "\nError loading frame " +
                                     str(frame_no) + ", header row " + str(row) +
                                     "\nProblem line: " + str(line))
        # self.set_avg_energies()
        return

    def write_rod(self, filename):

        try:
            self.p_i
        except AttributeError:
            self.p_i = self.get_p_i(self.current_r)

        # First, write header info
        rod_file = open(filename, "w")
        rod_file.write("format,ffea_rod\n")
        rod_file.write("version," + str(rod_creator_version) + "\n")
        rod_file.write("HEADER,ROD," + str(self.rod_no) + "\n")
        rod_file.write("num_elements," + str(self.num_elements) + "\n")
        rod_file.write("length," + str(self.length) + "\n")
        rod_file.write("num_rods," + str(self.num_rods) + "\n")
        rod_file.write("row1,equil_r\n")
        rod_file.write("row2,equil_m\n")
        rod_file.write("row3,current_r\n")
        rod_file.write("row4,current_m\n")
        rod_file.write("row5,internal_perturbed_x_energy_positive\n")
        rod_file.write("row6,internal_perturbed_y_energy_positive\n")
        rod_file.write("row7,internal_perturbed_z_energy_positive\n")
        rod_file.write("row8,twisted_energy_positive\n")
        rod_file.write("row9,internal_perturbed_x_energy_negative\n")
        rod_file.write("row10,internal_perturbed_y_energy_negative\n")
        rod_file.write("row11,internal_perturbed_z_energy_negative\n")
        rod_file.write("row12,twisted_energy_negative\n")
        rod_file.write("row13,material_params\n")
        rod_file.write("row14,B_matrix\n")
        rod_file.write("row15,steric_energy\n")
        rod_file.write("row16,steric_force\n")
        rod_file.write("row17,num_neighbours\n")

        # Connections (note: this is temporary, it might end up in the .ffea file)
        rod_file.write("CONNECTIONS,ROD,0\n")
        rod_file.write("[rodelement], [blobno], [blobelement]\n")

        rod_file.write("---END HEADER---\n")

        # Write trajectory

        try:
            StringIO
            shutil
        except NameError:
            import cStringIO as StringIO
#            import shutil

        def write_array(array, file_obj):
            sio = StringIO.StringIO()
            np.savetxt(sio, array, newline=",")
            str_to_write = sio.getvalue()[:-1] + "\n"
            file_obj.write(str_to_write)

        for frame in range(self.num_frames):

            rod_file.write("FRAME " + str(frame) +
                           " ROD " + str(self.rod_no) + "\n")
            write_array(self.equil_r[frame].flatten(), rod_file)
            write_array(self.equil_m[frame].flatten(), rod_file)
            write_array(self.current_r[frame].flatten(), rod_file)
            write_array(self.current_m[frame].flatten(), rod_file)
            write_array(
                self.internal_perturbed_x_energy_positive[frame].flatten(), rod_file)
            write_array(
                self.internal_perturbed_y_energy_positive[frame].flatten(), rod_file)
            write_array(
                self.internal_perturbed_z_energy_positive[frame].flatten(), rod_file)
            write_array(
                self.twisted_energy_positive[frame].flatten(), rod_file)
            write_array(
                self.internal_perturbed_x_energy_negative[frame].flatten(), rod_file)
            write_array(
                self.internal_perturbed_y_energy_negative[frame].flatten(), rod_file)
            write_array(
                self.internal_perturbed_z_energy_negative[frame].flatten(), rod_file)
            write_array(
                self.twisted_energy_negative[frame].flatten(), rod_file)
            write_array(self.material_params[frame].flatten(), rod_file)
            write_array(self.B_matrix[frame].flatten(), rod_file)
            write_array(self.steric_energy[frame].flatten(), rod_file)
            write_array(self.steric_force[frame].flatten(), rod_file)
            write_array(self.num_neighbours[frame].flatten(), rod_file)

        rod_file.close()

    def get_p_i(self, x):
        """
        The rod trajectory files are only given in terms of the node positions.
        If we want p_i, we have to reconstruct it.
        Params: x, the array containing node positions (e.g. current_r, equil_r)
        for each frame.
        Returns: an array containing the p_i vectors for each frame.
        """
        # x = either self.current_r, or self.equil_r!

        # Find p_i
        x_shifted = np.copy(x)

        for frame in range(len(x_shifted)):
            x_shifted[frame] = np.roll(x_shifted[frame], 1, axis=0)
        p_i = x - x_shifted  # get all the p_i values in one big giant array operation

        # Remove first p_i
        p_i_new = np.empty(
            [self.num_frames, self.num_elements - 1, len(self.current_r[0][0])])
        for frame in range(len(p_i_new)):
            p_i_new[frame] = np.delete(p_i[frame], 0, axis=0)
        return p_i_new

    def set_avg_energies(self):
        """
        The rod has 4 degrees of freedon (x, y, z and twist) and 3 energy types
        (stretch, bend and twist). This function creates 2 arrays, filtering
        the energy by either DOF or energy type. It uses the average of the
        positive and negative perturbations done by the numerical differentation
        of the energy function. So it's not any more wrong than the differentiation
        itself. Right?
        This function isn't actually used in the rod analysis anywhere, it's mostly
        for testing.
        Params: none, but make sure to initialize the rod properly beforehand.
        Returns: none, but populates the attributes unperturbed_energy_type and
        unperturbed_energy_dof.
        """
        # this is dragons, do not use
        self.unperturbed_energy_type = np.empty(
            [self.num_frames, self.num_elements, self.get_num_dimensions(5)])
        self.unperturbed_energy_dof = np.empty(
            [self.num_frames, self.num_elements, 4])
        for frame in range(len(self.internal_perturbed_x_energy_positive)):
            for node in range(len(self.internal_perturbed_x_energy_positive[frame])):
                self.unperturbed_energy_type[frame][node] = (self.internal_perturbed_x_energy_positive[frame][node]
                                                             + self.internal_perturbed_x_energy_negative[frame][node]
                                                             + self.internal_perturbed_y_energy_positive[frame][node]
                                                             + self.internal_perturbed_y_energy_negative[frame][node]
                                                             + self.internal_perturbed_z_energy_positive[frame][node]
                                                             + self.internal_perturbed_z_energy_negative[frame][node]
                                                             + self.twisted_energy_positive[frame][node]
                                                             + self.twisted_energy_negative[frame][node]) / 8

                self.unperturbed_energy_dof[frame][node][0] = (np.sum(self.internal_perturbed_x_energy_positive[frame][node])
                                                               + np.sum(self.internal_perturbed_x_energy_negative[frame][node])) / 2
                self.unperturbed_energy_dof[frame][node][1] = (np.sum(self.internal_perturbed_y_energy_positive[frame][node])
                                                               + np.sum(self.internal_perturbed_y_energy_negative[frame][node])) / 2
                self.unperturbed_energy_dof[frame][node][2] = (np.sum(self.internal_perturbed_z_energy_positive[frame][node])
                                                               + np.sum(self.internal_perturbed_z_energy_negative[frame][node])) / 2
                self.unperturbed_energy_dof[frame][node][3] = (np.sum(self.twisted_energy_positive[frame][node])
                                                               + np.sum(self.twisted_energy_negative[frame][node])) / 2
        self.unperturbed_energy = np.sum(self.unperturbed_energy_type, axis=2)

    def scale(self, scale_factor):
        self.current_m *= scale_factor
        self.current_r *= scale_factor
        self.equil_r *= scale_factor
        self.equil_m *= scale_factor

    def translate(self, shift, translate_curr=True, translate_equil=True):
        #        self.current_m+=shift
        if translate_curr:
            self.current_r += shift
        if translate_equil:
            self.equil_r += shift
#        self.equil_m+=shift

    def calc_centroid(self, equil=False):
        if equil:
            self.centroid_traj = np.average(self.equil_r, axis=1)
        else:
            self.centroid_traj = np.average(self.current_r, axis=1)
        self.initial_centroid = self.centroid_traj[0]
        return self.centroid_traj

    def rotate(self, xyz):

        equil_centroid = self.calc_centroid(equil=True)[0]
        self.translate(-1 * equil_centroid, translate_curr=False)

        current_centroid = self.calc_centroid()
        self.translate(-1 * current_centroid, translate_equil=False)

        Rx = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(xyz[0]), np.sin(xyz[0])], [
                      0.0, np.sin(xyz[0]), np.cos(xyz[0])]])
        Ry = np.array([[np.cos(xyz[1.0]), 0, np.sin(xyz[1])], [
                      0.0, 1.0, 0.0], [-1 * np.sin(xyz[1]), 0.0, np.cos(xyz[1])]])
        Rz = np.array([[np.cos(xyz[2]), -1.0 * np.sin(xyz[2]), 0.0],
                       [np.sin(xyz[2]), np.cos(xyz[2]), 0], [0.0, 0.0, 1.0]])

        RyRx = np.matmul(Ry, Rx)
        RzRyRx = np.matmul(Rz, RyRx)

        for i in range(len(self.current_r[0])):
            self.equil_r[0][i] = np.matmul(RzRyRx, self.equil_r[0][i])
            self.current_r[0][i] = np.matmul(RzRyRx, self.current_r[0][i])
            self.current_m[0][i] = np.matmul(RzRyRx, self.current_m[0][i])
            self.equil_m[0][i] = np.matmul(RzRyRx, self.equil_m[0][i])

        self.translate(-1 * equil_centroid, translate_curr=False)
        self.translate(-1 * current_centroid, translate_equil=False)


"""
                       _           _
     /\               | |         (_)
    /  \   _ __   __ _| |_   _ ___ _ ___
   / /\ \ | '_ \ / _` | | | | / __| / __|
  / ____ \| | | | (_| | | |_| \__ \ \__ \
 /_/    \_\_| |_|\__,_|_|\__, |___/_|___/
                          __/ |
                         |___/
"""


class anal_rod:
    """
    The rod analysis class. This is where the actual tests happen. It's initialised
    using a rod object. It doesn't modify that object, but it does a load of
    measurements on it. It also plots the results of those measurements.

    Much like the rod itself, many quantities here are given as numpy arrays,
    of the form [frame][element][x,y,z] or sometimes [frame][element][quantity].

    For reference, the tests that it does are:
        - Equipartition of energy
        - Euler beam dynamics
        - Persistence length
        - Symmetry of energies and forces

    The attributes are:
        rod: rod object.
        rod_math: instance of the py_rod_math object, containing rod-specific
        math functions.
        analytical_P: value of the persistence length converted analytically.
        analytical_deflection: deflection calculated analytically (for the bending
        beam test).
        analytical_kbT: analytical thermal energy, for the equipartition test
        p_i, equil_p_i - array containing values of p_i for each frame and element
        deflections: value of deflection for each frame (bending beam test)
        EI_from_deflection: value of EI calculated from euler beam theory for
        each frame.
        P: persistence length (again, given for each frame)
        bending_energy: bending energy between each adjacent pair of elements
        for each frame. Note that this is not actually calculated during the
        simulation (only the derivatives are).
        twist_amount: the value of delta_theta - equil_delta_theta in the
        twist energy formula.
        twist_energy: the value of twist energy between pairs of elements
        in the rod.
        p_i_extension: similar to the above, but this time for the quantity of
        r-r_equil for the stretching energy. There's also p_i_extension_x, y and
        z.
        average_twist_energy: twist energy averaged over the entire rod, given
        by frame.
        EI: the analytical value of EI, also known as the bending modulus B.
        starting_length: the length of the rod (along the contours) on the
        first frame.
        average_extension: the average p_i_extension for each frame. Also given
        in terms of x, y and z.
        delta_x, delta_y, delta_z, delta_twist - how much the rod has moved
        between each frame!
    """

    def __init__(self, rod):
        """
        Initialize the rod analysis. This doens't do much, other than setting
        a few internal variables and firing up an instane of py_rod_math.
        Params: rod, an instance of FFEA_rod.
        Returns: nothing.
        """
        self.rod = rod
        self.py_rod_math = py_rod_math()
        self.analytical_P = None
        self.analytical_deflection = None
        self.analytical_kbT = None

        self.get_p_i = self.rod.get_p_i

    """
   _____                         _
  / ____|                       | |
 | |     ___  _ __ _ __ ___  ___| |_ _ __   ___  ___ ___
 | |    / _ \| '__| '__/ _ \/ __| __| '_ \ / _ \/ __/ __|
 | |___| (_) | |  | | |  __/ (__| |_| | | |  __/\__ \__ \
  \_____\___/|_|  |_|  \___|\___|\__|_| |_|\___||___/___/

    """

    def get_deflection(self):
        """
        Get the vertical deflection away from the equilibrium at the end of the rod.
        Used in the euler beam test.
        Params: none.
        Returns: none, but it sets self.deflections, an array.
        """
        frames = len(self.rod.current_r)
        nodes = len(self.rod.current_r[0]) - 1
        self.deflections = np.zeros([frames])
        for frame in range(frames):
            self.deflections[frame] = self.rod.current_r[frame][nodes][1] - \
                self.rod.equil_r[frame][nodes][1]

    def get_EI_from_deflection(self, force_applied):
        """
        This gets the approximate value of EI (B) for rods with isotropic
        bending response. W=(PL^3)/(3EI) where W is the deflection, P is the
        force applied to the end of the beam, L is the beam length.
        Params: force_applied, the force applied to the end of the beam in N (a
        float).
        Returns: none, but sets self.EI_from_deflection, an array.
        """
        try:
            self.deflections
        except AttributeError:
            self.get_deflection()

        frames = len(self.rod.current_r)
        nodes = len(self.rod.current_r[0]) - 1
        beam_length = rod_math.get_length(
            self.rod.equil_r[0][1] - self.rod.equil_r[0][nodes])
        self.EI_from_deflection = np.zeros([frames])
        for frame in range(frames):
            self.EI_from_deflection[frame] = (
                force_applied * (beam_length**3)) / (3 * rod_math.get_length(self.deflections[frame]))

    def get_persistence_length(self):
        """
        Get the value of the persistence length P from the current state of the
        rod.
        P = -L/log(<cos \theta>) where P is the persistence length, L is the
        length of the rod, and cos \theta is the cosine of the angle between
        the element at the start of the rod and the element at the end.
        Params: none.
        Returns: none, but populates self.P, an array containing a value for
        the persistence length for each frame.
        """
        try:
            self.p_i
        except AttributeError:
            self.p_i = self.get_p_i(self.rod.current_r)

        self.P = np.zeros(len(self.rod.current_r))
        cos_theta = np.zeros(len(self.rod.current_r))

        for frame in range(len(self.rod.current_r)):
            a = rod_math.normalize(self.p_i[frame][0])
            b = rod_math.normalize(self.p_i[frame][-1])
            cos_theta[frame] = rod_math.get_cos_theta(a, b)

        avg_cos_theta = np.average(cos_theta)
#        avg_cos_theta = np.sqrt(np.mean(cos_theta**2)) #rms?

        for frame in range(len(self.rod.current_r)):
            # -1 for there being 1 less element than nodes, -1 for the length being measured from the center of both elements.
            L = self.get_absolute_length(
                0, len(self.rod.current_r[0]) - 2, frame)
            self.P[frame] = -L / np.log(avg_cos_theta)

        if self.P[frame] > 0.0001:
            print("L = " + str(L))
            print("a = " + str(a))
            print("b = " + str(b))
            print("cos_theta = " + str(cos_theta))
            print("P = " + str(self.P[frame]))

    def get_entire_rod_stretch_energy(self):

        self.p_i = self.get_p_i(self.rod.current_r)
        self.equil_p_i = self.get_p_i(self.rod.equil_r)

        self.stretching_energy = np.zeros(
            [len(self.p_i), len(self.p_i[0]) - 1])

        for frame in range(len(self.p_i)):
            for element in range(len(self.stretching_energy[0])):
                self.stretching_energy[frame][element] = self.py_rod_math.get_stretch_energy(
                    self.get_constant_parameter(0), self.p_i[frame][element], self.equil_p_i[frame][element])

    def get_bending_response(self):
        """
        Get the bending energy between adjacent rod elements for every pair of
        rod elements in every frame. Mathematically, this is identical to what
        goes on inside the rod simulation proper (although it happens 8x more there).
        Params: none
        Returns: none, but it populates self.bending_energy, an array.
        """
        try:
            self.p_i
        except AttributeError:
            self.p_i = self.get_p_i(self.rod.current_r)
        try:
            self.equil_p_i
        except AttributeError:
            self.equil_p_i = self.get_p_i(self.rod.equil_r)

        try:
            self.get_constant_EI()
        except ValueError:
            import warnings
            warnings.warn(
                "EI is not constant. If this is for an equipartition test, it won't work.")

        self.bending_energy = np.zeros([len(self.p_i), len(self.p_i[0]) - 1])

        for frame in range(len(self.p_i)):
            for element in range(len(self.bending_energy[0])):

                # NOTE TO SELF: ADD NORMALISATION!

                B_j = np.matrix([[self.rod.B_matrix[frame][element][0], self.rod.B_matrix[0][0][1]], [
                                self.rod.B_matrix[0][0][2], self.rod.B_matrix[0][0][3]]])
                B_i = np.matrix([[self.rod.B_matrix[frame][element + 1][0], self.rod.B_matrix[0]
                                  [0][1]], [self.rod.B_matrix[0][0][2], self.rod.B_matrix[0][0][3]]])

                pim1 = self.p_i[frame][element]
                pi = self.p_i[frame][element + 1]
                mim1 = self.py_rod_math.normalize(
                    self.rod.current_m[frame][element])
                nim1 = np.cross(mim1, pim1 / np.linalg.norm(pim1))
                mi = self.py_rod_math.normalize(
                    self.rod.current_m[frame][element + 1])
                ni = np.cross(mi, pi / np.linalg.norm(pi))
                omega_i, equil_kbi = self.py_rod_math.omega(
                    pi, pim1, nim1, mim1)
                omega_j, equil_kbj = self.py_rod_math.omega(pi, pim1, ni, mi)

                equil_pim1 = self.equil_p_i[frame][element]
                equil_pi = self.equil_p_i[frame][element + 1]
                equil_mim1 = self.py_rod_math.normalize(
                    self.rod.equil_m[frame][element])
                equil_nim1 = np.cross(
                    equil_mim1, equil_pim1 / np.linalg.norm(equil_pim1))
                equil_mi = self.py_rod_math.normalize(
                    self.rod.equil_m[frame][element + 1])
                equil_ni = np.cross(equil_mi, equil_pi /
                                    np.linalg.norm(equil_pi))

                equil_omega_i, equil_kbi = self.py_rod_math.omega(
                    equil_pi, equil_pim1, equil_nim1, equil_mim1)
                equil_omega_j, equil_kbj = self.py_rod_math.omega(
                    equil_pi, equil_pim1, equil_ni, equil_mi)

                delta_omega_i = omega_i - equil_omega_i
                delta_omega_j = omega_j - equil_omega_j

                inner_i = np.dot(np.transpose(delta_omega_i),
                                 np.dot(B_j, delta_omega_i))
                inner_j = np.dot(np.transpose(delta_omega_j),
                                 np.dot(B_i, delta_omega_j))

                energy = 0.5 * (inner_i + inner_j) * (1 / (2 * (self.py_rod_math.get_length(
                    equil_pi) + self.py_rod_math.get_length(equil_pim1))))

                # shouldn't this be regular pi and pim1?
                l_i = self.py_rod_math.get_length(
                    equil_pi) + self.py_rod_math.get_length(equil_pim1)

                self.bending_energy[frame][element] = energy

    def get_bending_response_mutual(self, rotate=False):
        """
        Get the bending energy between adjacent rod elements for every pair of
        rod elements in every frame, using the mutual parallel transport method.
        Params: none
        Returns: none, but it populates self.bending_energy, an array.
        """

        def get_weights(a, b):
            a_length = rod_math.get_length(a)
            b_length = rod_math.get_length(b)
            weight1 = a_length / (a_length + b_length)
            return weight1

        def get_mutual_frame_inverse(pim1, pi, weights=None):
            if not weights:
                weights = get_weights(pim1, pi)
            mutual_element = (1 / weights) * rod_math.normalize(pim1) + \
                (1 / (1 - weights)) * rod_math.normalize(pi)
            return rod_math.normalize(mutual_element)

        def get_mutual_axes_inverse(mim1, mi, pim1=None, pi=None, weight=None):
            if not weight:
                weight = get_weights(pim1, pi)
            m_mutual = (mim1 * (1.0 / weight) + mi * (1.0 / (1 - weight))) / \
                (rod_math.get_length(mi) + rod_math.get_length(mim1))
            return rod_math.normalize(m_mutual)

        try:
            self.p_i
        except AttributeError:
            self.p_i = self.get_p_i(self.rod.current_r)
        try:
            self.equil_p_ir
        except AttributeError:
            self.equil_p_i = self.get_p_i(self.rod.equil_r)

        try:
            self.get_constant_EI()
        except ValueError:
            import warnings
            warnings.warn(
                "EI is not constant. If this is for an equipartition test, it won't work.")

        self.bending_energy = np.zeros([len(self.p_i), len(self.p_i[0]) - 1])
        #L_i = np.zeros( [len(self.p_i), len(self.p_i[0])-1] )
        #delta_omega_all = np.zeros( [len(self.p_i), len(self.p_i[0])-1, 2] )

        for frame in range(len(self.p_i)):
            for element in range(len(self.bending_energy[0])):

                # NOTE TO SELF: ADD NORMALISATION!

                B = np.matrix([[self.rod.B_matrix[frame][element + 1][0], self.rod.B_matrix[0]
                                [0][1]], [self.rod.B_matrix[0][0][2], self.rod.B_matrix[0][0][3]]])

                #mutual_l = get_mutual_frame_inverse(self.p_i[frame][element], self.p_i[frame][element+1])
                #equil_mutual_l = get_mutual_frame_inverse(self.equil_p_i[frame][element], self.equil_p_i[frame][element+1])

                pim1 = self.p_i[frame][element]
                pi = self.p_i[frame][element + 1]
                mim1 = self.py_rod_math.normalize(
                    self.rod.current_m[frame][element])
                nim1 = np.cross(mim1, pim1 / np.linalg.norm(pim1))
                mi = self.py_rod_math.normalize(
                    self.rod.current_m[frame][element + 1])
                ni = np.cross(mi, pi / np.linalg.norm(pi))

                # get weighting between elements adjacent to bending node!
                weight = get_weights(pim1, pi)

                # get weighted element at that node
                mutual_l = get_mutual_frame_inverse(pim1, pi, weight)

                # parallel transport our two existing material axes onto our new element
                mutual_mi = np.array(rod_math.parallel_transport(
                    mi,  rod_math.normalize(pi), mutual_l))[0]
                mutual_mim1 = np.array(rod_math.parallel_transport(
                    mim1,  rod_math.normalize(pim1), mutual_l))[0]

                # get the weighted average to create the mutual material axes
                mutual_m_rotated = get_mutual_axes_inverse(
                    mutual_mim1, mutual_mi, pim1, pi, weight=weight)

                #angle_between = np.arccos(np.dot(mutual_mi, mutual_mim1))
                #mutual_angle = get_mutual_angle_inverse(mutual_mi, mutual_mim1, angle_between)
                #mutual_m_rotated = rod_math.rodrigues(mutual_mim1, mutual_l, mutual_angle)
                if rotate:
                    mutual_m_rotated = rod_math.rodrigues(
                        mutual_m_rotated, mutual_l, rotate)
                omega, kb = rod_math.omega(pi, pim1, np.cross(
                    mutual_m_rotated, mutual_l), mutual_m_rotated)

                equil_pim1 = self.equil_p_i[frame][element]
                equil_pi = self.equil_p_i[frame][element + 1]
                equil_mim1 = self.py_rod_math.normalize(
                    self.rod.equil_m[frame][element])
                equil_nim1 = np.cross(
                    equil_mim1, equil_pim1 / np.linalg.norm(equil_pim1))
                equil_mi = self.py_rod_math.normalize(
                    self.rod.equil_m[frame][element + 1])
                equil_ni = np.cross(equil_mi, equil_pi /
                                    np.linalg.norm(equil_pi))

                # get weighting between elements adjacent to bending node!
                equil_weight = get_weights(equil_pim1, equil_pi)

                # get weighted element at that node
                equil_mutual_l = get_mutual_frame_inverse(
                    equil_pim1, equil_pi, equil_weight)

                # parallel transport our two existing material axes onto our new element
                equil_mutual_mi = np.array(rod_math.parallel_transport(
                    equil_mi, rod_math.normalize(equil_pi), equil_mutual_l))[0]
                equil_mutual_mim1 = np.array(rod_math.parallel_transport(
                    equil_mim1, rod_math.normalize(equil_pim1), equil_mutual_l))[0]

                # get the weighted average to create the mutual material axes
                equil_mutual_m_rotated = get_mutual_axes_inverse(
                    equil_mutual_mim1, equil_mutual_mi, equil_pim1, equil_pi, weight=equil_weight)

                equil_omega, equil_kb = rod_math.omega(equil_pi, equil_pim1, np.cross(
                    equil_mutual_m_rotated, equil_mutual_l), equil_mutual_m_rotated)

                delta_omega = omega - equil_omega

                inner = np.dot(np.transpose(delta_omega),
                               np.dot(B, delta_omega))

                energy = 0.5 * \
                    (inner) * (1 / ((rod_math.get_length(equil_pi) +
                                     rod_math.get_length(equil_pim1)) / 2.0))

                #L_i[frame][element] = (rod_math.get_length(equil_pi)+rod_math.get_length(equil_pim1))/(2.0)

                self.bending_energy[frame][element] = energy

                #delta_omega_all[frame][element] = np.array(omega - equil_omega).transpose()[0]

        # return delta_omega_all, L_i

    def get_twist_amount(self, set_twist_amount=False):
        """
        Get the twist energy between adjacent rod elements for every pair of
        rod elements in every frame. Mathematically, this is identical to what
        goes on inside the rod simulation proper (although it happens 8x more there).
        Params: none
        Returns: none, but it populates self.twist_energy, an array.
        """

        try:
            self.p_i
        except AttributeError:
            self.p_i = self.get_p_i(self.rod.current_r)
            self.equil_p_i = self.get_p_i(self.rod.equil_r)

        self.twist_amount = np.zeros([len(self.p_i), len(self.p_i[0]) - 1])
        self.twist_energy = np.zeros([len(self.p_i), len(self.p_i[0]) - 1])

        for frame in range(len(self.p_i)):
            for element in range(len(self.p_i[frame]) - 1):

                mim1 = self.rod.current_m[frame][element]
                mim1_equil = self.rod.equil_m[frame][element]
                mi = self.rod.current_m[frame][element + 1]
                mi_equil = self.rod.equil_m[frame][element + 1]
                pi = self.p_i[frame][element + 1]
                pim1 = self.p_i[frame][element]
                pim1_equil = self.equil_p_i[frame][element]
                pi_equil = self.equil_p_i[frame][element + 1]
                beta = self.rod.material_params[frame][element][1]
                twist_energy = self.py_rod_math.get_twist_energy(
                    mim1, mim1_equil, mi, mi_equil, pim1, pim1_equil, pi, pi_equil, beta)

                if set_twist_amount:
                    mim1_transported = rod_math.parallel_transport(rod_math.normalize(
                        mim1), rod_math.normalize(pim1), rod_math.normalize(pi))
                    delta_theta = rod_math.get_signed_angle(
                        mim1_transported, rod_math.normalize(mi), rod_math.normalize(pi))
                    self.twist_amount[frame][element] = delta_theta

                self.twist_energy[frame][element] = twist_energy

    def get_equipartition(self):
        """
        This gets all the quantities needed for an equipartition test, including
        stretching, bending and twisting energy. Basically, run this and then
        plot(temp=300) for the full equipartition experience.
        """
        self.p_i = self.get_p_i(self.rod.current_r)
        self.equil_p_i = self.get_p_i(self.rod.equil_r)

        self.get_stretch_energy()

        self.get_bending_response_mutual()

        self.get_twist_amount()

        self.average_twist_energy = np.zeros(len(self.twist_energy))
        for frame in range(len(self.average_twist_energy)):
            self.average_twist_energy[frame] = np.average(
                self.twist_energy[frame])

    def get_constant_EI(self):
        """
        This returns the value of EI (basically the value of the diagonal elements
        of B) but only if the rod is isotropic and the value of B is constant.
        Parameters: none
        Returns: none, but sets self.EI, a float.
        """
        if self.rod.B_matrix[0][0][1] != 0 and self.rod.B_matrix[0][0][2] != 0:
            raise ValueError(
                "Can't do that test for this rod. B matrix isn't diagonal, y'see.")
        for node in range(len(self.rod.B_matrix[0]) - 1):
            if self.rod.B_matrix[0][node][0] != self.rod.B_matrix[0][node][3]:
                raise ValueError(
                    "Can't do that test for this rod. B matrix isn't uniform, y'see.")
            if self.rod.B_matrix[0][node][3] != self.rod.B_matrix[0][node + 1][0]:
                raise ValueError(
                    "Can't do that test for this rod. B matrix isn't uniform, y'see.")
        self.EI = self.rod.B_matrix[0][0][0]

    def get_constant_parameter(self, constant_index):
        """
        Get the value of a material parameter based on the index (0 for stretch,
        1 for twist). Only if that parameter is constant for the whole rod and
        the whole trajectory.
        Params: constant index, an int.
        Returns: the value of that parameter in SI units, a float.
        """
        for element in range(len(self.rod.material_params[0]) - 1):
            if self.rod.material_params[0][element].all() != self.rod.material_params[0][element + 1].all():
                raise ValueError(
                    "Can't do that test for this rod. material parameters ain't uniform, y'see.")
        return self.rod.material_params[0][element][constant_index]

    def get_starting_length(self):
        """
        Get the end to end length of the rod on the first frae.
        Params: none.
        Returns: none, but sets self.starting_length, a float.
        """
        endtoend = self.rod.current_r[0][0] - self.rod.current_r[0][-1]
        self.starting_length = rod_math.get_length(endtoend)

    def get_average_quantities(self):
        """
        Sets the average quantities as they are used in the equipartition test.
        Basically identical to the equipartition test, only it also sets
        self.average_extension_sq, which we need.
        Params:none
        Returns: none, but sets self.average_extension_sq, and also
        average_extension_sq_x, average_extension_sq_y, average_extension_sq_z.
        """
        self.get_equipartition()
        #self.average_extension_sq = np.average(self.p_i_extension**2, 1)
        #self.average_extension_sq_x = np.average(self.p_i_extension_x**2, 1)
        #self.average_extension_sq_y = np.average(self.p_i_extension_y**2, 1)
        #self.average_extension_sq_z = np.average(self.p_i_extension_z**2, 1)

    def get_stretch_energy(self):
        """
        Computes the stretch energy for the rod using the new discretisation-
        independent formula.
        Params: none
        Returns: none, but sets self.stretch_energy
        """
        try:
            self.p_i
        except AttributeError:
            self.p_i = self.get_p_i(self.rod.current_r)
            self.equil_p_i = self.get_p_i(self.rod.equil_r)

        self.stretch_energy = np.zeros([len(self.p_i), len(self.p_i[0])])

        for frame_no in range(len(self.stretch_energy)):
            for element_no in range(len(self.p_i[frame_no])):
                k = self.get_constant_parameter(
                    0) / rod_math.get_length(self.equil_p_i[frame_no][element_no])
                delta_e = rod_math.get_length(
                    self.p_i[frame_no][element_no]) - rod_math.get_length(self.equil_p_i[frame_no][element_no])
                self.stretch_energy[frame_no][element_no] = 0.5 * \
                    k * (delta_e**2)

    def plot(self, force=None, temp=None):
        """
        For any tests you've done, this will plot the results.
        If you've run self.get_average_quantities(), it will plot the
        equipartition results. If you've run get_persistence_length, it plots
        those results, if you've run get_deflection, it plots that.
        The plots are saved as PDFs whose names are based on the filename of
        the loaded rod.
        Params:
            force, the force applied to the end node in the bending beam test
            (a float)
            temp, the temperature of the rod in kelvin (needed for the equipartition
            and persistence length tests.
        """

        def do_a_plot(xlabel, ylabel, xthing, ything, analytical_thing, filename):
            try:
                plot_fig = plt.figure()
                plot_ax = plot_fig.add_subplot(1, 1, 1)
                plot_ax.plot(xthing, ything)
                plot_ax.set_ylabel(ylabel)
                plot_ax.set_xlabel(xlabel)
                if analytical_thing != None:
                    plot_ax.axhline(y=analytical_thing,
                                    color='r', linestyle='--')
                plot_fig.savefig(filename, bbox_inches='tight')
                plt.close(plot_fig)
            except AttributeError:
                pass

        figname = self.rod.filename.split(".")[0]

        try:
            self.analytical_kbT = temp * 1.38064852 * 10**-23
        except TypeError:
            self.analytical_kbT = None

        self.get_constant_EI()
        self.get_starting_length()

        try:
            self.analytical_deflection = rod_math.get_analytical_deflection(
                force, self.get_absolute_length(1, self.rod.num_elements - 1, 0), self.EI)
        except TypeError:
            self.analytical_deflection = None

        try:
            half_kbT = self.analytical_kbT / 2
        except TypeError:  # if analytical_kbt = None
            half_kbT = None

        try:
            self.analytical_P = rod_math.get_analytical_persistence_length(
                self.EI, temp)
        except TypeError:
            self.analytical_P = None

        try:
            do_a_plot('Simulation steps', 'Persistence length (m)', range(len(self.rod.current_r))[
                      6:], self.P[6:], self.analytical_P, figname + "_persistence.pdf")
        except AttributeError:
            pass

        try:
            do_a_plot('Simulation steps', 'Deflection amount (m)', range(len(self.rod.current_r))[
                      6:], self.deflections[6:], self.analytical_deflection, figname + "_deflection.pdf")
        except AttributeError:
            pass

        try:
            do_a_plot(r'Simulation steps', r' $\frac{1}{2} k \langle x^2 \rangle $ (J)', range(len(
                self.rod.current_r)), self.stretch_energy, half_kbT, figname + "_equipartition_stretch.pdf")
        except AttributeError:
            pass

        try:
            do_a_plot(r'Simulation steps', r' $\frac{1}{2} k \langle x^2 \rangle $ (J)', range(len(
                self.rod.current_r)), self.average_twist_energy, half_kbT, figname + "_equipartition_twist.pdf")
        except AttributeError:
            pass

        try:
            do_a_plot(r'Simulation steps', r' $\frac{1}{2} B \langle \omega^2 \rangle $ (J)', range(len(
                self.rod.current_r)), self.whole_rod_avg(self.bending_energy), half_kbT, figname + "_equipartition_bend.pdf")
        except AttributeError:
            pass

    def persistence_segments(self):
        """
        Todo: this is an alternative version of the persistence length calulation
        that omits some of the worm-like chain assumptions.
        """
        try:
            self.p_i
        except AttributeError:
            self.p_i = self.get_p_i(self.rod.current_r)
            self.equil_p_i = self.get_p_i(self.rod.equil_r)

        s_p_fjc = np.zeros([self.rod.num_frames, self.rod.num_elements - 2])
        s_p_wlc = np.zeros([self.rod.num_frames, self.rod.num_elements - 2])
        cos_angle = np.zeros([self.rod.num_frames, self.rod.num_elements - 2])
        for frame in range(self.rod.num_frames):
            for elem in range(self.rod.num_elements - 2):
                cos_angle[frame][elem] = rod_math.get_cos_theta(
                    self.p_i[frame][elem], self.p_i[frame][elem + 1])
                s_p_fjc[frame][elem] = -1 / (np.log(cos_angle[frame][elem]))
                s_p_wlc[frame][elem] = 2 / \
                    (np.arccos(cos_angle[frame][elem])**2)
        s_p_fjc_alt = -1 / (np.log(np.average(cos_angle, axis=0)))
        s_p_wlc_alt = 2 / (np.average(np.arccos(cos_angle), axis=0)**2)
        return s_p_fjc, s_p_wlc, s_p_fjc_alt, s_p_wlc_alt

        """
   _____                                _
  / ____|                              | |
 | (___  _   _ _ __ ___  _ __ ___   ___| |_ _ __ _   _
  \___ \| | | | '_ ` _ \| '_ ` _ \ / _ \ __| '__| | | |
  ____) | |_| | | | | | | | | | | |  __/ |_| |  | |_| |
 |_____/ \__, |_| |_| |_|_| |_| |_|\___|\__|_|   \__, |
          __/ |                                   __/ |
         |___/                                   |___/

        """

    def do_stretch_symmetry_test(self):
        """
        For a rod created with the symmetry test function in rod_tests, this
        checks that the energies and dynamics are symmetric about the central
        node, for stretch energy.
        Params:none
        returns: True if all tests pass, False otherwise, and prints some info.
        """
        passed_all_tests = True

        # ENERGIES\LOOKUPS

        if not rod_math.approximately_equal(self.rod.unperturbed_energy_type[1][0][1], self.rod.unperturbed_energy_type[1][4][1], 0.001):
            print("End stretch energies are not symmetric.")
            print("Here's some debug info:")
            print("Node energy at first node: " +
                  str(self.rod.unperturbed_energy_type[1][0][1]))
            print("Node energy at last node: " +
                  str(self.rod.unperturbed_energy_type[1][4][1]))
            passed_all_tests = False

        if not rod_math.approximately_equal(self.rod.unperturbed_energy_type[1][1][1], self.rod.unperturbed_energy_type[1][3][1], 0.001):
            print("Opposite stretch energies are not symmetric.")
            print("Here's some debug info:")
            print("Node energy at second node: " +
                  str(self.rod.unperturbed_energy_type[1][0][1]))
            print("Node energy at second to last node: " +
                  str(self.rod.unperturbed_energy_type[1][4][1]))
            passed_all_tests = False

        # DYNAMICS

        if not rod_math.approximately_equal(self.rod.current_r[1][2][0], self.rod.current_r[2][2][0], 0.001):
            print(
                "The node in the middle moves! So the stretch dynamics are not symmetric.")
            print("Here's some debug info:")
            print("Initial x position of the central node: " +
                  str(self.rod.current_r[1][2]))
            print("Final x position of the central node: " +
                  str(self.rod.current_r[2][2]))

        self.p_i = self.get_p_i(self.rod.current_r)
        for element_no in range(self.rod.num_elements - 1):
            if np.dot(rod_math.normalize(self.p_i[2][element_no]), [1, 0, 0]) < 0.999:
                print("The rod doesn't stay straight.")
                print("Here's some debug info:")
                print("Rod state: " + str(self.rod.p_i[2][0]))
                passed_all_tests = False
                break

        leftmost_delta_n = self.rod.current_r[2][0] - self.rod.current_r[1][0]
        rightmost_delta_n = self.rod.current_r[2][4] - self.rod.current_r[1][4]
        dot_product = np.dot(rod_math.normalize(
            leftmost_delta_n), rod_math.normalize(rightmost_delta_n))
        if dot_product > 0.0001:
            print("The nodes at the far ends are not travelling in opposite directions.")
            print("Here's some debug info:")
            print("Leftmost node delta n: " + str(leftmost_delta_n))
            print("Rightmost node delta n: " + str(rightmost_delta_n))
            passed_all_tests = False

        return passed_all_tests

    def do_bend_symmetry_test(self):
        """
        For a rod created with the symmetry test function in rod_tests, this
        checks that the energies and dynamics are symmetric about the central
        node, for bend energy.
        Params:none
        returns: True if all tests pass, False otherwise, and prints some info.
        """
        passed_all_tests = True

        # ENERGIES\LOOKUPS

        symmetric_energies = rod_math.approximately_equal(
            self.rod.unperturbed_energy_type[1][1][1], self.rod.unperturbed_energy_type[1][3][1], 0.01)
        if not symmetric_energies:
            print("Energies are not symmetric - the energies are not the same on the let and right sides of the rod, even though it's reflected through the middle.")
            print("Here's some extra debug info:")
            print("Left energy (node 1): " +
                  str(self.rod.unperturbed_energy_type[1][1][1]))
            print("Right energy (node 3): " +
                  str(self.rod.unperturbed_energy_type[1][3][1]))
            passed_all_tests = False
        else:
            print("Passed symmetric energy tests.")

        # DYNAMICS
        #middle_angle = np.arccos(rod_math.get_cos_theta(self.rod.current_x[1][2] - self.rod.bar_x[0][2], [0,1,0] ))

        middle_dot = np.dot((rod_math.normalize(
            self.rod.current_r[1][2] - self.rod.current_r[0][2])), [0, -1, 0])

        if middle_dot < 0.999:
            print(
                "Failed bending dynamics symmetry test. The central node is not going straight down.")
            print("Here's some extra debug info:")
            print("Delta x = " +
                  str(self.rod.current_r[1][2] - self.rod.current_r[0][2]))
            passed_all_tests = False
        else:
            print("Passed central node dynamics test.")

        left = self.rod.current_r[2][1] - self.rod.current_r[1][1]
        right = self.rod.current_r[2][3] - self.rod.current_r[1][3]
        side_dynamics = rod_math.approximately_equal(
            left[0], right[0] * -1, 0.01)

        if side_dynamics == False:
            print("Failed bending dynamics symmetry test. 1st and 3rd nodes are not moving in opposite directions, even though they are at opposite positions.")
            print("Here's some extra debug info:")
            print("1st vector: " + str(left))
            print("3rd vector: " + str(right))
            passed_all_tests = False
        else:
            print("Passed reflection dynamics tests.")

        return passed_all_tests

        """
  _    _ _   _ _ _ _
 | |  | | | (_) (_) |
 | |  | | |_ _| |_| |_ _   _
 | |  | | __| | | | __| | | |
 | |__| | |_| | | | |_| |_| |
  \____/ \__|_|_|_|\__|\__, |
                        __/ |
                       |___/
        """

    def get_delta_r(self):
        """
        Get the value of delta_r from the rod dynamics by subtracting the
        positions of the nth nodes from the nth plus 1th nodes.
        Params: none
        Returns: none, but sets delta_r_x, delta_r_y, delta_r_z, delta_twist,
        all arrays.
        """
        self.delta_r_x = np.zeros(
            [len(self.rod.current_r), len(self.rod.current_r[0])])
        self.delta_r_y = np.zeros(
            [len(self.rod.current_r), len(self.rod.current_r[0])])
        self.delta_r_z = np.zeros(
            [len(self.rod.current_r), len(self.rod.current_r[0])])
        self.delta_twist = np.zeros(
            [len(self.rod.current_r), len(self.rod.current_r[0])])
        for frame in range(len(self.rod.current_r) - 1):
            for node in range(len(self.rod.current_r[frame])):
                self.delta_r_x[frame][node] = -1 * (
                    self.rod.current_r[frame][node][0] - self.rod.current_r[frame + 1][node][0])
                self.delta_r_y[frame][node] = -1 * (
                    self.rod.current_r[frame][node][1] - self.rod.current_r[frame + 1][node][1])
                self.delta_twist[frame][node] = np.dot(rod_math.normalize(
                    self.rod.current_m[frame][node]), rod_math.normalize(self.rod.current_m[frame + 1][node]))

    def get_L_i(self):
        """
        """
        try:
            self.equil_p_i
        except AttributeError:
            self.equil_p_i = self.get_p_i(self.rod.equil_r)

        L_i = np.zeros(
            [len(self.rod.current_r), len(self.rod.current_r[0]) - 1])
        for frame_no in range(len(L_i)):
            for index in range(len(L_i[0])):
                L_i[frame_no][index] = 0.5 * (rod_math.get_length(
                    self.equil_p_i[index + 1]) + rod_math.get_length(self.equil_p_i[index]))
        return L_i
        # todo!!!
        # note: L_i[0][0] is the L_i for the first node, but due to the way this is indexed this is actually L1 not L0.

    def get_absolute_length(self, starting_node, ending_node, frame, equil=False, path=True, non_path_axis=0):
        """
        Get the length along the rod (the contour length) of a subset of the
        rod, or the whole rod, in meters.
        Params: starting noe, ending node - indices of the start and end
        of the range of nodes to be measured.
        Frame: index of the frame to be used.
        Returns: length, a float.
        """

        if equil:
            try:
                self.equil_p_i
            except AttributeError:
                self.equil_p_i = self.get_p_i(self.rod.equil_r)
            p_i = self.equil_p_i
        else:
            try:
                self.p_i
            except AttributeError:
                self.p_i = self.get_p_i(self.rod.current_r)
            p_i = self.p_i

        if path:
            absolute_length = 0
            for i in range(ending_node - starting_node):
                absolute_length += rod_math.get_length(
                    p_i[frame][i + starting_node])
        else:
            absolute_length = self.rod.current_r[frame][starting_node][non_path_axis] - \
                self.rod.current_r[frame][ending_node][non_path_axis]
        return abs(absolute_length)

    def whole_rod_avg(self, thing):
        """
        Given an array of form [frame][element][thing], this function returns
        the average [thing] over every element for each frame.
        Params: thing, an array.
        Returns: the average over the whole rod, also an array.
        """
        return np.average(thing, 1)

    def filter_x(self, thing, dimension):
        """
        Given an array of form [frame][element][x,y,z], this returns an array
        of the form [frame][element][thing], where thing is the quantitiy for
        the dimension specified as a parameter.
        Params: thing, the array to filter, and dimension, the dimension to
        filter to.
        Returns: the filtered array.
        """
        return thing[:, :, dimension]

    def apply_transformation_4x4(self, pos_array, T):
        """
        Apply a 4x4 transformation matrix to our 2-D array of 3-D points. Returns
        the newly translated array.
        """
        translated_array = np.zeros([len(pos_array), 3])
        for i in range(len(pos_array)):
            node_1 = [pos_array[i][0], pos_array[i][1], pos_array[i][2], 1]
            translated_array[i] = np.dot(T, node_1)[:3]
        return translated_array

    def align_to_equil(self):
        """
        Align the trajectory to the equilibrium configuration using the icp
        (iterative closest point) library.
        Params: none
        Retruns: none, but updates self.rod.current_r
        """
        try:
            import icp
        except ImportError:
            raise ImportError(
                "Please install the icp (iterative closest point) library!")

        for frame_index in range(len(self.rod.current_r)):  # aargh
            T, distances = icp.icp(
                self.rod.current_r[frame_index], self.rod.equil_r[frame_index], max_iterations=10000, tolerance=1e-9)
            self.rod.current_r[frame_index] = self.apply_transformation_4x4(
                self.rod.current_r[frame_index], T)

            if frame_index % 1000 == 0:
                print("Aligning frame " + str(frame_index) + "...")

        return

    def thin(self, target_num_frames):
        """
        Remove frames from the trajectory.
        Params: target_num_frames - the number of frames you desire.
        Returns: nothing, but it changes all the self.rod arrays.
        Also, like decimate, it might not give you exactly the number
        of frames you asked for, because it goes for an even interval
        instead.
        """
        interval = int(len(self.rod.current_r) / target_num_frames)

        self.rod.equil_r = self.rod.equil_r[::interval]
        self.rod.equil_m = self.rod.equil_m[::interval]
        self.rod.current_r = self.rod.current_r[::interval]
        self.rod.current_m = self.rod.current_m[::interval]
        self.rod.internal_perturbed_x_energy_positive = self.rod.internal_perturbed_x_energy_positive[
            ::interval]
        self.rod.internal_perturbed_y_energy_positive = self.rod.internal_perturbed_y_energy_positive[
            ::interval]
        self.rod.internal_perturbed_z_energy_positive = self.rod.internal_perturbed_z_energy_positive[
            ::interval]
        self.rod.twisted_energy_positive = self.rod.twisted_energy_positive[::interval]
        self.rod.internal_perturbed_x_energy_negative = self.rod.internal_perturbed_x_energy_negative[
            ::interval]
        self.rod.internal_perturbed_y_energy_negative = self.rod.internal_perturbed_y_energy_negative[
            ::interval]
        self.rod.internal_perturbed_z_energy_negative = self.rod.internal_perturbed_z_energy_negative[
            ::interval]
        self.rod.twisted_energy_negative = self.rod.twisted_energy_negative[::interval]
        self.rod.material_params = self.rod.material_params[::interval]
        self.rod.B_matrix = self.rod.B_matrix[::interval]
        self.rod.steric_energy = self.rod.steric_energy[::interval]
        self.rod.steric_force = self.rod.steric_force[::interval]
        self.rod.num_neighbours = self.rod.num_neighbours[::interval]
        self.rod.num_frames = len(self.rod.current_r)

        try:
            self.p_i = self.p_i[::interval]
        except AttributeError:
            pass
        try:
            self.equil_p_i = self.equil_p_i[::interval]
        except AttributeError:
            pass

    def get_node_rmsd(self, align=False):
        """
        Get the per-node RMSD for the rod. Otherwise known as the time-averaged
        RMSD for each node.
        Params: align - whether the trajectory has already been aligned to the
        equilibrium configuration. If not, it uses the iterative closest point
        algorithm to align them. Please note: this algorithm can get confused
        with very bendy rods. If you see a gigantic spike in the middle of your
        RMSD, then it's confused (easier to spot in get_time_rmsd).
        Returns: the rmsd, a 1d numpy array indexed by node.
        """

        rmsd = np.zeros(len(self.rod.current_r[0]))

        if align:
            try:
                import icp
            except ImportError:
                raise ImportError(
                    "Please install the icp (iterative closest point) library!")

        if align:
            for frame_index in range(len(self.rod.current_r)):  # aargh
                if frame_index % 15000 == 0:
                    # note: icp algorithm is confused if it has to start from scratch all the time,
                    # but we reset it every so often to mitigate against floating point errors
                    T, distances = icp.icp(
                        self.rod.equil_r[frame_index], self.rod.current_r[frame_index], max_iterations=10000, tolerance=1e-9)
                    self.rod.equil_r[frame_index] = self.apply_transformation_4x4(
                        self.rod.equil_r[frame_index], T)
                else:
                    T, distances = icp.icp(
                        self.rod.equil_r[frame_index - 1], self.rod.current_r[frame_index], max_iterations=10000, tolerance=1e-9)
                    self.rod.equil_r[frame_index] = self.apply_transformation_4x4(
                        self.rod.equil_r[frame_index - 1], T)

                if frame_index % 1000 == 0:
                    print("Aligning frame " + str(frame_index) + "...")
        #N = len(equil)

        for node_no in range(len(rmsd)):

            current = self.rod.current_r[:, node_no]
            equil = self.rod.equil_r[:, node_no]

            rmsd[node_no] = np.sqrt(np.average(np.linalg.norm(
                current - equil, axis=1)**2))  # need get length!!!!

        return rmsd

    def get_B_eigenvalues(self):
        """
        Get the eigenvalues of the B matrix for each node.
        Params: none
        Returns: an array of the eigenvalues. Each B value will have two
        eigenvalues, so the resulting array will be (n-2)x2, where n is the
        number of nodes (-2 because the nodes at the end of the rod do not) have
        bending energies associated with them!).
        """
        evals = []
        for B_matrix in self.rod.B_matrix[0]:
            eigval, evecs = np.linalg.eig(
                np.matrix(np.reshape(B_matrix, [2, 2])))
            evals.append(eigval)
        return np.array(evals)

    def get_time_rmsd(self, is_aligned=False, max_frame_index=-1):
        """
        Get the time RMSD for the rod. That is, get the average RMSD of all the
        nodes for each frame in the trajectory.
        Params:
            is_aligned: whether the equilibrium and current structures are
            aligned or not. If they aren't aligned, they need to be, or you
            won't get best-fit RMSD, which is really the only useful RMSD.
            max_frame_index: find the rmsd for all frames up to this one.
        Returns:
            the time rmsd, a 1d numpy array.
        """

        if not is_aligned:
            self.get_node_rmsd(align=True)

        current_r = self.rod.current_r[:max_frame_index]
        equil_r = self.rod.equil_r[:max_frame_index]

        rmsds = np.zeros(len(equil_r))

        for frame in range(len(rmsds)):
            rmsds[frame] = np.sqrt(np.average(np.linalg.norm(
                current_r[frame] - equil_r[frame], axis=1)**2))

        return rmsds

    def subdivide(self, iterations):
        """
        For a given rod (but not a trajectory!) create new nodes which are
        linearly interpolated to be between the positions of the existing nodes.
        Params: iterations - number of iterations to run. Each iteration
        doubles the number of nodes in the rod.
        Returns: nothing, but updates all the self.rod ndarrays.
        """
        def interp_r(node1, node2):
            return np.average([node1, node2], axis=0)

        if self.rod.num_frames > 1:
            print("Yo, don't subdivide a trajectory!")
            return
        for i in range(iterations):
            element_index = range(self.rod.num_elements)
            element_index.reverse()
            element_index = element_index[0:-1]

            for element in element_index:
                self.rod.equil_r = np.array([np.insert(self.rod.equil_r[0], element, interp_r(
                    self.rod.equil_r[0][element], self.rod.equil_r[0][element - 1]), axis=0)])
                self.rod.equil_m = np.array([np.insert(self.rod.equil_m[0], element, interp_r(
                    self.rod.equil_m[0][element], self.rod.equil_m[0][element - 1]), axis=0)])
                self.rod.current_r = np.array([np.insert(self.rod.current_r[0], element, interp_r(
                    self.rod.current_r[0][element], self.rod.current_r[0][element - 1]), axis=0)])
                self.rod.current_m = np.array([np.insert(self.rod.current_m[0], element, interp_r(
                    self.rod.current_m[0][element], self.rod.current_m[0][element - 1]), axis=0)])
                self.rod.internal_perturbed_x_energy_positive = np.array(
                    [np.insert(self.rod.internal_perturbed_x_energy_positive[0], element, 0, axis=0)])
                self.rod.internal_perturbed_y_energy_positive = np.array(
                    [np.insert(self.rod.internal_perturbed_y_energy_positive[0], element, 0, axis=0)])
                self.rod.internal_perturbed_z_energy_positive = np.array(
                    [np.insert(self.rod.internal_perturbed_z_energy_positive[0], element, 0, axis=0)])
                self.rod.twisted_energy_positive = np.array(
                    [np.insert(self.rod.twisted_energy_positive[0], element, 0, axis=0)])
                self.rod.internal_perturbed_x_energy_negative = np.array(
                    [np.insert(self.rod.internal_perturbed_x_energy_negative[0], element, 0, axis=0)])
                self.rod.internal_perturbed_y_energy_negative = np.array(
                    [np.insert(self.rod.internal_perturbed_y_energy_negative[0], element, 0, axis=0)])
                self.rod.internal_perturbed_z_energy_negative = np.array(
                    [np.insert(self.rod.internal_perturbed_z_energy_negative[0], element, 0, axis=0)])
                self.rod.twisted_energy_negative = np.array(
                    [np.insert(self.rod.twisted_energy_negative[0], element, 0, axis=0)])
                self.rod.material_params = np.array([np.insert(self.rod.material_params[0], element, interp_r(
                    self.rod.material_params[0][element], self.rod.material_params[0][element - 1]), axis=0)])
                self.rod.B_matrix = np.array([np.insert(self.rod.B_matrix[0], element, interp_r(
                    self.rod.B_matrix[0][element], self.rod.B_matrix[0][element - 1]), axis=0)])
            self.rod.num_elements = len(self.rod.equil_r[0])
            self.rod.length = 3 * self.rod.num_elements
        return

    def decimate(self, determine_simplification_func, target_length, margin=0.5):
        """
        For a given rod or rod trajectory, reduce the number of nodes in that
        trajectory, by averaging about nodes, grouped by their indices.
        Params:
            determine_simplification_func: the function that reduces the number
            of nodes in an ndarray, normally ndc_extractor.determine_simplification
            target_length: the desired number of nodes. Note: the actual
            rod you get out might be a number nearby this but not exactly it,
            to try and make the step size (e.g. the number of nodes that are
            averaged together) even.
            margin: how much to adjust the target length to make the step size
            even.
        Returns: noting. But updates all the self.rod arrays.
        """
        # note: use the determine_simplification function from ndc_extractor.py
        self.rod.equil_r = determine_simplification_func(
            self.rod.equil_r, target_length, margin)
        self.rod.equil_m = determine_simplification_func(
            self.rod.equil_m, target_length, margin)
        self.rod.current_r = determine_simplification_func(
            self.rod.current_r, target_length, margin)
        self.rod.current_m = determine_simplification_func(
            self.rod.current_m, target_length, margin)
        self.rod.internal_perturbed_x_energy_positive = determine_simplification_func(
            self.rod.internal_perturbed_x_energy_positive, target_length, margin)
        self.rod.internal_perturbed_y_energy_positive = determine_simplification_func(
            self.rod.internal_perturbed_y_energy_positive, target_length, margin)
        self.rod.internal_perturbed_z_energy_positive = determine_simplification_func(
            self.rod.internal_perturbed_z_energy_positive, target_length, margin)
        self.rod.twisted_energy_positive = determine_simplification_func(
            self.rod.twisted_energy_positive, target_length, margin)
        self.rod.internal_perturbed_x_energy_negative = determine_simplification_func(
            self.rod.internal_perturbed_x_energy_negative, target_length, margin)
        self.rod.internal_perturbed_y_energy_negative = determine_simplification_func(
            self.rod.internal_perturbed_y_energy_negative, target_length, margin)
        self.rod.internal_perturbed_z_energy_negative = determine_simplification_func(
            self.rod.internal_perturbed_z_energy_negative, target_length, margin)
        self.rod.twisted_energy_negative = determine_simplification_func(
            self.rod.twisted_energy_negative, target_length, margin)
        self.rod.material_params = determine_simplification_func(
            self.rod.material_params, target_length, margin)
        self.rod.B_matrix = determine_simplification_func(
            self.rod.B_matrix, target_length, margin)
        self.rod.steric_energy = determine_simplification_func(
            self.rod.steric_energy, target_length, margin)
        self.rod.steric_force = determine_simplification_func(
            self.rod.steric_force, target_length, margin)
        self.rod.num_elements = len(self.rod.equil_r[0])
        self.rod.length = 3 * self.rod.num_elements


"""
  __  __       _   _
 |  \/  |     | | | |
 | \  / | __ _| |_| |__  ___
 | |\/| |/ _` | __| '_ \/ __|
 | |  | | (_| | |_| | | \__ \
 |_|  |_|\__,_|\__|_| |_|___/

"""


class py_rod_math:
    """
    This object contains pure functions relating to mathematical formulas for
    the rods, simular to rod_math.cpp.
    """

    def get_twist_angle(self, m_i, m_j):
        """
        Get the angle between two vectors, assumed to have been parallel
        tranported into the same basis.
        Params: m_i, m_j, two vectors
        Returns: the angle between them.
        """
        return np.arccos(np.dot(self.normalize(m_i), self.normalize(m_j)))

    def get_signed_angle(self, m1, m2, l):
        """
        Same as the above, but the angle is signed (left hand rotation).
        Params: m1, m2 - the two material axis vectors being measured.
        l: the normalized element vector.
        returns: the angle between them, in radians.
        """
        # https://stackoverflow.com/a/33920320
        return np.arctan2(np.dot(np.cross(m2, m1), l), np.dot(m1, m2))

    def parallel_transport(self, m_i, a, b):
        """
        Given one of the material frame vectors m_i, and two vectors, a, and b,
        this function transports the material frame from the basis of a to the
        basis of b.
        """
        v = np.cross(a, b)
        c = np.dot(a, b)
        vx = np.matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        c_factor = 1. / (1 + c)
        R = np.identity(3) + vx + np.linalg.matrix_power(vx, 2) * c_factor
        m_i_prime = np.dot(R, m_i)
        return m_i_prime

    def get_length(self, vec):
        """
        Given a vector, vec, return the length of the vector, a float.
        """
        return np.linalg.norm(vec)

    def normalize(self, vec):
        """
        Given a vector, vec, return the normalized version of that vector, a
        vector (a 1-D numpy array containing 3 values.)
        """
        return vec / np.linalg.norm(vec)

    def get_cos_theta(self, a, b):
        """
        For two vectors, a and b, get the cosine of the angle between them.
        Returns a float.
        """
        return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))

    def get_dist(self, a, b):
        """
        Get the length of the vector between two vectors a and b. Returns a
        float.
        """
        distvec = a - b
        return self.get_length(distvec)

    def kb_i(self, pi, pim1):
        """
        For p_i and p_im1, return the value of the vector kb_i. For more
        information on this, please check the bend energy functions in rod_math.
        """
        return np.cross(2 * pim1, pi) / (np.linalg.norm(pi) * np.linalg.norm(pim1) + np.dot(pi, pim1))

    def omega(self, pi, pim1, m2j, m1j):  # n, m
        """
        For for two elements and two material frame bectors, return the 2-vector
        omega. For more information on this, please check the bend energy functions
        in rod_math.
        """
        kb_i = self.kb_i(pi, pim1)
        omega_horiz = np.matrix([np.dot(kb_i, m2j), -np.dot(kb_i, m1j)])
        return np.transpose(omega_horiz), kb_i

    def get_analytical_deflection(self, P, L, EI):
        """
        Get the excepcted deflection W of the rod based on P, the force acting
        upon the end node, L, the length, and EI (B), the bending modulus.
        W = (P*L^)3/(3*EI)
        """
        return (P * (L**3)) / (3 * EI)

    def get_analytical_persistence_length(self, EI, T):
        """
        Get the expected value of the persistence length based on the temperature
        T and the bending modulus EI (B).
        P = EI/(k_B*T)
        """
        k_B = 1.38064852 * 10**-23
        return EI / (k_B * T)

    def approximately_equal(self, a, b, factor_of_error):
        """
        Given two values and a margin of error, check if the two values are
        within that error of one another.
        Params: a, b, margin
        Returns: boolean, true if they are within the margin, false otherwise.
        """
        margin_of_error = factor_of_error * np.average([np.abs(a) + np.abs(b)])
        a_max = a + margin_of_error
        a_min = a - margin_of_error
        if b < a_min or b > a_max:
            return False
        elif b == a:
            return True
        else:
            return True

    def rodrigues(self, v, k, theta):
        """
        A python version of the rogrigues rotation of a vector V about an axis
        K, by an angle theta.
        Params: v, the vector (list, numpy, whatever), k, a vector, theta, a float
        Return: v_prime, the rotated vector (a numpy array)
        """
        k = self.normalize(k)
        try:
            return np.cos(theta) * v + np.cross(k, v) * np.sin(theta) + k * np.dot(k, v) * (1 - np.cos(theta))
        except TypeError:
            return np.cos(theta) * np.array(v) + np.cross(k, v) * np.sin(theta) + np.array(k) * np.dot(k, v) * (1 - np.cos(theta))

    def get_stretch_energy(self, k, p_i, p_i_equil):
        """
        Get the stretch energy for a particular element.
        Params:
            k, the stretching constant (not the spring constant!)
            p_i, a 1-d array of 3 points, the p_ith element
            p_i_equil, the equilibrium p_ith element
        Returns:
            the stretch energy.
        """
        diff = self.get_length(p_i) - self.get_length(p_i_equil)

        stretch_energy = (diff * diff * 0.5 * k) / self.get_length(p_i_equil)

        return stretch_energy

    def get_twist_energy(self, mim1, mim1_equil, mi, mi_equil, pim1, pim1_equil, pi, pi_equil, beta):
        """
        Get the energy associated with a particular twist. Note: this uses the
        new version of the twist energy calculation that allows for intrinsically
        twisted equilibrium configurations.
        Params:
            mim1: 1-d numpy array with xyz co-ordinates for mim1, the i-1th material axis.
            mim1_equil: the equilibrium configuration of the i-1th material axis.
            mi: the ith material axis
            mi_equil: the equilibrium ith material axis
            pim1: 1-d numpy array with xyz co-ordinates for pim1, the i-1th element.
            pim1_equil: the equilibrium configuration of the i-1th element.
            pi: the ith element
            pi_equil: the equilibrium ith element
            beta: a float, the twisting constant
        Returns:
            the twist energy
        """
        mim1_transported = rod_math.parallel_transport(rod_math.normalize(
            mim1), rod_math.normalize(pim1), rod_math.normalize(pi))
        delta_theta = rod_math.get_signed_angle(
            mim1_transported, rod_math.normalize(mi), rod_math.normalize(pi))
        equil_mim1_transported = rod_math.parallel_transport(rod_math.normalize(
            mim1_equil), rod_math.normalize(pim1_equil), rod_math.normalize(pi_equil))
        equil_delta_theta = rod_math.get_signed_angle(
            equil_mim1_transported, rod_math.normalize(mi_equil), rod_math.normalize(pi_equil))
        Li = (rod_math.get_length(pim1_equil) +
              rod_math.get_length(pi_equil)) / 2
        twist_energy = beta / (2 * Li) * np.power(np.mod(delta_theta -
                                                         equil_delta_theta + np.pi, 2 * np.pi) - np.pi, 2)
        return float(twist_energy)

    def get_bend_energy(self, p_im1, p_i, p_im1_equil, p_i_equil, n_im1, m_im1, n_im1_equil, m_im1_equil, n_i, m_i, n_i_equil, m_i_equil, B_i_equil, B_im1_equil):
        raise Exception("Not implemented yet!")

    def perpendicularize(self, m_i, p_i):
        """
        Make two vectors perpendicular. Usually used to make the material axes and elements parallel.
         \f[\widetilde{m_{1 i}}' = \widetilde{m_{1 i}} - ( \widetilde{m_{1 i}} \cdot \widetilde{l_i}) \widetilde{\hat{l_i}}\f]
         where \f$l\f$ is the normalized tangent, \f$m\f$ is the current material frame and \f$m'\f$ is the new one.
        Params:
             m_i: material axis (to be perpendicularized), a 1d array
             p_i: element, a 1d array
        Returns:
            m_i_prime, the perpendicularized material axis, 1 1d array.
        """
        t_i = self.normalize(p_i)
        m_i_dot_t_i = np.dot(m_i, t_i)
        m_i_prime = m_i - m_i_dot_t_i * t_i
        return m_i_prime

    def get_euler_angles(self, rm):

        if rm[6] is not 1 and rm[6] is not -1:
            theta_1 = -np.arcsin(rm[6])
            psi_1 = np.arctan2(rm[7] / np.cos(theta_1),
                               rm[8] / np.cos(theta_1))
            phi_1 = np.arctan2(rm[3] / np.cos(theta_1),
                               rm[0] / np.cos(theta_1))

        else:
            print("Gimbal lock! Nice!")
            if rm[6] == -1:
                theta_1 = np.pi / 2
                psi_1 = phi_1 + np.arctan2(rm[1], rm[2])
            else:
                theta_1 = -np.pi / 2
                psi_1 = -phi_1 + np.arctan2(-rm[1], -rm[2])

        return np.array([phi_1, theta_1, psi_1])


class fast_rod_math(py_rod_math):
    def __init__(this):
        this.get_strech_energy = rod_math_core.get_stretch_32
        this.get_bend_energy = rod_math_core.get_bend_32
        this.get_twist_energy = rod_math_core.get_twist_32


try:
    import rod.rod_math_core as rod_math_core
    rod_math = fast_rod_math()
    print("Using fast math functions.")
    #py_rod_math.get_twist_energy = rod_math_core.get_twist_32
    #py_rod_math.get_stretch_energy = rod_math_core.get_stretch_32
    #py_rod_math.get_bend_energy = rod_math_core.get_bend_32
    rod_math_core_status = True
except Exception as e:
    print("Could not import fast math functions from c. Using slow ones instead. Hey, check this error!")
    print(e)
    rod_math = py_rod_math()
    rod_math_core_status = False
    pass


test_threejs = """[new THREE.Vector3(289.76843686945404, 452.51481137238443, 56.10018915737797),
	new THREE.Vector3(5.170639673700947, 66.80453441905567, 28.658512034683916),
	new THREE.Vector3(-206.8243653217158, 400.79421993281676, -201.7518997762312),
	new THREE.Vector3(-379.2121784101623, 450.151231810768, 137.16480548796358)]"""

"""
   _____                _
  / ____|              | |
 | |     _ __ ___  __ _| |_ ___  _ __
 | |    | '__/ _ \/ _` | __/ _ \| '__|
 | |____| | |  __/ (_| | || (_) | |
  \_____|_|  \___|\__,_|\__\___/|_|

"""


class rod_creator:
    """
    The rod_creator class is just a collection of functions that create rods,
    it doesn't have any state or anything. In fact, immediately after this
    script is loaded, the rod_creator class definition is replaced by an
    instance of itself.

    All of the creation functions here output rod frames, so in order to use
    the rod creator, you'll need to initialize an empty rod object, and then
    run either create_rod_spline or create_rod_parametric to create
    rod.equil_r[0] and rod.current_r[0].  Then, use create_material_frame and possibly
    rotate_material_frame to create rod.equil_m[0] and rod.current_m[0].
    """

    def parse_threejs(self, threejs):
        """
        Parse a threejs spline object, generated and exported via this tool:
        https://threejs.org/examples/?q=spline#webgl_geometry_spline_editor
        Parameters: a string containing the contents of the 'exportSpline' dialog.
        Returns 3 arrays containing the x, y and z co-ordinates of the spline
        control points.
        """
        rod_spline = []
        for line in threejs.split("\n"):
            line = line.replace("new THREE.Vector3(", "")
            line = line.replace("),", "")
            line = line.replace(")", "")
            line = line.replace("[", "")
            line = line.replace("]", "")
            line = line.split(", ")
            rod_spline.append([float(line[0]), float(line[1]), float(line[2])])
        # kinda crappy
        return np.array(rod_spline)[:, 0], np.array(rod_spline)[:, 1], np.array(rod_spline)[:, 2]

    def parse_pdb(self, pdb_file, atom_names=["CA"]):
        try:
            import FFEA_pdb
        except ImportError:
            try:
                import ffeatools.FFEA_pdb as FFEA_pdb
            except ImportError:
                raise ImportError(
                    "Couldn't find FFEA_pdb, is FFEA not installed?")

        pos = []
        pdb = FFEA_pdb.FFEA_pdb(pdb_file)

        if len(atom_names) == 0:  # get all atoms
            for chain in pdb.chain:
                for frame in chain.frame:
                    pos.append(frame.pos)
            return np.concatenate(pos) * 1e-10  # angrstroms to m
        else:
            atoms_to_do = []  # filter atoms based on string user gives
            for chain_no in range(len(pdb.chain)):
                for atom_no in range(len(pdb.chain[chain_no].atom)):
                    if pdb.chain[chain_no].atom[atom_no].res in atom_names or pdb.chain[chain_no].atom[atom_no].name.replace(' ', '') in atom_names:
                        atoms_to_do.append([atom_no, chain_no])
            for atoms_and_chains in atoms_to_do:
                atom = atoms_and_chains[0]
                chain = atoms_and_chains[1]
                pos.append(pdb.chain[chain].frame[0].pos[atom])
            return np.array(pos) * 1e-10  # angrstroms to m

    def preview_rod(self, rod):
        """
        Use matplotlib to plot a frame, as generated by the create_rod_spline
        or create_rod_parametric functions. Useful if you want to see what
        your rod will look like without having to write the trajectory and
        open up pymol.
        """
        fig2 = plt.figure(2)
        ax3d = fig2.add_subplot(111, projection='3d')
        ax3d.plot(rod[:, 0], rod[:, 1], rod[:, 2], 'g')
        plt.show()

    def create_spline(self, x_points, y_points, z_points, num_samples, degree=3, smoothness=0):
        """
        Given a set of spline control points and some spline parameters, make
        a spline, and sample it at a bunch of points, returning a 2-d array
        containing the x, y, and z co-ordinates of the spline at those points.

        Parameters:
            x_points, y_points, z_points: spline control points
            num_samples: number of times to sample the spline equation
            degree: degree of the polynomial
            smoothness: number of continuous derivatives of the spline function
        Returns:
            rod, a 2-d array containing the co-ordinates of the spline samples
            knots_tck - a tuple (t,c,k) containing the vector of knots, b-spline
            coefficients, and the degree of the spline
            spline_params: an array of the values of the parameter
        """
        # degree = 3 # degree if the polynomial
        # smoothness = 2 # number of continuous derivatives of the spline function

        sample_points = np.linspace(0, 1, num_samples)
        knots_tck, spline_params = interpolate.splprep(
            [x_points, y_points, z_points], s=smoothness, k=degree)
        x_interped, y_interped, z_interped = interpolate.splev(
            sample_points, knots_tck)
        rod = np.empty([num_samples, 3])
        rod[:, 0] = x_interped
        rod[:, 1] = y_interped
        rod[:, 2] = z_interped
        return rod, knots_tck, spline_params

    def create_rod_spline(self, data, num_nodes, degree=3, smoothness=2):
        """
        Given a 2-d array of spline control points (or a string containing the
        output from the threejs spline creator), the number of nodes desired, the
        smoothness and degree of the desired spline, return a rod frame.
        Use this to fill in the current_r[0] or equil_r[0] of a object.

        Params:
            data: either a 2-d array containing spline control points in 3 dimensions,
            or a string containing the output from this tool:
            https://threejs.org/examples/?q=spline#webgl_geometry_spline_editor
            num_nodes: number of times to sample the spline equation
            degree: degree of the polynomial
            smoothness: number of continuous derivatives of the spline function
        Returns:
            a rod frame, for either current_r or equil_r (a 2-d numpy array of points)
        """
        if type(data) == type(np.array([])):
            transposed = np.transpose(data)
            current_r, knots_tck, spline_params = self.create_spline(
                transposed[0], transposed[1], transposed[2], num_nodes, degree=degree, smoothness=smoothness)

        if type(data) == type("string"):
            x_points, y_points, z_points = self.parse_threejs(data)
            current_r, knots_tck, spline_params = self.create_spline(
                x_points, y_points, z_points, num_nodes, degree=degree, smoothness=smoothness)

        if callable(data):
            raise ValueError("Use create_rod_parametric!")

        return current_r

    def create_rod_parametric(self, x_eqn, y_eqn, z_eqn, t_max, t_min, num_nodes):
        """
        Create a rod frame (a 2-d array of 3-d points) using 3 parametric
        equations, one for each dimension.

        Params:
            x_eqn, y_eqn, z_eqn: python functions with one parameter (t)
            t_max, t_min: max and min values for t (min at node 0, max at n)
            num_nodes: number of nodes
        Returns:
            a 2-d numpy array containing the positions of the rod nodes,
            that you can use as rod.current_r[0] and rod.equil_r[0]
        """
        t_values = np.linspace(t_min, t_max, num=num_nodes)
        rod = np.empty([num_nodes, 3])
        for i in range(len(t_values)):
            rod[i][0] = x_eqn(t_values[i])
            rod[i][1] = y_eqn(t_values[i])
            rod[i][2] = z_eqn(t_values[i])
        return rod

    def create_material_frame(self, rod):
        """
        For a given rod, create a material frame.
        For the rod to be at equilibrium twist, all of the m and n vectors
        must be pointing the same way relative to the rod segments. So, to
        create an equilibrium material frame, you create a material frame
        arbitrarily for one element, and then parallel transport it up the
        entire length of the rod.

        Params:
            rod - a rod object with a current_r

        Returns:
            current_material_frame: material frame for the rod in its current state
            equil_material_frame: material frame for the rod in its
        """
        try:
            rod.p_i
            rod.equil_p_i
        except AttributeError:
            rod.p_i = rod.get_p_i(rod.current_r)
            rod.equil_p_i = rod.get_p_i(rod.equil_r)

        if rod.num_frames > 1:
            print(
                "Are you creating a material frame for a rod that already has a trajectory?")

        current_material_frame = np.zeros([1, rod.num_elements, 3])
        equil_material_frame = np.zeros([1, rod.num_elements, 3])

        # The selection of the first material frame is arbitrary. If our first two p_i elements are different, we can use their
        # cross product to give us a value. If not, we use a somewhat arbitrary value.
        if np.array_equal(rod_math.normalize(rod.p_i[0][0]), rod_math.normalize(rod.p_i[0][1])):
            current_first_material_frame = rod_math.normalize(np.cross(rod_math.normalize(
                rod.p_i[0][0]), np.array([0.48154341,  0.84270097, -0.24077171])))
        else:
            current_first_material_frame = rod_math.normalize(
                np.cross(rod_math.normalize(rod.p_i[0][0]), rod_math.normalize(rod.p_i[0][1])))

        if np.array_equal(rod_math.normalize(rod.equil_p_i[0][0]), rod_math.normalize(rod.equil_p_i[0][1])):
            equil_first_material_frame = rod_math.normalize(np.cross(rod_math.normalize(
                rod.equil_p_i[0][0]), np.array([0.48154341,  0.84270097, -0.24077171])))
        else:
            equil_first_material_frame = rod_math.normalize(np.cross(rod_math.normalize(
                rod.equil_p_i[0][0]), rod_math.normalize(rod.equil_p_i[0][1])))

        current_material_frame[0][0] = current_first_material_frame
        equil_material_frame[0][0] = equil_first_material_frame

        for element in range(len(rod.p_i[0])):
            if element == 0:
                continue
            current_material_frame[0][element] = rod_math.parallel_transport(
                current_material_frame[0][0], rod_math.normalize(rod.p_i[0][0]), rod_math.normalize(rod.p_i[0][element]))
            equil_material_frame[0][element] = rod_math.parallel_transport(equil_material_frame[0][0], rod_math.normalize(
                rod.equil_p_i[0][0]), rod_math.normalize(rod.equil_p_i[0][element]))

        return current_material_frame, equil_material_frame

    # todo: add a function that will parallel transform the material frames from equil_m to current_m and vice versa

    def rotate_material_frame(self, rod, function, material_frame, function_range=(0, 2 * np.pi)):
        """
        Rotate the material frame. Note: you should not use this to rotate the
        equilibrium material frame, that should always be straight. If you
        want to create an intrinsically twisted rod with an anisotropic bending
        modulus, rotate the B matrix instead!

        Params:
            rod - a rod object
            function: the function that determines what the bend angle will be
            material_frame: rod.current_m
            function_range: a tuple containing the parameter of the function at
            the beginning and end of the rod. Default is (0, 2*pi).

        Returns:
            a material frame frame.
        """
        try:
            rod.p_i
            rod.equil_p_i
        except AttributeError:
            rod.p_i = rod.get_p_i(rod.current_r)
            rod.equil_p_i = rod.get_p_i(rod.current_r)

        for i in range(rod.num_elements - 1):
            func_value = function((i / rod.num_elements - 1) * 2 * np.pi)
            material_frame[0][i] = rod_math.rodrigues(
                rod.p_i[0][i], material_frame[0][i], func_value)

        return material_frame

    def set_params(self, rod, stretch_constant, torsion_constant, radius, bending_modulus=None, rod_segments=None, E=None, r=None, frames=None):
        """
        Set the material parameters uniformally. By default, this is for
        the entire rod, although optionally you can set some sub-region of
        the rod.
        Params:
            rod: rod object.
            stretch_constant: stretch constant, in SI units (this is normally
            the spring constant multiplied by the measured spring length!)
            torsion_constnat: torsion constant, in SI units (this is the same deal)
            radius: radius of the rod (for computing the mobility). Half of the equilibrium rod length is about right.
            bending modulus: this is normally EI, the young's modulus multiplied
            by the second moment of area. Optionally, you can instead supply
            E, I: a young's modulus and a radius.
            rod_segment (optional): a range object containing the indices
            of nodes to modify.
        Returns:
            nothing, it modifies the rod object that's given to it.
        """

        def second_moment_of_area_circle(r):
            return (np.pi / 4.0) * (r**4)

        if rod_segments == None:
            rod_segments = range(0, len(rod.current_r[0]))

        if frames == None:
            frames = range(0, len(rod.current_r))

        if bending_modulus == None:
            if E or r == None:
                raise ValueError(
                    "Need either a bending modulus OR E and r, chief.")
            bending_modulus = E * second_moment_of_area_circle(r)

        for frame in frames:
            for segment in rod_segments:
                rod.material_params[frame][segment] = np.array(
                    [stretch_constant, torsion_constant, radius])
                rod.B_matrix[frame][segment] = np.array(
                    [bending_modulus, 0, 0, bending_modulus])

    def get_connection_info(self, script, blob_no, face_no=None, element_no=None, nodes=None, use_first=True):
        """
        Used to acquire the information necessary to set up a rod-blob coupling.
        Params:
            script - an instance of an FFEA_script
            blob_no - ID of the blob with the connection (int)
            nodes - a one-dimensional list\array containing the indices (ints) of the three nodes on the face of the element.
        Returns:
            the element id containing those nodes (int)
        """
        def get_face_node_element(surface, face_no):
            return surface.face[face_no].n, surface.face[face_no].elindex

        def get_element_nodes_face(top, element_no):
            return top.element[element_no].n, [top.element[element_no].get_linear_face(0),
                                               top.element[element_no].get_linear_face(
                                                   1),
                                               top.element[element_no].get_linear_face(
                                                   2),
                                               top.element[element_no].get_linear_face(3)]

        surface = script.load_surface(blob_no)
        top = script.load_topology(blob_no)

        if face_no:
            return get_face_node_element(surface, face_no)

        if element_no:
            print("Not implemented")
            return

        if nodes:
            for element in range(len(top.element)):
                elefaces = top.element[element].n
                if all(elem in elefaces for elem in nodes):
                    attachment_element_index = element
            return attachment_element_index

        print("Supply a face_no or an element_no.")
        return

    def get_euler_angles_from_pdb(self, surface_nodes, r0, r1):
        surf_vec_1 = np.array(surface_nodes[0]) - np.array(surface_nodes[1])
        surf_vec_2 = np.array(surface_nodes[0]) - np.array(surface_nodes[2])
        surf_vec_1 = rod_math.normalize(surf_vec_1)
        surf_vec_2 = rod_math.normalize(surf_vec_2)
        surf_norm = np.cross(surf_vec_1, surf_vec_2)
        p_a_norm = rod_math.normalize(np.array(r1) - np.array(r0))

        direction_dotprod = np.dot(surf_norm, p_a_norm)
        if direction_dotprod < 0:
            surf_norm *= -1

        a = surf_norm
        b = p_a_norm

        v = np.cross(a, b)
        c = np.dot(a, b)
        vx = np.matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        c_factor = 1. / (1 + c)
        R = np.identity(3) + vx + np.linalg.matrix_power(vx, 2) * c_factor

        rm = np.squeeze(np.asarray(R)).flatten()

        euler_angles = rod_math.get_euler_angles(rm)
        return euler_angles


rod_creator = rod_creator()


def load_grace(filename, comment="#"):
    # load header
    header_columns = []
    with open(filename) as grace_file:
        for line in grace_file:
            if comment not in line:
                break
            else:
                line = line.replace("\n", "")
                line = line.replace("#", "")
                line = line.replace("@  ", "")
                line = line.replace("@", "")
                line = line.split(" ")
                header_columns.append(line)

    data = np.genfromtxt(filename, comments="#",
                         skip_header=len(header_columns))

    return header_columns, data
