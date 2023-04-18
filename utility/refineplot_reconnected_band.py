#!/usr/bin/env python
# Author: BaoTing@cmt.tsinghua for RcBand, ScBand, plot_broken_y and compare band
# Author: LiYang@cmt.tsinghua for BandData
# Date: 2021.04.08
# Descripution: This python code is designed for calculate the fermi energy
#                and plot the band structure use the *.Band file.
#              ::Input File::
#                - *.Band
#              ::Output File::
#                - *.Band.png
#                - *.Band.conv
#
# Usage: python3 plot_openmx_band.py -h

import sys
import argparse
import math
import numpy as np
import re
import matplotlib.pyplot as plt
import json
import copy


BOHR = 0.529177210903  # A
HARTREE = 27.21138602  # eV
GREEK_SYMBOLS = ['Gamma', 'Delta', 'Theta', 'Lambda', 'Xi',
                 'Pi', 'Sigma', 'Phi', 'Psi', 'Omega']


class BandData():
    def __init__(self, data_type):
        self.data_type = data_type
        self.band_data = {}
        return

    def file_read(self, filename):
        with open(filename) as frp:
            self.file_lines = frp.readlines()
        return

    def __read_openmx_basic_info(self):
        '''Read int the basic information of the band'''
        # Read in the basic info.
        line = self.file_lines[0].split()
        band_num_each_spin = int(line[0])
        spin_num = int(line[1]) + 1
        band_num = band_num_each_spin * spin_num
        fermi_energy = float(line[2]) * HARTREE

        # Record the info.
        self.band_data["band_num_each_spin"] = band_num_each_spin
        self.band_data["spin_num"] = spin_num
        self.band_data["band_num"] = band_num
        self.band_data["fermi_energy"] = fermi_energy

        # Variours check
        if spin_num not in [1, 2]:
            print("[error] Spin Number ERROR!")
            sys.exit(1)
        return

    def __read_openmx_rlv(self):
        '''Read in the openmx reciprocal lattice vector (rlv)'''
        # Reciprocal lattice vector :: rlv
        # Each row of the 'rlv' matrix is a basis vector
        line = self.file_lines[1].split()
        rlv = np.zeros((3, 3))
        rlv[0, 0] = float(line[0])
        rlv[0, 1] = float(line[1])
        rlv[0, 2] = float(line[2])
        rlv[1, 0] = float(line[3])
        rlv[1, 1] = float(line[4])
        rlv[1, 2] = float(line[5])
        rlv[2, 0] = float(line[6])
        rlv[2, 1] = float(line[7])
        rlv[2, 2] = float(line[8])
        self.band_data["rlv"] = rlv

    def __read_openmx_kpath(self):
        '''Read the openmx band kpath info'''
        # Read in the kpath number
        line = self.file_lines[2]
        kpath_num = int(line)
        self.band_data["kpath_num"] = kpath_num
        # Read in the kpath block
        kpath_densities = [1 for _ in range(kpath_num)]
        kpath_vectors = [[[], []] for _ in range(kpath_num)]
        kpath_symbols = [[] for _ in range(kpath_num)]
        line_i = 2
        for kpath_i in range(kpath_num):
            line_i += 1
            line = self.file_lines[line_i].split()
            kpath_densities[kpath_i] = int(line[0])
            kpath_vectors[kpath_i][0] = [float(line[1]), float(line[2]),
                                         float(line[3])]
            kpath_vectors[kpath_i][1] = [float(line[4]), float(line[5]),
                                         float(line[6])]
            kpath_symbols[kpath_i] = [line[7], line[8]]
        kpoints_num = sum(kpath_densities)
        # Record the data
        self.band_data["kpath_densities"] = kpath_densities
        self.band_data["kpath_vectors"] = kpath_vectors
        self.band_data["kpath_symbols"] = kpath_symbols
        self.band_data["kpoints_num"] = kpoints_num
        return

    def __sort_energys(self, energys):
        '''sort the energys'''
        kpoints_num = len(energys[0])
        band_num = len(energys)
        sorted_energys = np.zeros((band_num, kpoints_num))
        for kpoints_i in range(kpoints_num):
            current_k_column = energys[:, kpoints_i]
            sorted_current_k_column = np.sort(current_k_column)
            sorted_energys[:, kpoints_i] = sorted_current_k_column
        return sorted_energys

    def __read_openmx_energys(self):
        '''Read in the OPENMX energys'''
        band_num_each_spin = self.band_data["band_num_each_spin"]
        kpoints_num = self.band_data["kpoints_num"]
        spin_num = self.band_data["spin_num"]
        # Prepare the data array
        spin_up_energys = np.zeros((band_num_each_spin, kpoints_num))
        if spin_num == 2:
            spin_dn_energys = np.zeros((band_num_each_spin, kpoints_num))
        else:
            spin_dn_energys = []
        kpoint_vectors = [[] for _ in range(kpoints_num)]
        # Read in the data
        line_i = 2 + self.band_data["kpath_num"]
        for kpoint_i in range(kpoints_num):
            # The kpoints line
            line_i += 1
            line = self.file_lines[line_i].split()
            kpoint_vectors[kpoint_i] = [float(line[1]),
                                        float(line[2]),
                                        float(line[3])]
            # The (spin-up) energys line
            line_i += 1
            line = self.file_lines[line_i].split()
            for band_i in range(band_num_each_spin):
                spin_up_energys[band_i, kpoint_i] = float(line[band_i]) * HARTREE
            # The spin-down energys line
            if spin_num == 2:
                line_i += 2
                line = self.file_lines[line_i].split()
                for band_i in range(band_num_each_spin):
                    spin_dn_energys[band_i, kpoint_i] = float(line[band_i]) * HARTREE
        # Post Process the Band energys
        sorted_energys = spin_up_energys
        if spin_num == 2:
            sorted_energys = np.concatenate((spin_up_energys, spin_up_energys),
                                            axis=0)
        # Sort the band
        sorted_energys = self.__sort_energys(sorted_energys)
        # Record the data
        self.band_data["kpoint_vectors"] = kpoint_vectors
        self.band_data["spin_up_energys"] = spin_up_energys
        self.band_data["spin_dn_energys"] = spin_dn_energys
        self.band_data["sorted_energys"] = sorted_energys
        return

    def __cal_k_distance(self, rlv, beg_kpt_frac, end_kpt_frac, distance_unit=1):
        '''Calcualte the k-cooridnate in Cartesian indicator'''
        beg_kpt_cart = np.array([0.0, 0.0, 0.0])
        end_kpt_cart = np.array([0.0, 0.0, 0.0])
        #
        # k_cart = k_frac * rlv
        #                               __              __
        #                               | b_1x b_1y b_1z |
        #        = (kf_1, kf_2, kf_3) * | b_2x b_2y b_2z |
        #                               | b_3x b_3y b_3z |
        #                               --              --
        #        = (kc_x, kc_y, kc_z)
        #
        for xyz in range(3):
            beg_kpt_cart[xyz] = 0.0
            for b_i in range(3):
                beg_kpt_cart[xyz] += beg_kpt_frac[b_i] * rlv[b_i, xyz]
                end_kpt_cart[xyz] += end_kpt_frac[b_i] * rlv[b_i, xyz]
        # Calculate the k distance of the two kpoints
        k_distance = math.sqrt(sum((beg_kpt_cart-end_kpt_cart) ** 2))
        k_distance /= distance_unit
        return k_distance

    def __get_kpt_coords(self, distance_unit):
        '''Get the coords of each kpoints in k-space'''
        kpath_num = self.band_data["kpath_num"]
        kpath_vectors = self.band_data["kpath_vectors"]
        kpoint_vectors = self.band_data["kpoint_vectors"]
        kpath_densities = self.band_data["kpath_densities"]
        rlv = self.band_data["rlv"]
        kpoints_num = self.band_data["kpoints_num"]
        # Prepare the data list
        hsk_distance_list = np.zeros(kpath_num)
        sum_hsk_distance_list = np.zeros(kpath_num)
        kpoints_coords = np.zeros(kpoints_num)
        hsk_coords = np.zeros(kpath_num+1)
        # Get the distance for high-symmetry kpoints
        for kpath_i in range(kpath_num):
            start_hsk = kpath_vectors[kpath_i][0]
            end_hsk = kpath_vectors[kpath_i][1]
            hsk_distance_list[kpath_i] = \
                self.__cal_k_distance(rlv, start_hsk, end_hsk, distance_unit)
            sum_hsk_distance_list[kpath_i] = \
                sum(hsk_distance_list[0:kpath_i+1])
        hsk_coords[1:] = sum_hsk_distance_list
        # Get the distance in k-space of k-points on the k-path
        kpoints_i = -1
        for kpath_i in range(kpath_num):
            # Count the Previous kpath distance
            pre_path_distance = hsk_coords[kpath_i]
            # Calculate the kpoints' distance in current kpath
            for _ in range(kpath_densities[kpath_i]):
                kpoints_i += 1
                start_hsk = kpath_vectors[kpath_i][0]
                end_hsk = kpoint_vectors[kpoints_i]
                # The total distance equals to (pre_dis + curr_dis)
                kpoints_coords[kpoints_i] = pre_path_distance + \
                    self.__cal_k_distance(rlv, start_hsk, end_hsk, distance_unit)
        # Record the data
        self.band_data["hsk_coords"] = hsk_coords
        self.band_data["kpoints_coords"] = kpoints_coords
        return

    def __refine_fermi_energy(self):
        '''Refine the fermi energy and the center of HOMO and LUMO'''
        fermi_energy = self.band_data["fermi_energy"]
        energys = self.band_data["sorted_energys"]
        band_num = self.band_data["band_num"]
        kpoints_num = self.band_data["kpoints_num"]
        spin_num = self.band_data["spin_num"]
        # find the LUMO and HOMO
        min_homo_ediff = fermi_energy - energys[0, 0]
        homo_band_index = 0
        homo_kpt_index = 0
        min_lumo_ediff = energys[band_num-1, 0] - fermi_energy
        lumo_band_index = band_num-1
        lumo_kpt_index = 0
        for band_i in range(band_num):
            for kpoint_i in range(kpoints_num):
                curr_energy = energys[band_i, kpoint_i]
                lumo_ediff = curr_energy - fermi_energy
                homo_ediff = fermi_energy - curr_energy
                if (lumo_ediff >= 0) and (lumo_ediff < min_lumo_ediff):
                    min_lumo_ediff = lumo_ediff
                    lumo_band_index = band_i
                    lumo_kpt_index = kpoint_i
                elif (homo_ediff > 0) and (homo_ediff < min_homo_ediff):
                    min_homo_ediff = homo_ediff
                    homo_band_index = band_i
                    homo_kpt_index = kpoint_i
        lumo_energy = energys[lumo_band_index, lumo_kpt_index]
        homo_energy = energys[homo_band_index, homo_kpt_index]
        refined_fermi_energy = (lumo_energy + homo_energy) / 2

        # baot: when using read from json, for gapped system, use refined energy, ignore lihe's comment
        # lihe: no need to do refine
        # refined_fermi_energy = fermi_energy

        # Shift the zero energy to the fermi level
        self.band_data["origin_spin_up_energys"] = self.band_data["spin_up_energys"]
        self.band_data["origin_spin_dn_energys"] = self.band_data["spin_dn_energys"]
        self.band_data["origin_sorted_energys"] = self.band_data["sorted_energys"]
        self.band_data["spin_up_energys"] -= refined_fermi_energy
        self.band_data["sorted_energys"] -= refined_fermi_energy
        if spin_num == 2:
            self.band_data["spin_up_energys"] -= refined_fermi_energy  # might be spin_dn??? -baot
        # Record the data
        self.band_data["refined_fermi_energy"] = refined_fermi_energy
        self.band_data["lumo_energy"] = lumo_energy
        self.band_data["homo_energy"] = homo_energy
        self.band_data["lumo_band_index"] = lumo_band_index
        self.band_data["lumo_kpt_index"] = lumo_kpt_index
        self.band_data["homo_band_index"] = homo_band_index
        self.band_data["homo_kpt_index"] = homo_kpt_index
        self.band_data["min_lumo_ediff"] = min_lumo_ediff
        self.band_data["min_homo_ediff"] = min_homo_ediff
        return

    def __prepare_plot_kpt_symbol(self):
        '''Prepare the kpoints symbols for the plot'''
        kpath_num = self.band_data["kpath_num"]
        kpath_symbols = self.band_data["kpath_symbols"]
        # Prepare the symbol of k-axis (xtics)
        hsk_symbols = ['' for _ in range(kpath_num+1)]
        # Set
        hsk_symbols[0] = kpath_symbols[0][0]
        for kpath_i in range(1, kpath_num):
            if kpath_symbols[kpath_i][0] == kpath_symbols[kpath_i-1][1]:
                hsk_symbols[kpath_i] = kpath_symbols[kpath_i][0]
            else:
                hsk_symbols[kpath_i] = "%s|%s" % (kpath_symbols[kpath_i - 1][1],
                                                  kpath_symbols[kpath_i][0])
        hsk_symbols[-1] = kpath_symbols[-1][1]
        # Plot the Band
        plot_hsk_symbols = []
        for symbol in hsk_symbols:
            symbol = symbol.replace("\\", "")
            for greek_symbol in GREEK_SYMBOLS:
                if greek_symbol == 'Gamma':
                    latex_greek_symbol = 'Γ'
                else:
                    latex_greek_symbol = "$\\" + greek_symbol + "$"
                symbol = re.sub(greek_symbol, "orz", symbol,
                                flags=re.I)
                symbol = symbol.replace("orz", latex_greek_symbol)
            symbol = re.sub(r'_\d+', lambda x: '$'+x[0]+'$', symbol)
            plot_hsk_symbols.append(symbol)
        # Record the data
        self.band_data["hsk_symbols"] = hsk_symbols
        self.band_data["plot_hsk_symbols"] = plot_hsk_symbols
        return

    def get_band_data(self):
        '''Get the band data'''
        # Read in data
        if self.data_type == 'openmx':
            distance_unit = 2 * math.pi * BOHR
            self.__read_openmx_basic_info()
            self.__read_openmx_rlv()
            self.__read_openmx_kpath()
            self.__read_openmx_energys()
        elif self.data_type == 'vasp':
            print("[TODO]")
            sys.exit(0)
        else:
            print("[TODO]")
            sys.exit(0)
        # Post process
        self.__get_kpt_coords(distance_unit)
        self.__refine_fermi_energy()
        self.__prepare_plot_kpt_symbol()
        return


class RcBandData():
    '''
    a class for reconnected band data processing
    '''

    def __init__(self, data_type):
        self.data_type = data_type
        self.band_data = {}
        self.dft_data = {}
        return

    def refined_reconnected_fermi_energy(self):
        '''by baot, tested on black P(P4) only
        used for band_reconnected.json from connect_interpolate.py
        The band data here is about the fermi_energy=0.0 already, 
        here we move the fermi level to  the middle of the gap
        '''

        band_num = self.band_data["band_num"]
        kpoints_num = self.band_data["kpoints_num"]
        spin_num = self.band_data["spin_num"]
        lumo_energy = self.band_data["lumo_energy"]
        homo_energy = self.band_data["homo_energy"]
        self.band_data['refined_fermi_energy'] = (lumo_energy + homo_energy) / 2

        # baot: the energy shift required to center the fermi energy for reconnected band data
        energy_shift = self.band_data["fermi_energy"] - self.band_data['refined_fermi_energy']
        self.band_data["spin_up_energys"] += energy_shift
        if spin_num == 2:
            self.band_data["spin_dn_energys"] += energy_shift

        # find the LUMO and HOMO
        energys = self.band_data["spin_up_energys"]  # band_data["sorted_energys"]
        homo_band_index = 0
        homo_kpt_index = 0
        lumo_band_index = band_num-1
        lumo_kpt_index = 0
        lumo_energy = 1e5
        homo_energy = -1e5
        for band_i in range(band_num):
            for kpoint_i in range(kpoints_num):
                curr_energy = energys[band_i, kpoint_i]
                if (curr_energy >= 0) and (curr_energy < lumo_energy):
                    lumo_energy = curr_energy
                    lumo_band_index = band_i
                    lumo_kpt_index = kpoint_i
                if (curr_energy < 0) and (curr_energy > homo_energy):
                    homo_energy = curr_energy
                    homo_band_index = band_i
                    homo_kpt_index = kpoint_i

        # Record the data
        self.band_data["lumo_energy"] = lumo_energy
        self.band_data["homo_energy"] = homo_energy
        self.band_data["lumo_band_index"] = lumo_band_index
        self.band_data["lumo_kpt_index"] = lumo_kpt_index
        self.band_data["homo_band_index"] = homo_band_index
        self.band_data["homo_kpt_index"] = homo_kpt_index
        return

    def get_reconnected_bandwidth(self, band_index=0, gap_tol=0.005, flat_tol=0.006):
        '''by baot,
        used for band_reconnected.json from connect_interpolate.py
        count all the bandwidth, work for spin_num = 1 only,
        for system with clean dirac cone like graphene, the gap tol is important to identify whether gapped
        '''
        band_num = self.band_data["band_num"]
        kpoints_num = self.band_data["kpoints_num"]
        energys = self.band_data["spin_up_energys"]  # band_data["sorted_energys"]
        gapped = True
        bandwidth_info = {}
        cross_fermi_index = []
        flat_index = []
        for band_i in range(band_num):
            band_energy = energys[band_i, :]
            temp_dict = {}
            temp_dict['max'] = max(band_energy)
            temp_dict['min'] = min(band_energy)
            temp_dict['width'] = temp_dict['max']-temp_dict['min']
            bandwidth_info[band_i] = temp_dict
            if temp_dict['width'] < flat_tol:
                flat_index.append(band_i)
            if temp_dict['max'] > 0 and temp_dict['min'] < 0:
                gapped = False
                cross_fermi_index.append(band_i)
                continue
            '''
      # TODO
      elif temp_dict['max']<0 and temp_dict['max']> 0-gap_tol:
        # check 3 neighbour bands to find dirac cone type system
      '''
        self.band_data['gapped'] = gapped
        self.band_data['gap_size'] = self.band_data['lumo_energy'] - self.band_data['homo_energy']
        self.band_data['cross_fermi_index'] = cross_fermi_index
        self.band_data['flat_index'] = flat_index
        self.band_data['bandwidth_info'] = bandwidth_info
        return

    def advise_plot_range(self, plot_args, cover_band=10):
        '''
        advise_plot_range will change plot_args["min/max_plot_energy"] to cover at least 10 bands above/under E_F
        '''
        energys = self.band_data["spin_up_energys"]
        cross_fermi_index = self.band_data['cross_fermi_index']
        homo_index = self.band_data['homo_band_index']
        bandwidth_info = self.band_data['bandwidth_info']

        if len(cross_fermi_index) == 0:
            index_up = homo_index + cover_band
            index_dn = homo_index - cover_band + 1
        else:
            index_up = max(cross_fermi_index) + cover_band
            index_dn = min(cross_fermi_index) - cover_band + 1
        if index_dn < 0:
            index_dn = 0
        up = bandwidth_info[index_up]['max']
        dn = bandwidth_info[index_dn]['min']
        up = (up//0.5)*0.5 + 0.5
        dn = (dn//0.5)*0.5
        boundary = max(up, 0-dn)
        plot_args["min_plot_energy"] = 0-boundary
        plot_args["max_plot_energy"] = boundary
        self.band_data['advise_plot_range'] = [-boundary, boundary]
        return plot_args


class ScBandData():
    '''scatter band data class
    '''

    def __init__(self, data_type):
        self.data_type = data_type
        self.band_data = {}
        return

    def file_read(self, filename):
        with open(filename, encoding='utf-8') as frp:
            self.file_lines = frp.readlines()
        return

    def __read_openmx_basic_info(self):
        '''Read int the basic information of the band'''
        # Read in the basic info.
        line = self.file_lines[0].split()
        band_num_each_spin = int(line[0])
        spin_num = int(line[1]) + 1
        band_num = band_num_each_spin * spin_num
        fermi_energy = float(line[2]) * HARTREE

        # Record the info.
        self.band_data["band_num_each_spin"] = band_num_each_spin
        self.band_data["spin_num"] = spin_num
        self.band_data["band_num"] = band_num
        self.band_data["fermi_energy"] = fermi_energy

        # Variours check
        if spin_num not in [1, 2]:
            print("[error] Spin Number ERROR!")
            sys.exit(1)
        return

    def __read_openmx_rlv(self):
        '''Read in the openmx reciprocal lattice vector (rlv)'''
        # Reciprocal lattice vector :: rlv
        # Each row of the 'rlv' matrix is a basis vector
        line = self.file_lines[1].split()
        rlv = np.zeros((3, 3))
        rlv[0, 0] = float(line[0])
        rlv[0, 1] = float(line[1])
        rlv[0, 2] = float(line[2])
        rlv[1, 0] = float(line[3])
        rlv[1, 1] = float(line[4])
        rlv[1, 2] = float(line[5])
        rlv[2, 0] = float(line[6])
        rlv[2, 1] = float(line[7])
        rlv[2, 2] = float(line[8])
        self.band_data["rlv"] = rlv

    def __read_openmx_kpath(self):
        '''Read the openmx band kpath info'''
        # Read in the kpath number
        line = self.file_lines[2]
        kpath_num = int(line)
        self.band_data["kpath_num"] = kpath_num
        # Read in the kpath block
        kpath_densities = [1 for _ in range(kpath_num)]
        kpath_vectors = [[[], []] for _ in range(kpath_num)]
        kpath_symbols = [[] for _ in range(kpath_num)]
        line_i = 2
        for kpath_i in range(kpath_num):
            line_i += 1
            line = self.file_lines[line_i].split()
            kpath_densities[kpath_i] = int(line[0])
            kpath_vectors[kpath_i][0] = [float(line[1]), float(line[2]),
                                         float(line[3])]
            kpath_vectors[kpath_i][1] = [float(line[4]), float(line[5]),
                                         float(line[6])]
            kpath_symbols[kpath_i] = [line[7], line[8]]
        kpoints_num = sum(kpath_densities)
        # Record the data
        self.band_data["kpath_densities"] = kpath_densities
        self.band_data["kpath_vectors"] = kpath_vectors
        self.band_data["kpath_symbols"] = kpath_symbols
        self.band_data["kpoints_num"] = kpoints_num
        return

    def __sort_energys(self, energys):
        '''sort the energys'''
        kpoints_num = len(energys[0])
        band_num = len(energys)
        sorted_energys = np.zeros((band_num, kpoints_num))
        for kpoints_i in range(kpoints_num):
            current_k_column = energys[:, kpoints_i]
            sorted_current_k_column = np.sort(current_k_column)
            sorted_energys[:, kpoints_i] = sorted_current_k_column
        return sorted_energys

    def __read_openmx_energys(self):
        '''Read in the OPENMX energys'''
        band_num_each_spin = self.band_data["band_num_each_spin"]
        kpoints_num = self.band_data["kpoints_num"]
        spin_num = self.band_data["spin_num"]
        # Prepare the data array
        spin_up_energys = np.zeros((band_num_each_spin, kpoints_num))
        if spin_num == 2:
            spin_dn_energys = np.zeros((band_num_each_spin, kpoints_num))
        else:
            spin_dn_energys = []
        kpoint_vectors = [[] for _ in range(kpoints_num)]
        # Read in the data
        line_i = 2 + self.band_data["kpath_num"]
        for kpoint_i in range(kpoints_num):
            # The kpoints line
            line_i += 1
            line = self.file_lines[line_i].split()
            kpoint_vectors[kpoint_i] = [float(line[1]),
                                        float(line[2]),
                                        float(line[3])]
            # The (spin-up) energys line
            line_i += 1
            line = self.file_lines[line_i].split()
            for band_i in range(band_num_each_spin):
                spin_up_energys[band_i, kpoint_i] = float(line[band_i]) * HARTREE
            # The spin-down energys line
            if spin_num == 2:
                line_i += 2
                line = self.file_lines[line_i].split()
                for band_i in range(band_num_each_spin):
                    spin_dn_energys[band_i, kpoint_i] = float(line[band_i]) * HARTREE
        # Post Process the Band energys
        sorted_energys = spin_up_energys
        if spin_num == 2:
            sorted_energys = np.concatenate((spin_up_energys, spin_up_energys),
                                            axis=0)

        # lihe: Read band from each calc
        assert spin_num == 1
        print(spin_up_energys.shape)
        for root, folders, files in os.walk(self.egval_path):
            if 'egval.dat' in files:
                index_k = int(os.path.split(root)[-1]) - 1
                spin_up_energys[:, index_k] = np.loadtxt(os.path.join(root, 'egval.dat'))

        # Sort the band
        sorted_energys = self.__sort_energys(sorted_energys)
        # Record the data
        self.band_data["kpoint_vectors"] = kpoint_vectors
        self.band_data["spin_up_energys"] = spin_up_energys
        self.band_data["spin_dn_energys"] = spin_dn_energys
        self.band_data["sorted_energys"] = sorted_energys
        return

    def __cal_k_distance(self, rlv, beg_kpt_frac, end_kpt_frac, distance_unit=1):
        '''Calcualte the k-cooridnate in Cartesian indicator'''
        beg_kpt_cart = np.array([0.0, 0.0, 0.0])
        end_kpt_cart = np.array([0.0, 0.0, 0.0])
        #
        # k_cart = k_frac * rlv
        #                               __              __
        #                               | b_1x b_1y b_1z |
        #        = (kf_1, kf_2, kf_3) * | b_2x b_2y b_2z |
        #                               | b_3x b_3y b_3z |
        #                               --              --
        #        = (kc_x, kc_y, kc_z)
        #
        for xyz in range(3):
            beg_kpt_cart[xyz] = 0.0
            for b_i in range(3):
                beg_kpt_cart[xyz] += beg_kpt_frac[b_i] * rlv[b_i, xyz]
                end_kpt_cart[xyz] += end_kpt_frac[b_i] * rlv[b_i, xyz]
        # Calculate the k distance of the two kpoints
        k_distance = math.sqrt(sum((beg_kpt_cart-end_kpt_cart) ** 2))
        k_distance /= distance_unit
        return k_distance

    def __get_kpt_coords(self, distance_unit):
        '''Get the coords of each kpoints in k-space'''
        kpath_num = self.band_data["kpath_num"]
        kpath_vectors = self.band_data["kpath_vectors"]
        kpoint_vectors = self.band_data["kpoint_vectors"]
        kpath_densities = self.band_data["kpath_densities"]
        rlv = self.band_data["rlv"]
        kpoints_num = self.band_data["kpoints_num"]
        # Prepare the data list
        hsk_distance_list = np.zeros(kpath_num)
        sum_hsk_distance_list = np.zeros(kpath_num)
        kpoints_coords = np.zeros(kpoints_num)
        hsk_coords = np.zeros(kpath_num+1)
        # Get the distance for high-symmetry kpoints
        for kpath_i in range(kpath_num):
            start_hsk = kpath_vectors[kpath_i][0]
            end_hsk = kpath_vectors[kpath_i][1]
            hsk_distance_list[kpath_i] = \
                self.__cal_k_distance(rlv, start_hsk, end_hsk, distance_unit)
            sum_hsk_distance_list[kpath_i] = \
                sum(hsk_distance_list[0:kpath_i+1])
        hsk_coords[1:] = sum_hsk_distance_list
        # Get the distance in k-space of k-points on the k-path
        kpoints_i = -1
        for kpath_i in range(kpath_num):
            # Count the Previous kpath distance
            pre_path_distance = hsk_coords[kpath_i]
            # Calculate the kpoints' distance in current kpath
            for _ in range(kpath_densities[kpath_i]):
                kpoints_i += 1
                start_hsk = kpath_vectors[kpath_i][0]
                end_hsk = kpoint_vectors[kpoints_i]
                # The total distance equals to (pre_dis + curr_dis)
                kpoints_coords[kpoints_i] = pre_path_distance + \
                    self.__cal_k_distance(rlv, start_hsk, end_hsk, distance_unit)
        # Record the data
        self.band_data["hsk_coords"] = hsk_coords
        self.band_data["kpoints_coords"] = kpoints_coords
        return

    def __refine_fermi_energy(self):
        '''Refine the fermi energy and the center of HOMO and LUMO'''
        fermi_energy = self.band_data["fermi_energy"]
        energys = self.band_data["sorted_energys"]
        band_num = self.band_data["band_num"]
        kpoints_num = self.band_data["kpoints_num"]
        spin_num = self.band_data["spin_num"]
        # find the LUMO and HOMO
        min_homo_ediff = fermi_energy - energys[0, 0]
        homo_band_index = 0
        homo_kpt_index = 0
        min_lumo_ediff = energys[band_num-1, 0] - fermi_energy
        lumo_band_index = band_num-1
        lumo_kpt_index = 0
        for band_i in range(band_num):
            for kpoint_i in range(kpoints_num):
                curr_energy = energys[band_i, kpoint_i]
                lumo_ediff = curr_energy - fermi_energy
                homo_ediff = fermi_energy - curr_energy
                if (lumo_ediff >= 0) and (lumo_ediff < min_lumo_ediff):
                    min_lumo_ediff = lumo_ediff
                    lumo_band_index = band_i
                    lumo_kpt_index = kpoint_i
                elif (homo_ediff > 0) and (homo_ediff < min_homo_ediff):
                    min_homo_ediff = homo_ediff
                    homo_band_index = band_i
                    homo_kpt_index = kpoint_i
        lumo_energy = energys[lumo_band_index, lumo_kpt_index]
        homo_energy = energys[homo_band_index, homo_kpt_index]
        refined_fermi_energy = (lumo_energy + homo_energy) / 2

        # lihe: no need to do refine
        refined_fermi_energy = fermi_energy

        # Shift the zero energy to the fermi level
        self.band_data["origin_spin_up_energys"] = self.band_data["spin_up_energys"]
        self.band_data["origin_spin_dn_energys"] = self.band_data["spin_dn_energys"]
        self.band_data["origin_sorted_energys"] = self.band_data["sorted_energys"]
        self.band_data["spin_up_energys"] -= refined_fermi_energy
        self.band_data["sorted_energys"] -= refined_fermi_energy
        if spin_num == 2:
            self.band_data["spin_up_energys"] -= refined_fermi_energy
        # Record the data
        self.band_data["refined_fermi_energy"] = refined_fermi_energy
        self.band_data["lumo_energy"] = lumo_energy
        self.band_data["homo_energy"] = homo_energy
        self.band_data["lumo_band_index"] = lumo_band_index
        self.band_data["lumo_kpt_index"] = lumo_kpt_index
        self.band_data["homo_band_index"] = homo_band_index
        self.band_data["homo_kpt_index"] = homo_kpt_index
        self.band_data["min_lumo_ediff"] = min_lumo_ediff
        self.band_data["min_homo_ediff"] = min_homo_ediff
        return

    def __prepare_plot_kpt_symbol(self):
        '''Prepare the kpoints symbols for the plot'''
        kpath_num = self.band_data["kpath_num"]
        kpath_symbols = self.band_data["kpath_symbols"]
        # Prepare the symbol of k-axis (xtics)
        hsk_symbols = ['' for _ in range(kpath_num+1)]
        # Set
        hsk_symbols[0] = kpath_symbols[0][0]
        for kpath_i in range(1, kpath_num):
            if kpath_symbols[kpath_i][0] == kpath_symbols[kpath_i-1][1]:
                hsk_symbols[kpath_i] = kpath_symbols[kpath_i][0]
            else:
                hsk_symbols[kpath_i] = "%s|%s" % (kpath_symbols[kpath_i - 1][1],
                                                  kpath_symbols[kpath_i][0])
        hsk_symbols[-1] = kpath_symbols[-1][1]
        # Plot the Band
        plot_hsk_symbols = []
        for symbol in hsk_symbols:
            symbol = symbol.replace("\\", "")
            for greek_symbol in GREEK_SYMBOLS:
                if greek_symbol == 'Gamma':
                    latex_greek_symbol = 'Γ'
                else:
                    latex_greek_symbol = "$\\" + greek_symbol + "$"
                symbol = re.sub(greek_symbol, "orz", symbol,
                                flags=re.I)
                symbol = symbol.replace("orz", latex_greek_symbol)
            symbol = re.sub(r'_\d+', lambda x: '$'+x[0]+'$', symbol)
            plot_hsk_symbols.append(symbol)
        # Record the data
        self.band_data["hsk_symbols"] = hsk_symbols
        self.band_data["plot_hsk_symbols"] = plot_hsk_symbols
        return

    def get_band_data(self, read_path='./egval_k'):
        '''Get the band data'''
        self.egval_path = read_path
        # Read in data
        if self.data_type == 'openmx':
            distance_unit = 2 * math.pi * BOHR
            self.__read_openmx_basic_info()
            self.__read_openmx_rlv()
            self.__read_openmx_kpath()
            self.__read_openmx_energys()
        elif self.data_type == 'vasp':
            print("[TODO]")
            sys.exit(0)
        else:
            print("[TODO]")
            sys.exit(0)
        # Post process
        self.__get_kpt_coords(distance_unit)
        self.__refine_fermi_energy()
        self.__prepare_plot_kpt_symbol()
        return


def get_command_line_input():
    '''Read in the command line parameters'''
    parser = argparse.ArgumentParser("Basic band plot parameters")
    parser.add_argument('-t', '--type', dest='data_type',
                        default='openmx', type=str, choices=['openmx', 'vasp'],
                        help='Type of the band calculation.')
    parser.add_argument('-d', '--ymin', dest='min_plot_energy',
                        default=-3, type=float,
                        help='Minimal plot energy windows.')
    parser.add_argument('-u', '--ymax', dest='max_plot_energy',
                        default=3, type=float,
                        help='Maximal plot energy windows.')
    parser.add_argument('-f', '--format', dest='plot_format',
                        default='png', type=str, choices=['png', 'eps', 'pdf'],
                        help='Plot format.')
    parser.add_argument('-i', '--dpi', dest='plot_dpi',
                        default=400, type=int,
                        help='Plot resolution (dpi).')
    parser.add_argument('-a', '--datafile', dest='data_filename',
                        default='openmx.Band', type=str,
                        help='Input data filename.')
    parser.add_argument('-o', '--output', dest='file_tag',
                        default='band', type=str,
                        help='Output file tag name.')
    parser.add_argument('-x', '--no-plot', dest='no_plot', action='store_const',
                        const=True, default=False,
                        help='Do not plot the band.')
    args = parser.parse_args()
    plot_args = {"data_type": args.data_type,
                 "min_plot_energy": args.min_plot_energy,
                 "max_plot_energy": args.max_plot_energy,
                 "plot_format": args.plot_format,
                 "plot_dpi": args.plot_dpi,
                 "file_tag": args.file_tag,
                 "data_filename": args.data_filename,
                 "no_plot": args.no_plot}
    return plot_args


def band_save_to_json(data, file_tag):
    '''Save the band data to json file'''
    filename = "%s.json" % file_tag
    data_save = copy.deepcopy(data)
    for key, vals in data_save.items():
        if type(vals) is np.ndarray:
            data_save[key] = vals.tolist()
    with open(filename, 'w') as jfwp:
        json.dump(data_save, jfwp)
    return


def dft_band_plot(band_data_obj, plot_args):
    '''band plot function'''
    hsk_coords = band_data_obj.band_data["hsk_coords"]
    plot_hsk_symbols = band_data_obj.band_data["plot_hsk_symbols"]
    kpath_num = band_data_obj.band_data["kpath_num"]
    band_num_each_spin = band_data_obj.band_data["band_num_each_spin"]
    kpoints_coords = band_data_obj.band_data["kpoints_coords"]
    spin_num = band_data_obj.band_data["spin_num"]
    spin_up_energys = band_data_obj.band_data["spin_up_energys"]
    spin_dn_energys = band_data_obj.band_data["spin_dn_energys"]
    min_plot_energy = plot_args["min_plot_energy"]
    max_plot_energy = plot_args["max_plot_energy"]
    file_tag = plot_args["file_tag"]
    plot_format = plot_args["plot_format"]
    plot_dpi = plot_args["plot_dpi"]
    # Design the Figure
    # For GUI less server
    plt.switch_backend('agg')
    # Set the Fonts
    # plt.rcParams.update({'font.size': 14,
    #                      'font.family': 'STIXGeneral',
    #                      'mathtext.fontset': 'stix'})
    plt.rcParams.update({'font.size': 22,
                         'font.family': 'Arial',
                         'mathtext.fontset': 'cm'})
    # Set the spacing between the axis and labels
    plt.rcParams['xtick.major.pad'] = '6'
    plt.rcParams['ytick.major.pad'] = '6'
    # Set the ticks 'inside' the axis
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    # Create the figure and axis object
    fig = plt.figure(figsize=(5.5, 5.5))
    band_plot = fig.add_subplot(1, 1, 1)
    # Set the range of plot
    x_min = 0.0
    x_max = hsk_coords[-1]
    y_min = min_plot_energy
    y_max = max_plot_energy
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # Set the label of x and y axis
    plt.xlabel('')
    plt.ylabel('Energy (eV)')
    # Set the Ticks of x and y axis
    plt.xticks(hsk_coords)
    band_plot.set_xticklabels(plot_hsk_symbols)
    plt.yticks(size=14)
    # Plot the solid lines for High symmetic k-points
    for kpath_i in range(kpath_num+1):
        plt.vlines(hsk_coords[kpath_i], y_min, y_max, colors="black", linewidth=0.7)
    # Plot the fermi energy surface with a dashed line
    plt.hlines(0.0, x_min, x_max, colors="black",
               linestyles="dashed", linewidth=0.7)
    # Plot the Band Structure
    for band_i in range(band_num_each_spin):
        x = kpoints_coords
        y = spin_up_energys[band_i]
        band_plot.plot(x, y, 'r-', linewidth=1.5)
    if spin_num == 2:
        for band_i in range(band_num_each_spin):
            x = kpoints_coords
            y = spin_dn_energys[band_i]
            band_plot.plot(x, y, '-', color='#0564c3', linewidth=1)
    # Save the figure
    plot_filename = "%s.%s" % (file_tag, plot_format)
    plt.tight_layout()
    plt.savefig(plot_filename, format=plot_format, dpi=plot_dpi, transparent=True)
    plt.savefig('band.svg', transparent=True)
    return


def band_plot(band_data_obj, plot_args, savepath):
    '''band plot function'''
    hsk_coords = band_data_obj.band_data["hsk_coords"]
    plot_hsk_symbols = band_data_obj.band_data["plot_hsk_symbols"]
    kpath_num = band_data_obj.band_data["kpath_num"]
    band_num_each_spin = band_data_obj.band_data["band_num_each_spin"]
    kpoints_coords = band_data_obj.band_data["kpoints_coords"]
    spin_num = band_data_obj.band_data["spin_num"]
    spin_up_energys = band_data_obj.band_data["spin_up_energys"]
    spin_dn_energys = band_data_obj.band_data["spin_dn_energys"]
    min_plot_energy = plot_args["min_plot_energy"]
    max_plot_energy = plot_args["max_plot_energy"]
    file_tag = plot_args["file_tag"]
    plot_format = plot_args["plot_format"]
    plot_dpi = plot_args["plot_dpi"]
    # Design the Figure
    # For GUI less server
    plt.switch_backend('agg')
    # Set the Fonts
    # plt.rcParams.update({'font.size': 14,
    #                      'font.family': 'STIXGeneral',
    #                      'mathtext.fontset': 'stix'})
    plt.rcParams.update({'font.size': 14,
                         'font.family': 'Arial',
                         'mathtext.fontset': 'cm'})
    # Set the spacing between the axis and labels
    plt.rcParams['xtick.major.pad'] = '6'
    plt.rcParams['ytick.major.pad'] = '6'
    # Set the ticks 'inside' the axis
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    # Create the figure and axis object
    fig = plt.figure(figsize=(5.5, 5.5))
    band_plot = fig.add_subplot(1, 1, 1)
    # Set the range of plot
    x_min = 0.0
    x_max = hsk_coords[-1]
    y_min = min_plot_energy
    y_max = max_plot_energy
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # Set the label of x and y axis
    plt.xlabel('')
    plt.ylabel('Energy (eV)')
    # Set the Ticks of x and y axis
    plt.xticks(hsk_coords)
    band_plot.set_xticklabels(plot_hsk_symbols)
    plt.yticks(size=10)
    # Plot the solid lines for High symmetic k-points
    for kpath_i in range(kpath_num+1):
        plt.vlines(hsk_coords[kpath_i], y_min, y_max, colors="black", linewidth=0.7)
    # Plot the fermi energy surface with a dashed line
    plt.hlines(0.0, x_min, x_max, colors="black",
               linestyles="dashed", linewidth=0.7)

    # Plot the Band Structure
    for band_i in range(band_num_each_spin):
        x = kpoints_coords
        y = spin_up_energys[band_i]
        # band_plot.plot(x, y, 'b-', linewidth=1.5)
        band_plot.scatter(x, y, c='blue', s=0.8)

    if spin_num == 2:
        for band_i in range(band_num_each_spin):
            x = kpoints_coords
            y = spin_dn_energys[band_i]
            band_plot.plot(x, y, '-', color='#0564c3', linewidth=1)  # color='#0564c3'
    # Save the figure
    plot_filename = "%s.%s" % (file_tag, plot_format)
    plt.tight_layout()
    plt.savefig(savepath+'/band_reconnected.png', format=plot_format, dpi=plot_dpi, transparent=True)
    plt.savefig(savepath+'/band_reconnected.svg', transparent=True)
    plt.close()
    return


def band_plot_brokeny(band_data_obj, plot_args, cut_energy, highlight_num=1, savepath='./'):
    '''
    baot: band plot function for P4, 
    hightlight the flat conduction band,
    set the broken y axis
    '''

    # at most 4 flat bands to show
    # https://www.sioe.cn/yingyong/yanse-rgb-16/
    colorlist = ['red', 'royalblue', 'orange', 'forestgreen']
    colorlist_RGB = [(255, 0, 0), (65, 105, 225), (255, 165, 0), (34, 139, 34)]
    homo_index = band_data_obj.band_data['homo_band_index']
    highlight_band = [i for i in range(homo_index, homo_index-highlight_num, -1)]  # band in special color

    hsk_coords = band_data_obj.band_data["hsk_coords"]
    plot_hsk_symbols = band_data_obj.band_data["plot_hsk_symbols"]
    kpath_num = band_data_obj.band_data["kpath_num"]
    band_num_each_spin = band_data_obj.band_data["band_num_each_spin"]
    kpoints_coords = band_data_obj.band_data["kpoints_coords"]
    spin_num = band_data_obj.band_data["spin_num"]
    spin_up_energys = band_data_obj.band_data["spin_up_energys"]
    spin_dn_energys = band_data_obj.band_data["spin_dn_energys"]
    min_plot_energy = plot_args["min_plot_energy"]
    max_plot_energy = plot_args["max_plot_energy"]
    file_tag = plot_args["file_tag"]
    plot_format = plot_args["plot_format"]
    plot_dpi = plot_args["plot_dpi"]
    # Design the Figure
    # For GUI less server
    plt.switch_backend('agg')
    # Set the Fonts
    # plt.rcParams.update({'font.size': 14,
    #                      'font.family': 'STIXGeneral',
    #                      'mathtext.fontset': 'stix'})
    plt.rcParams.update({'font.size': 14,
                         'font.family': 'Arial',
                         'mathtext.fontset': 'cm'})
    # Set the spacing between the axis and labels
    plt.rcParams['xtick.major.pad'] = '4'
    plt.rcParams['ytick.major.pad'] = '4'
    # Set the ticks 'inside' the axis
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # Create the figure and axis object
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 7.5), sharex=True)
    fig.subplots_adjust(hspace=0.01)  # adjust space between axes
    ax1.spines.bottom.set_visible(False)
    ax2.spines.top.set_visible(False)
    ax1.xaxis.tick_top()

    # Set the range of plot
    x_min = 0.0
    x_max = hsk_coords[-1]
    y_min = min_plot_energy
    y_max = max_plot_energy
    # plt.xlim(x_min, x_max)
    # plt.ylim(y_min, y_max)
    # Set the Ticks of x and y axis
    # plt.xticks(hsk_coords)

    ax1.set_xticks(hsk_coords)
    ax1.set_xticklabels(plot_hsk_symbols)
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(cut_energy[1], y_max)
    ax1.set_yticks([0.1, 0.15, 0.2])

    ax2.set_xlim(x_min, x_max)
    ax2.set_ylim(y_min, cut_energy[0])
    ax2.set_yticks([-0.2, -0.15, -0.1])

    for ax in [ax1, ax2]:
        # Plot the solid lines for High symmetic k-points
        for kpath_i in range(kpath_num+1):
            #plt.vlines(hsk_coords[kpath_i], y_min, y_max, colors="black", linewidth=0.7)
            ax.vlines(hsk_coords[kpath_i], y_min, y_max, colors="black", linewidth=0.7)

        # Plot the fermi energy surface with a dashed line
        # plt.hlines(0.0, x_min, x_max, colors="black", linestyles="dashed", linewidth=0.7)
        ax.hlines(0.0, x_min, x_max, colors="black", linestyles="dashed", linewidth=0.7)

        # Plot the Band Structure
        for band_i in range(band_num_each_spin):
            x = kpoints_coords
            y = spin_up_energys[band_i]
            ax.plot(x, y, ls='-', color='grey', linewidth=1.5)

        # Plot the hightlight Band Structure
        for i in range(highlight_num):
            band_i = highlight_band[i]
            x = kpoints_coords
            y = spin_up_energys[band_i]
            ax.plot(x, y, ls='-', color=colorlist[i], linewidth=2)

        if spin_num == 2:
            for band_i in range(band_num_each_spin):

                x = kpoints_coords
                y = spin_dn_energys[band_i]
                ax.plot(x, y, '-', color='#0564c3', linewidth=1)  # color='#0564c3'

    # 创建轴断刻度线，d用于调节其偏转角度
    d = 0.5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                  linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
    ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

    #plt.ylabel('Energy (eV)')

    # Save the figure
    plt.tight_layout()
    plt.savefig(savepath+'/band_reconnected_brokeny.png', format=plot_format, dpi=plot_dpi, transparent=True)
    plt.savefig(savepath+'/band_reconnected_brokeny.svg', transparent=True)
    return


def compare_plot(band_data_obj, plot_args, savepath='./', auto_align=True, manual_shift=None,
                 highlight_band=[], use_scatter=False):
    '''band plot function
    the input band_data_obj should be RcBandData
    manual align is to align the reconnected band to the DFT band, which is prior to auto shift
    hightlight_band gives the band index of RC band to highlight(rg. flatband)
    '''
    shift = 0

    if 'usescatter' in plot_args.keys():
        sc_data = band_data_obj.sc_data["spin_up_energys"]
    if 'shift_to_match_dft' in plot_args.keys():
        auto_align = False
        shift = plot_args['shift_to_match_dft']

    hsk_coords = band_data_obj.band_data["hsk_coords"]
    plot_hsk_symbols = band_data_obj.band_data["plot_hsk_symbols"]
    kpath_num = band_data_obj.band_data["kpath_num"]
    band_num_each_spin = band_data_obj.band_data["band_num_each_spin"]
    kpoints_coords = band_data_obj.band_data["kpoints_coords"]
    spin_num = band_data_obj.band_data["spin_num"]
    spin_up_energys = band_data_obj.band_data["spin_up_energys"]
    spin_dn_energys = band_data_obj.band_data["spin_dn_energys"]
    min_plot_energy = plot_args['min_plot_energy']
    max_plot_energy = plot_args['max_plot_energy']

    dft_data = band_data_obj.dft_data["spin_up_energys"]

    # for the fermi energy unmatch, here we have no choice but to mannually align it
    # we will search for the matched band near homo, to avoid band index error, 7 bands will be checked
    # this works only if the kmesh are the same
    if auto_align:
        diff_list = []
        dft_homo = dft_data[band_data_obj.dft_data['homo_band_index']][:]
        for i in range(-3, 4):
            rcband_homo = spin_up_energys[band_data_obj.band_data['homo_band_index']+i][:]
            sum_abs_diff = 0
            for x1, x2 in zip(dft_homo, rcband_homo):
                sum_abs_diff += abs(x1 - x2)
            diff_list.append(sum_abs_diff)
        match_index = np.argmin(np.array(diff_list)) - 3
        rcband_homo = spin_up_energys[band_data_obj.band_data['homo_band_index']+match_index][:]
        shift = sum([item[0]-item[1] for item in zip(dft_homo, rcband_homo)])
        shift = shift/len(dft_homo)
        if abs(shift) > 0.5:
            shift = 0.0

    # Design the Figure
    # For GUI less server
    plt.switch_backend('agg')
    plt.rcParams.update({'font.size': 10,
                         'font.family': 'Arial',
                         'mathtext.fontset': 'cm'})
    # Set the spacing between the axis and labels
    plt.rcParams['xtick.major.pad'] = '5'
    plt.rcParams['ytick.major.pad'] = '5'
    # Set the ticks 'inside' the axis
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    
    
    x_min = 0.0
    x_max = hsk_coords[-1]
    y_min = min_plot_energy
    y_max = max_plot_energy
    yticks_candidate = list(np.arange(-10,10,0.5))

    def plot_use_ax(ax,y_min,y_max):
      # Set the range of plot
      ax.set_xlim(x_min, x_max)
      ax.set_ylim(y_min, y_max)
      yticks = [ i for i in yticks_candidate if i <= y_max and i>= y_min]
      # Set the label of x and y axis
      # ax.ylabel('Energy (eV)')
      # Set the Ticks of x and y axis
      ax.set_xticks(hsk_coords)
      ax.set_xticklabels(plot_hsk_symbols)
      ax.set_yticks(yticks)
      # Plot the solid lines for High symmetic k-points
      for kpath_i in range(kpath_num+1):
          ax.vlines(hsk_coords[kpath_i], y_min, y_max, colors="black", linewidth=0.7)
      # Plot the fermi energy surface with a dashed line
      ax.hlines(0.0, x_min, x_max, colors="black",
                linestyles="dashed", linewidth=0.7)

      # Plot the DFT Band Structure
      dft_band_num = band_data_obj.dft_data["band_num_each_spin"]
      dft_kpoints_coords = band_data_obj.dft_data["kpoints_coords"]
      for band_i in range(dft_band_num):
          x = dft_kpoints_coords
          y = dft_data[band_i]
          ax.plot(x, y, 'r-', linewidth=1.0)

      # Plot the Reconnected Band Structure
      for band_i in range(band_num_each_spin):
          x = kpoints_coords
          y = spin_up_energys[band_i] + shift
          ax.scatter(x, y, c='royalblue', s=1.3, zorder=1e4)

      # Plot the highlighted band structure
      for band_i in range(band_num_each_spin):
          if band_i in highlight_band:
            x = kpoints_coords
            y = spin_up_energys[band_i] + shift
            ax.scatter(x, y, c='orange', s=1.8, zorder=1e5)
    
    if 'break_y' in plot_args.keys():
      fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5.5, 5.5), sharex=True)
      ax1.spines.bottom.set_visible(False)
      ax2.spines.top.set_visible(False)
      plot_use_ax(ax1, y_min=plot_args['break_y'], y_max=y_max)
      plot_use_ax(ax2, y_min=y_min, y_max=-plot_args['break_y'])
      fig.text(0.01, 0.5, 'Energy (eV)', fontsize=12, va='center', rotation='vertical')
      
      # 创建轴断刻度线，d用于调节其偏转角度
      d = 0.5  # proportion of vertical to horizontal extent of the slanted line
      kwargs = dict(marker=[(-1, -d), (1, d)], markersize=8,
                    linestyle="none", color='k', mec='k', mew=1, clip_on=False)
      ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
      ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
    else:
      # Create the figure and axis object
      fig = plt.figure(figsize=(5.5, 5.5))
      ax = fig.add_subplot(1, 1, 1)
      plot_use_ax(ax, y_min=y_min, y_max=y_max)
      fig.text(0.01, 0.5, 'Energy (eV)', fontsize=12, va='center', rotation='vertical')

      

    # Save the figure
    # plt.tight_layout()
    plt.savefig(savepath+'/band_compare.png', dpi=600, transparent=True)
    plt.savefig(savepath+'/band_compare.svg', transparent=True)
    plt.close()
    return shift


def band():
    '''band functions'''
    plot_args = get_command_line_input()
    band_data_obj = BandData(plot_args["data_type"])
    band_data_obj.file_read(plot_args["data_filename"])
    band_data_obj.get_band_data()
    if not plot_args["no_plot"]:
        band_plot(band_data_obj, plot_args)


def load_process_band(readpath, savepath, refine_fermi=True, get_bandwidth=True,
                      advise_plot_range=True, plot=True, compare_dft=True, dft_path=None,
                      usescatter=False, egval_path='egval_k/', plot_para={}):
    '''
    by baot:
    load and process data from reconnected json
    add path argument for file input and output
    the options are not indepedent, see the source code
    When usescatter is not none, the scatter path should be offered, eg. egval_k/
    '''
    # load plot_args, the apointed will be prior to get_cmd
    plot_args = get_command_line_input()
    for key in plot_para.keys():
        plot_args[key] = plot_para[key]

    band_data_obj = RcBandData(plot_args["data_type"])
    with open(readpath+'/band_reconnect.json', 'r') as f:
        data_new = json.load(f)
    for key, val in data_new.items():
        if type(val) is list:
            data_new[key] = np.array(val)
    band_data_obj.band_data = data_new

    if refine_fermi:
        # baot: make sure refine fermi energy works
        band_data_obj.refined_reconnected_fermi_energy()
    if get_bandwidth:
        band_data_obj.get_reconnected_bandwidth()
    if advise_plot_range:
        plot_args = band_data_obj.advise_plot_range(plot_args=plot_args, cover_band=15)
    # ==
    if not plot_args["no_plot"] and plot:
        band_plot(band_data_obj, plot_args, savepath=savepath)
        # band_plot_P4(band_data_obj, plot_args, cut_energy=[-0.09,0.09], highlight_num= 2 )

    shift = None
    if compare_dft:
        dft_data_obj = BandData(plot_args["data_type"])
        dft_data_obj.file_read(dft_path+'/openmx.Band')
        dft_data_obj.get_band_data()
        band_data_obj.dft_data = dft_data_obj.band_data
        # give the dft data to the reconnected data obj
        if usescatter:
            sc_data_obj = ScBandData(plot_args["data_type"])
            sc_data_obj.get_band_data(egval_path)  # usescatter = 'egval_k/'
            band_data_obj.sc_data = sc_data_obj.band_data
        if len(band_data_obj.band_data['cross_fermi_index']) > 0:
            temp = max(band_data_obj.band_data['cross_fermi_index'])
        else:
            temp = band_data_obj.band_data["homo_band_index"]
        highlight = [i for i in band_data_obj.band_data['flat_index'] if i > temp-10 and i < temp+11]
        shift = compare_plot(band_data_obj, plot_args, savepath=savepath, highlight_band=highlight)

    return band_data_obj, shift
# +----------------+
# |  Main Process  |
# +----------------+


def main():
    load_process_band()
    return


if __name__ == '__main__':
    main()
