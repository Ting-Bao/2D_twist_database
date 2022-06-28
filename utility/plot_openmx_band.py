# Author: Liyang@CMT
# Date: 2020.4.17
# Description: This python code is designed for calculate the fermi energy
#                and plot the band structure use the *.Band file.
#              ::Input File::
#                - *.Band
#              ::Output File::
#                - *.Band.png
#                - *.Band.conv
# Usage: python3.7 band_plot.py -u <max_energy> -d <min_energy>
#

import os
import fnmatch
import sys
import getopt
import math
import numpy as np
import re
import matplotlib.pyplot as plt
import json

def band_data_sort(band_data):
    #kpoints_quantity = len(band_data[0])
    #average_band_energy_array = \
    #    np.sum(band_data, axis=1) / kpoints_quantity
    #return band_data[np.argsort(average_band_energy_array)]
    kpoints_quantity = len(band_data[0])
    band_quantity = len(band_data)
    sorted_band_data = np.zeros((band_quantity, kpoints_quantity))
    for kpoints_index in range(kpoints_quantity):
        current_k_column = band_data[:, kpoints_index]
        sorted_current_k_column = np.sort(current_k_column)
        sorted_band_data[:, kpoints_index] = sorted_current_k_column
    return sorted_band_data

def cal_k_distance(local_rlv, start_kpoint_frac, end_kpoint_frac):
    # Calculate the distance of two different kpoints in k-space
    # Reciprocal lattice vector but localy uesd :: local_rlv
    # Start kpoints in frac k-cooridnate :: start_kpoint
    # End kpoints in frac k-cooridnate :: end_kpoint
    # 1 bohr = 0.53 angstrom
    one_bohr = 0.529177210903 #A

    # Calcualte the k-cooridnate in Cartesian indicator
    start_kpoints_cart = np.array([0.0, 0.0, 0.0])
    end_kpoints_cart = np.array([0.0, 0.0, 0.0])
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
        start_kpoints_cart[xyz] = 0.0
        for b_index in range(3):
            start_kpoints_cart[xyz] += start_kpoint_frac[b_index] * \
                local_rlv[b_index, xyz]
            end_kpoints_cart[xyz] += end_kpoint_frac[b_index] * \
                local_rlv[b_index, xyz]
    # Calculate the k distance of the two kpoints
    k_distance = math.sqrt(sum((start_kpoints_cart-end_kpoints_cart) ** 2))
    k_distance /= one_bohr
    k_distance /= (2*math.pi)
    return k_distance

def main(argv):
    ############################
    ### Some useful constant ###
    ############################
    # 1Ha = 27.21138602eV
    kUnitScale = 27.21138602
    kDefaultPlotEnergyRange = [-6, 6]
    kDefaultPlotFormat = 'png'
    kDefaultDPI = 1000

    ##################################################
    ### Basic Parameters read in from command line ###
    ##################################################
    band_data_file = ''
    export_band_plot=True
    min_plot_energy = kDefaultPlotEnergyRange[0]
    max_plot_energy = kDefaultPlotEnergyRange[1]
    plot_format = kDefaultPlotFormat
    plot_dpi = kDefaultDPI
    try:
        opts, args = getopt.getopt(argv, "hn:d:u:f:r:b",
                                   ["name=", "min=", "max=",
                                    "format=", "dpi=","no-plot"])
    except getopt.GetoptError:
        print('band_plot.py -d <E_min> -u <E_max> -f <PlotFormat>')
        sys.exit(2)
    del args
    for opt, arg in opts:
        if opt == '-h':
            print('band_plot.py -d <E_min> -u <E_max> -f <PlotFormat>')
            sys.exit()
        elif opt in ("-n", "--name"):
            band_data_file = arg
        elif opt in ("-d", "--min"):
            min_plot_energy = float(arg)
        elif opt in ("-u", "--max"):
            max_plot_energy = float(arg)
        elif opt in ("-f", "--format"):
            plot_format = arg
        elif opt in ("-r", "--dpi"):
            plot_dpi = int(arg.split(".")[0])
        elif opt in ("-b", "--no-plot"):
            export_band_plot=False
        else:
            pass

    ##############################
    ### Band Data File Read In ###
    ##############################
    # Get the filename contants the band data
    # If the user do not appoint any data file, then find one!
    if band_data_file == '':
        for filename in os.listdir('.'):
            if fnmatch.fnmatch(filename, '*.Band'):
                band_data_file = filename
                del filename
                break
    # If do not find any valid file, then
    if band_data_file == '':
        print("[error] Can not find *.Band file in current folder!")
        sys.exit(0)
    # Open & Read file
    read_line_index = -1  # Read line index mark
    with open(band_data_file, 'r') as frp:
        str_file_lines = frp.readlines()
    # Read the 1st line :: Basic info.
    read_line_index += 1
    line = str_file_lines[read_line_index].split()
    band_quantity_each_spin = int(line[0])                     
    spin_num = int(line[1]) + 1                        # Spin Number
    if spin_num not in [1, 2]:
        print('[error] Spin Number Must be 1 or 2!!!')
        sys.exit(0)
    band_quantity = band_quantity_each_spin * spin_num # Band Quantity
    fermi_energy = float(line[2]) * kUnitScale         # Fermi Energy
    # Read the 2nd line :: rlv
    read_line_index += 1
    line = str_file_lines[read_line_index].split()
    # Reciprocal lattice vector :: rlv
    # Each row of the 'rlv' matrix is a basis vector
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
    # Read the 3rd line :: Kpath Number
    read_line_index += 1
    line = str_file_lines[read_line_index].split()
    kpath_quantity = int(line[0])
    # Read the K-Path
    kpath_density_list = [1 for path in range(kpath_quantity)]
    kpath_symbol_list = [[] for path in range(kpath_quantity)]
    kpath_vector_list = [[[], []] for path in range(kpath_quantity)]
    for kpath_index in range(kpath_quantity):
        read_line_index += 1
        line = str_file_lines[read_line_index].split()
        kpath_density_list[kpath_index] = int(line[0])
        kpath_vector_list[kpath_index][0] = [
            float(line[1]), float(line[2]), float(line[3])]
        kpath_vector_list[kpath_index][1] = [
            float(line[4]), float(line[5]), float(line[6])]
        kpath_symbol_list[kpath_index] = [line[7], line[8]]
    # Read band data
    # Establish the band data matrix using numpy array
    #   total kpoints number & Band Quantity are needed
    kpoints_quantity = sum(kpath_density_list)
    spin_up_band_data = np.zeros((band_quantity_each_spin, kpoints_quantity))
    if spin_num == 2:
        spin_down_band_data = \
            np.zeros((band_quantity_each_spin, kpoints_quantity))
    kpoints_vector_list = [[] for kpoints in range(kpoints_quantity)]
    # Record the band data started line
    data_start_line_index = read_line_index + 1
    # Re-initialize the line pointer
    read_line_index = data_start_line_index - 1
    for kpoints_index in range(kpoints_quantity):
        # Read Kpoints Vec. (spin up)
        read_line_index += 1
        line = str_file_lines[read_line_index].split()
        kpoints_vector_list[kpoints_index] = [float(line[1]), 
                                              float(line[2]), 
                                              float(line[3])]
        # Read Kpoints band. (spin up)
        read_line_index += 1
        line = str_file_lines[read_line_index].split()
        for band_index in range(band_quantity_each_spin):
            spin_up_band_data[band_index, kpoints_index] = \
                float(line[band_index]) * kUnitScale
        # If has the spin down
        if spin_num == 2: 
            # Read Kpoints band. (spin down)
            read_line_index += 2
            line = str_file_lines[read_line_index].split()
            for band_index in range(band_quantity_each_spin):
                spin_down_band_data[band_index, kpoints_index] = \
                    float(line[band_index]) * kUnitScale
    # Combine the band
    if spin_num == 1:
        sorted_band_data = spin_up_band_data
    else:
        sorted_band_data = \
            np.concatenate((spin_up_band_data, spin_down_band_data),axis=0)
    # Sort the band data
    sorted_band_data = band_data_sort(sorted_band_data)
    #
    # Summary [Part I] #
    # After the file reading, we get:
    # - Band Quantity                                :: band_quantity
    # - Spin Number                                  :: spin_num
    # - Band Quantity Each Spin                      :: band_quantity_each_spin
    # - Fermi Energy                                 :: fermi_energy
    # - Reciprocal lattice vector                    :: rlv
    #              +-              -+
    #              | b_1x b_1y b_1z |
    #        rlv = | b_2x b_2y b_2z |
    #              | b_3x b_3y b_3z |
    #              +-              -+
    # - Number of the higy symmetric k-path          :: kpath_quantity
    # - List of the number of kpoints in each k-path :: kpath_density_list
    # - Start and end points for each k-path         :: kpath_vector_list
    # - List of high symmetric points symbol         :: kpath_symbol_list
    # - Number of the kpoints in all Kpath           :: kpoints_quantity
    # - Coordinate of each k-points on the k-path    :: kpoints_vector_list
    # - Band energy of each kpoints on each band     :: sorted_band_data
    #                                    band_index(n)
    #                          +-                  -+
    #                          | E_11 E_12 ... E_1k |  kpoints_index(k)
    #       sorted_band_data = | E_21 E_22 ... E_2k |
    #                          | .             .    |
    #                          | E_n1 E_n2 ... E_nk |
    #                          +-                  -+
    # - Spin up Band data                         :: spin_up_band_data
    #        (hsa the same form as 'sorted_band_data')
    # - Spin down Band data (if exist, ISPIN==2)  :: spin_down_band_data
    #        (hsa the same form as 'sorted_band_data') 
    #

    ####################
    ### Data Process ###
    ####################
    # Get the length of each k-path in k-space
    hsk_distance_list = [0.0 for kpath_index in range(kpath_quantity)]
    sum_hsk_distance_list = [0.0 for kpath_index in range(kpath_quantity)]
    for kpath_index in range(kpath_quantity):
        start_hsk = kpath_vector_list[kpath_index][0]
        end_hsk = kpath_vector_list[kpath_index][1]
        hsk_distance_list[kpath_index] = cal_k_distance(rlv, start_hsk, end_hsk)
        sum_hsk_distance_list[kpath_index] = \
            sum(hsk_distance_list[0:kpath_index+1])
    hsk_corrdinate_list = [0.0] + sum_hsk_distance_list
    hsk_corrdinate_list = np.array(hsk_corrdinate_list)
    # Get the distance in k-space of k-points on the k-path
    kpoints_corrdinate_list = np.zeros(kpoints_quantity)
    kpoints_index = -1
    for kpath_index in range(kpath_quantity):
        # Count the Previous kpath distance
        pre_path_distance = hsk_corrdinate_list[kpath_index]
        # Calculate the kpoints' distance in current kpath
        for kpath_kpoints_index in range(kpath_density_list[kpath_index]):
            kpoints_index += 1
            start_hsk = kpath_vector_list[kpath_index][0]
            end_hsk = kpoints_vector_list[kpoints_index]
            # The total distance equals to (pre_dis + curr_dis)
            kpoints_corrdinate_list[kpoints_index] = \
                pre_path_distance + cal_k_distance(rlv, start_hsk, end_hsk)
    del kpath_kpoints_index
    # Calculate the Fermi Energy
    conut_under_fermi_list = [0 for kpoints_index in range(kpoints_quantity)]
    for kpoints_index in range(kpoints_quantity):
        conut_under_fermi = 0
        for band_index in range(band_quantity):
            if (fermi_energy > sorted_band_data[band_index][kpoints_index]):
                conut_under_fermi += 1
        conut_under_fermi_list[kpoints_index] = conut_under_fermi
    # Count the appearance of each element in the list
    ele_count = {} # {element:appearing_times}
    for conut_under_fermi in conut_under_fermi_list:
        ele_count[conut_under_fermi] = ele_count.get(conut_under_fermi, 0) + 1
    # Judge the distribution of the band under fermi level
    if len(ele_count) == 1:
        is_metal = False
    else:
        appear_times_list = []
        for key in ele_count.keys():
            appear_times_list.append(ele_count[key])
        confirm_number = 0
        for appear_times in appear_times_list:
            check_number = int(0.05*kpoints_quantity/len(ele_count)) + 1
            if appear_times > check_number:
                confirm_number += 1
        if confirm_number == 0:
            print("[error] Band kpoints under fermi counting error...")
            sys.exit(1)
        elif confirm_number == 1:
            is_metal = False
        else:
            is_metal = True
    # If is a insulator:
    # - the fermi energy is inside the gap
    # - the fermi energy is on the LUMO or HOMO
    if not is_metal:
        max_val = -1
        key_max = -1
        for key, val in ele_count.items():
            if val > max_val:
                key_max = key
                max_val = val
        if key_max == -1:
            print("[error] homo finding error!!!")
            sys.exit()
        homo_index = key_max - 1
        lumo_index = homo_index + 1
        homo_energy = max(sorted_band_data[homo_index])
        lumo_energy = min(sorted_band_data[lumo_index])
        fermi_energy = (homo_energy + lumo_energy) / 2
        band_gap = lumo_energy - homo_energy
    # If is a conductor
    else:
        band_gap = 0.0
        homo_index = band_quantity - 2
        band_k_average=np.sum(sorted_band_data, axis=1) / kpoints_quantity
        for band_index in range(band_quantity):
            if (fermi_energy <= band_k_average[band_index]):
                lumo_index = band_index
                homo_index = lumo_index - 1
                break
    # Translation the Energy to E_f = 0
    sorted_band_data = sorted_band_data - fermi_energy
    spin_up_band_data -= fermi_energy
    if spin_num == 2:
        spin_down_band_data -= fermi_energy
    #
    # Summary [Part II] #
    # After the Data Process, we get:
    # - Kpoints distance list                   :: kpoints_corrdinate_list
    # - Sum of previous kpath nodes distance    :: sum_hsk_distance_list
    # - Fermi Energy                            :: fermi_energy
    # - The High symmetric k-points' coordiante :: hsk_corrdinate_list
    # - An energy translated matrix             :: sorted_band_data
    # Use it together with sorted_band_data, we can finally plot
    #   the band structure.
    #

    #################
    ### Band Plot ###
    #################
    # Prepare for the plot
    # Special Symbol List
    greek_symbol_list = ['Gamma','Delta','Theta','Lambda','Xi',
                         'Pi','Sigma','Phi','Psi','Omega']
    # Prepare the symbol of k-axis (xtics)
    hsk_symbol_list = ['' for kpath_index in range(kpath_quantity+1)]
    hsk_symbol_list[0] = kpath_symbol_list[0][0]
    for kpath_index in range(1, kpath_quantity):
        if kpath_symbol_list[kpath_index][0] == \
           kpath_symbol_list[kpath_index-1][1]:
            hsk_symbol_list[kpath_index] = \
                kpath_symbol_list[kpath_index][0]
        else:
            hsk_symbol_list[kpath_index] = \
                kpath_symbol_list[kpath_index - 1][1] + \
                '|' + kpath_symbol_list[kpath_index][0]
    hsk_symbol_list[-1] = kpath_symbol_list[-1][1]
    ## Plot the Band
    if export_band_plot:
        plot_hsk_symbol_list = []
        for symbol in hsk_symbol_list:
            symbol = symbol.replace("\\", "")
            for greek_symbol in greek_symbol_list:
                latex_greek_symbol = "$\\" + greek_symbol + "$"
                symbol = re.sub(greek_symbol, "orz", symbol, 
                                flags=re.I)
                symbol = symbol.replace("orz", latex_greek_symbol)
            symbol = re.sub(r'_\d+', lambda x:'$'+x[0]+'$', symbol)
            plot_hsk_symbol_list.append(symbol)
        ## Design the Figure
        # Set the DPI and reslution of the figure
        # >>plt.rcParams['savefig.dpi'] = plot_dpi
        # >>plt.rcParams['figure.dpi'] = plot_dpi
        # Set the Fonts
        plt.rcParams.update({'font.size': 14,
                            'font.family': 'STIXGeneral',
                            'mathtext.fontset': 'stix'})
        # Set the spacing between the axis and labels
        plt.rcParams['xtick.major.pad']='6'
        plt.rcParams['ytick.major.pad']='6'
        # Set the ticks 'inside' the axis
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        # Create the figure and axis object
        fig = plt.figure()
        band_plot = fig.add_subplot(1, 1, 1)
        # Set the range of plot
        x_min = 0.0
        x_max = hsk_corrdinate_list[-1]
        y_min = min_plot_energy
        y_max = max_plot_energy
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        # Set the label of x and y axis
        # >>font = {'family': 'Arial', 'weight': 'normal', 'size': 18}
        plt.xlabel('')
        plt.ylabel('Energy (eV)')
        # Set the Ticks of x and y axis
        plt.xticks(hsk_corrdinate_list)
        band_plot.set_xticklabels(plot_hsk_symbol_list)
        # >>plt.tick_params(labelsize=18)
        # >>labels = band_plot.get_xticklabels() + band_plot.get_yticklabels()
        # >>for label in labels:
        # >>    label.set_fontname('Arial')
        # Plot the solid lines for High symmetic k-points
        for kpath_index in range(kpath_quantity+1):
            plt.vlines(hsk_corrdinate_list[kpath_index],
                    y_min, y_max, colors="black", linewidth=0.7)
        # Plot the fermi energy surface with a dashed line
        plt.hlines(0.0, x_min, x_max, colors="black",
                linestyles="dashed", linewidth=0.7)
        # Plot the Band Structure
        for band_index in range(band_quantity_each_spin):
            x = kpoints_corrdinate_list
            y = spin_up_band_data[band_index]
            band_plot.plot(x, y, 'r-', linewidth=1.2)
        if spin_num == 2:
            for band_index in range(band_quantity_each_spin):
                x = kpoints_corrdinate_list
                y = spin_down_band_data[band_index]
                band_plot.plot(x, y, '-', color='#0564c3', linewidth=0.8)
        # Save the figure
        plot_band_file_name = band_data_file + '.' + plot_format
        plt.savefig(plot_band_file_name, format=plot_format, dpi=plot_dpi)
    #
    # Summary [Part III] #
    # - High symmetric k-points' symbol list :: hsk_symbol_list
    # *.png is created!
    #

    ########################
    ### Band Data Output ###
    ########################
    # Write gamma to GAMMA
    store_hsk_symbol_list = []
    for symbol in hsk_symbol_list:
        symbol = symbol.replace("\\", "")
        for greek_symbol in greek_symbol_list:
            symbol = re.sub(greek_symbol, "orz", symbol, 
                            flags=re.I)
            symbol = symbol.replace("orz", greek_symbol)
        store_hsk_symbol_list.append(symbol)
    # Establish the stored data json block
    stored_data_josn = {
        "spin_num":spin_num,
        "rlv" : rlv.tolist(),
        "fermi_energy" : fermi_energy,
        "kpath_quantity" : kpath_quantity,
        "kpath_density_list" : kpath_density_list,
        "kpoints_quantity" : kpoints_quantity,
        "band_quantity_each_spin" : band_quantity_each_spin,
        "band_quantity" : band_quantity,
        "band_gap" : band_gap,
        "homo_band_order" : homo_index+1,
        "hsk_symbol_list" : store_hsk_symbol_list,
        "hsk_corrdinate_list" : hsk_corrdinate_list.tolist(),
        "kpath_symbol_list" : kpath_symbol_list,
        "kpath_vector_list" : kpath_vector_list,
        "kpoints_corrdinate_list" : kpoints_corrdinate_list.tolist(),
        "spin_band_data" : {},
        "sorted_band_data" : sorted_band_data.tolist(),
    }
    if spin_num == 2:
        stored_data_josn["spin_band_data"]["spin_up"] = \
            spin_up_band_data.tolist()
        stored_data_josn["spin_band_data"]["spin_down"] = \
            spin_down_band_data.tolist()
    # Save the josn data to file
    band_data_file_json = band_data_file + '.json'
    with open(band_data_file_json, 'w') as jfwp:
        json.dump(stored_data_josn, jfwp)
    #
    # Summary [Part IV] #
    # *.conv is created!
    # You can plot the Band with '*.conv' by using Gnuplot or Origin.
    #


if __name__ == "__main__":
    main(sys.argv[1:])


