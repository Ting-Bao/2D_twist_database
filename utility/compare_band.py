# modified by TingBAO @ THU
#
# default settings used in xyz program:
#   mv self.mp_band_pkl_conv.json mp.Band.json
#   python compare_band.py --n1 mp.Band.json --n2 openmx.Band.json -u 6.0 -d -6.0 -f png
# 
# 注意文件名中首个为蓝色，第二个为红色


# Author: Liyang@CMT
# Date: 2019.7.1
# Description: This python code is designed for compare two different band.
#              ::Input File::
#                - *.Band.json
#                - *.Band.json
#              ::Output File::
#                - *.Band.compare.eps
#                - *.Band.compare.json
# Usage: python3.7 band_plot.py -u <max_energy> -d <min_energy>
#

import os
import json
import fnmatch
import sys
import getopt
import math
import numpy as np
import re
import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt


def read_command_line(argv, kDefaultPlotEnergyRange,
                      kDefaultPlotFormat,
                      kDefaultDPI):
    band_data_file_1 = ''
    band_data_file_2 = ''
    export_band_plot = True
    min_plot_energy = kDefaultPlotEnergyRange[0]
    max_plot_energy = kDefaultPlotEnergyRange[1]
    plot_format = kDefaultPlotFormat
    plot_dpi = kDefaultDPI
    try:
        opts, args = getopt.getopt(argv, "hd:u:f:r:b",
                                   ["n1=", "n2=", "min=", "max=",
                                    "format=", "dpi=", "no-plot"])
    except getopt.GetoptError:
        print('band_plot.py --n1 <band1> --n2 <band2> -f <PlotFormat>')
        sys.exit(2)
    del args
    for opt, arg in opts:
        if opt == '-h':
            print('band_plot.py --n1 <band1> --n2 <band2> -f <PlotFormat>')
            sys.exit()
        elif opt in ("--n1"):
            band_data_file_1 = arg
        elif opt in ("--n2"):
            band_data_file_2 = arg
        elif opt in ("-d", "--min"):
            min_plot_energy = float(arg)
        elif opt in ("-u", "--max"):
            max_plot_energy = float(arg)
        elif opt in ("-f", "--format"):
            plot_format = arg
        elif opt in ("-r", "--dpi"):
            plot_dpi = int(arg.split(".")[0])
        elif opt in ("-b", "--no-plot"):
            export_band_plot = False
        else:
            pass
    return band_data_file_1, band_data_file_2, min_plot_energy, \
        max_plot_energy, plot_format, plot_dpi, export_band_plot

def read_json_file(band_data_file):
    with open(band_data_file, 'r') as jfrp:
        band_data = json.load(jfrp)
    return band_data

def calc_std_deviation(x1, x2, y1, y2):
    variance = 0.0
    if (len(x1) > len(x2)):
        y2 = np.interp(x1, x2, y2)
    else:
        y1 = np.interp(x2, x1, y1)
    for index in range(len(y1)):
        variance += abs(y1[index] - y2[index])
    std_deviation = variance / len(y1)
    return std_deviation

def main(argv):
    ############################
    ### Some useful constant ###
    ############################
    
    kDefaultPlotEnergyRange = [-2, 2]
    kDefaultPlotFormat = 'png'
    kDefaultDPI = 1000

    ##################################################
    ### Basic Parameters read in from command line ###
    ##################################################
    band_data_file_1, band_data_file_2, min_plot_energy, \
        max_plot_energy, plot_format, plot_dpi, export_band_plot\
        = read_command_line(argv, kDefaultPlotEnergyRange,
                            kDefaultPlotFormat,
                            kDefaultDPI)

    ##############################
    ### Band Data File Read In ###
    ##############################
    # Get the filename contants the band data
    # If the user do not appoint any data file, then find one!
    if (band_data_file_1 == ''):
        for filename in os.listdir('.'):
            if fnmatch.fnmatch(filename, '*.Band.json') and \
               (filename != band_data_file_2):
                band_data_file_1 = filename
                del filename
                break
    if (band_data_file_2 == ''):
        for filename in os.listdir('.'):
            if fnmatch.fnmatch(filename, '*.Band.json') and \
               (filename != band_data_file_1):
                band_data_file_2 = filename
                del filename
                break
    # If do not find any valid file, then
    if (band_data_file_1 == '') or (band_data_file_2 == ''):
        print("[error] Can not find enough *.Band.json file in current folder!")
        exit(0)
    band_mark_1 = '.'.join(band_data_file_1.split('.')[:-2])
    band_mark_2 = '.'.join(band_data_file_2.split('.')[:-2])
    ## Read File Data
    # File 1 : Read Json
    json_band_data_1 = read_json_file(band_data_file_1)
    # File 1 : Get necessary data
    homo_order_1 = json_band_data_1['homo_band_order']
    homo_index_1 = homo_order_1 - 1
    kpath_quantity_1 = json_band_data_1['kpath_quantity']
    kpoints_quantity_1 = json_band_data_1['kpoints_quantity']
    kpoints_corrdinate_list_1 = json_band_data_1['kpoints_corrdinate_list']
    hsk_corrdinate_list_1 = json_band_data_1['hsk_corrdinate_list']
    hsk_symbol_list_1 = json_band_data_1['hsk_symbol_list']
    band_quantity_1 = json_band_data_1['band_quantity']
    band_data_1 = json_band_data_1['sorted_band_data']
    homo_band_1 = band_data_1[homo_index_1]
    lumo_band_1 = band_data_1[homo_index_1+1]
    # File 2 : Read Json
    json_band_data_2 = read_json_file(band_data_file_2)
    # File 2 : Get necessary data
    homo_order_2 = json_band_data_2['homo_band_order']
    homo_index_2 = homo_order_2 - 1
    #kpath_quantity_2 = json_band_data_2['kpath_quantity']
    kpoints_quantity_2 = json_band_data_2['kpoints_quantity']
    kpoints_corrdinate_list_2 = json_band_data_2['kpoints_corrdinate_list']
    hsk_corrdinate_list_2 = json_band_data_2['hsk_corrdinate_list']
    #hsk_symbol_list_2 = json_band_data_2['hsk_symbol_list']
    band_quantity_2 = json_band_data_2['band_quantity']
    band_data_2 = json_band_data_2['sorted_band_data']
    homo_band_2 = band_data_2[homo_index_2]
    lumo_band_2 = band_data_2[homo_index_2+1]

    #using for compare mp ,vasp and openmx, may not be universal for any band compare
    if band_data_file_1 == 'openmx.Band.json' or band_data_file_1 ==  'vasp.Band.json':
        kpoints_corrdinate_list_1=list(map(lambda x: x*2*np.pi,kpoints_corrdinate_list_1))
        hsk_corrdinate_list_1=list(map(lambda x: x*2*np.pi,hsk_corrdinate_list_1))
        print("change openmx k coordinate to std_form")
    if band_data_file_2 == 'openmx.Band.json' or band_data_file_1 == 'vasp.Band.json':
        kpoints_corrdinate_list_2=list(map(lambda x: x*2*np.pi,kpoints_corrdinate_list_2))
        hsk_corrdinate_list_2=list(map(lambda x: x*2*np.pi,hsk_corrdinate_list_2))
        print("change openmx k coordinate to std_form")
    #compare kpath setting and rlv
    print(hsk_corrdinate_list_1,' ',band_data_file_1)
    print(hsk_corrdinate_list_2,' ',band_data_file_2)
    error = np.max(np.abs(np.array(hsk_corrdinate_list_1)-np.array(hsk_corrdinate_list_2)))
    assert error < 2.5e-2, "kpath of mp and refine_str are not coincident" + str(error)

    ####################
    ### Band Compare ###
    ####################
    homo_diff = calc_std_deviation(kpoints_corrdinate_list_1,
                                   kpoints_corrdinate_list_2,
                                   homo_band_1,
                                   homo_band_2)
    lumo_diff = calc_std_deviation(kpoints_corrdinate_list_1,
                                   kpoints_corrdinate_list_2,
                                   lumo_band_1,
                                   lumo_band_2)
    #################
    ### Band Plot ###
    #################
    # Design the Figure
    if export_band_plot:
        # Special Symbol List
        greek_symbol_list = ['Gamma','Delta','Theta','Lambda','Xi',
                             'Pi','Sigma','Phi','Psi','Omega']
        plot_hsk_symbol_list = []
        for symbol in hsk_symbol_list_1:
            symbol = symbol.replace("\\", "")
            for greek_symbol in greek_symbol_list:
                latex_greek_symbol = "$\\" + greek_symbol + "$"
                symbol = re.sub(greek_symbol, "orz", symbol, 
                                flags=re.I)
                symbol = symbol.replace("orz", latex_greek_symbol)
            symbol = re.sub(r'_\d+', lambda x:'$'+x[0]+'$', symbol)
            plot_hsk_symbol_list.append(symbol)
        # Set the DPI and reslution of the figure
        # >>plt.rcParams['savefig.dpi'] = plot_dpi
        # >>plt.rcParams['figure.dpi'] = plot_dpi
        # Set the Fonts
        plt.rcParams.update({'font.size': 14,
                             'font.family': 'STIXGeneral',
                             'mathtext.fontset': 'stix'})
        # Set the spacing between the axis and labels
        plt.rcParams['xtick.major.pad'] = '6'
        plt.rcParams['ytick.major.pad'] = '6'
        # Set the ticks 'inside' the axis
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        # Create the figure and axis object
        fig = plt.figure()
        band_plot = fig.add_subplot(1, 1, 1)
        # Set the range of plot
        x_min = 0.0
        x_max = hsk_corrdinate_list_1[-1]
        y_min = min_plot_energy
        y_max = max_plot_energy
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        # Set the label of x and y axis
        # >>font = {'family': 'Arial', 'weight': 'normal', 'size': 18}
        plt.xlabel('')
        plt.ylabel('Energy (eV)')
        # Set the Ticks of x and y axis
        plt.xticks(hsk_corrdinate_list_1)
        band_plot.set_xticklabels(plot_hsk_symbol_list)
        # >>plt.tick_params(labelsize=18)
        # >>labels = band_plot.get_xticklabels() + band_plot.get_yticklabels()
        # >>for label in labels:
        # >>    label.set_fontname('Arial')
        # Plot the solid lines for High symmetic k-points
        for kpath_index in range(kpath_quantity_1):
            plt.vlines(hsk_corrdinate_list_1[kpath_index],
                       y_min, y_max, colors="black", linewidth=0.7)
        # Plot the fermi energy surface with a dashed line
        plt.hlines(0.0, x_min, x_max, colors="black",
                   linestyles="dashed", linewidth=0.7)
        # Plot the Band Structure
        x1 = kpoints_corrdinate_list_1
        x2 = kpoints_corrdinate_list_2
        for band_index in range(band_quantity_1):
            y = band_data_1[band_index]
            band_plot.plot(x1, y, 'b-', linewidth=1.4)
        for band_index in range(band_quantity_2):
            y = band_data_2[band_index]
            band_plot.plot(x2, y, 'r-', linewidth=0.8)
        # Plot the band structure
        plot_band_file_name = band_mark_1 + '-' + band_mark_2 \
            + '.compare.' + plot_format
        plt.savefig(plot_band_file_name, format=plot_format, dpi=plot_dpi)
        # >>Show plot
        # >>plt.show()
        #

    ###########################
    ### Compare Data Output ###
    ###########################
    # Store Data Json 
    stored_data_josn = {
        "homo_diff" : homo_diff,
        "lumo_diff" : lumo_diff,
        "kpoints_quantity" : {
            "1" : kpoints_quantity_1,
            "2" : kpoints_quantity_2,
        },
        "band_quantity" : {
            "1" : band_quantity_1,
            "2" : band_quantity_2,
        },
        "homo_order" : {
            "1" : homo_order_1,
            "2" : homo_order_2,
        },
    }
    # Define the data filename and save the data
    band_compare_diff_file = band_mark_1 + '-' + band_mark_2 + '.compare.json'
    with open(band_compare_diff_file, 'w') as jfwp:
        json.dump(stored_data_josn, jfwp)

if __name__ == "__main__":
    main(sys.argv[1:])
