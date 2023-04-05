# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.07
# this code is to reconect and calculate and do statics on band width under VBM
import json
import os
import pandas as pd
from tqdm import tqdm
from utility.connect_interpolate_func import connect_band
from utility.refineplot_reconnected_band import *

datapath = 'data/infer_calculation/'


def main():
    matlist = os.listdir(datapath)
    matlist = [i for i in matlist if os.path.isdir(datapath+i)]
    pathlist = []
    for mat in matlist:
        twlist = os.listdir(os.path.join(datapath, mat))
        twlist = [i for i in twlist if os.path.isdir(os.path.join(datapath, mat, i))]
        for tw in twlist:
            path = os.path.join(datapath, mat, tw, 'band')
            pathlist.append(path)

    data = pd.DataFrame(columns=['Material formula', 'C2DB id', 'Twisted index', '#(half) occupied flat bands',
                                 'occbandwidth', '#occupied flat bands', 'uoccbandwidth'])
    index = 0
    for path in tqdm(pathlist):#
        print(path)
        connect_band(inputpath=path, outputpath=path, spline_interpolate=False)
        # write band_reconnected.json, the fermi energy need further refinement (implemented above)
        band_obj = load_process_band(readpath=path, savepath=path, dft_path=path+'/../dft_band/')

        # count the flat band with condition:
        # at most 10 (half) occupied bands
        # at most 10 unoccupied bands
        if len(band_obj.band_data['cross_fermi_index'])>0:
            temp = max(band_obj.band_data['cross_fermi_index'])
        else:
            temp = band_obj.band_data["homo_band_index"]
        uocp_flat_index = [i for i in band_obj.band_data['flat_index'] if i > temp and i < temp+11]
        ocp_flat_index = [i for i in band_obj.band_data['flat_index'] if i > temp-11 and i <= temp]
        if len(uocp_flat_index)+len(ocp_flat_index) == 0:
            continue
        else:
            ocp_bandwidth = [str(int(1000*band_obj.band_data['bandwidth_info'][i]['width'])) for i in ocp_flat_index]
            uocp_bandwidth = [str(int(1000*band_obj.band_data['bandwidth_info'][i]['width'])) for i in uocp_flat_index]

            matinfo = path.split('/')
            formula = ''
            for i in matinfo[2].split('-')[0]:
                if i.isdigit():
                    formula += '$_{}$'.format(i)
                else:
                    formula += i
            if len(ocp_bandwidth)==0:
                ocp_bw_info='N/A'
            else:
                ocp_bw_info=', '.join(ocp_bandwidth)
            if len(uocp_bandwidth)==0:
                uocp_bw_info='N/A'
            else:
                uocp_bw_info=', '.join(uocp_bandwidth)

            data.loc[index] = {'Material formula': formula,
                               'C2DB id': matinfo[2].split('-')[1], 
                               'Twisted index':matinfo[-2],
                               '#(half) occupied flat bands':len(ocp_bandwidth),
                               'occbandwidth':ocp_bw_info,
                               '#occupied flat bands': len(uocp_bandwidth),
                               'uoccbandwidth':uocp_bw_info}
            index += 1

    data.to_excel('./flatband.xlsx')

    print('finished!')


if __name__ == '__main__':
    main()
