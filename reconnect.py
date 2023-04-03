# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.07
# this code is to reconect and calculate and do statics on band width under VBM
import json
import os
from utility.connect_interpolate_func import connect_band
from utility.refineplot_reconnected_band import *

datapath = 'data/infer_calculation/'



def main():
    matlist = os.listdir(datapath)
    matlist = [i for i in matlist if os.path.isdir(datapath+i)]
    pathlist=[]
    for mat in matlist:
        twlist = os.listdir(os.path.join(datapath, mat))
        twlist = [i for i in twlist if os.path.isdir(os.path.join(datapath, mat,i))]
        for tw in twlist:
            path = os.path.join(datapath, mat, tw,'band')
            pathlist.append(path)
    
    for path in ['data/infer_calculation/Bi2Te3-61e21cf8e36d/2-1/band/']: #pathlist:#
        connect_band(inputpath = path, outputpath = path, spline_interpolate=False) 
        # write band_reconnected.json, the fermi energy need further refinement (implemented above)
        band_obj = load_process_band(readpath=path, savepath=path)
        band_obj.get_reconnected_bandwidth()
    print('finished!')
        
    

if __name__ == '__main__':
    main()

