# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2023.03
# this code is to plot the detailed mae fig

import json
import os
import shutil
import pandas as pd
from tqdm import tqdm
from utility.connect_interpolate_func import connect_band
from utility.refineplot_reconnected_band import *

datapath = 'data/infer_calculation/'


def get_path():
    matlist = os.listdir(datapath)
    matlist = [i for i in matlist if os.path.isdir(datapath+i)]
    pathlist = []
    for mat in matlist:
        twlist = os.listdir(os.path.join(datapath, mat))
        twlist = [i for i in twlist if os.path.isdir(os.path.join(datapath, mat, i))]
        for tw in twlist:
            path = os.path.join(datapath, mat, tw, 'band')
            pathlist.append(path)
    return pathlist
