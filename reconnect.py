# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.07

import json
import os

datapath = 'data/infer_calculation/'


def main():
    matlist = os.listdir(datapath)
    matlist = [i for i in matlist if os.path.isdir(datapath+i)]
    pathlist=[]
    for mat in matlist:
        twlist = os.listdir(os.path.join(datapath, mat))
        for tw in twlist:
            path = os.path.join(datapath, mat, tw)
            pathlist.append(path)

if __name__ == '__main__':
    main()
