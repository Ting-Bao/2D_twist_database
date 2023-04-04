# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.07
# this code is to reconect and calculate and do statics on band width under VBM
import json
import os

datapath = 'data/infer_calculation/'
remote_path = 'xurz@w001:/home/xurz/temp_baot/twisted_band_opmx/data/twisted_band_opmx/'


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
            dft_path = os.path.join(datapath, mat, tw,'dft_band')
            os.makedirs(dft_path, exist_ok = True)
            os.system('scp {}{}_{}/openmx.Band {}/'.format(remote_path,mat,tw,dft_path))

if __name__ == '__main__':
    main()

