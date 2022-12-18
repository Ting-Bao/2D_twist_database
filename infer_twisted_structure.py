# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.12
#
# this file is condcuted after calc_twisted_structure.py, 
#   which is used for compare the MAE of Hamiltonian matrix element and band structrue
#
# This file should be run on 403ubuntu

from utility.utils import *
import os
import shutil
import json

results_collection='/home/tingbao/Desktop/TS216NAS/Public/baot/results_collection/'
model_set='/home/tingbao/Desktop/TS216NAS/Public/baot/modelset/'

opmx_path='./data/twisted_band_opmx/'
model_path = model_set
infer_path = './data/infered_results/'

def select_model(source,target):
    '''
    find human-collect models and put into a uniform format
    '''
    models=[ i for i in os.listdir(source) if os.path.isdir(source+i) ]
    #print(len(model))
    for i in models:
        print(i,len(os.listdir(source+i+'/result/')))
    print("{} Models found!",format(len(models)))



if __name__=='__main__':
    select_model(source=results_collection, target= model_set)

