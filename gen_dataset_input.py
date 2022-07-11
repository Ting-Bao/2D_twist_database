# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.07
#
# this file generate input file of dataset, the calculation should be done in w001
# 

from utility.utils import *
from utility.get_data_batch import *

namelist='data/finalchoice_calculation/part1_prioirty/namelist.json'
#namelist='data/finalchoice_calculation/part2_normal/namelist.json'
outfile='data/finalchoice_calculation/part1_prioirty/'
#outfile='data/finalchoice_calculation/part2_normal/'


def main():
    name=read_namelist(namelist)
    for name in namelist:
        get_batch()

if __name__=='__main__':
    main()