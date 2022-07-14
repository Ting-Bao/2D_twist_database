# from xurz

import math
from gen_opmx_input_func import *

def cutoff_radius(element_name, accuracy):
    r_born_to_ang = 0.52917721  # convert Bohr to Ang
    element_info = find_basis(element_name, accuracy)
    basis = element_info[3]
    element_and_radius = basis.split('-')[0]
    radius_tmp = filter(lambda ch: ch in '0123456789.', element_and_radius)

    radius_str = ''
    for i in radius_tmp:
        radius_str += i
    radius_float = float(radius_str)

    radius_ang = radius_float * r_born_to_ang * 2
    radius_out = float(math.ceil(2*radius_ang)/2) # this will have a step of 0.5 --Ting BAO

    return radius_out  # in unit of Ang

#print(cutoff_radius('Nb',3))