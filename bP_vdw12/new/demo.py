"""Twist2D Demo."""
#%%
# +------------+
# | Usage Demo |
# +------------+
from twist2d import *

# Create an object for t2d
twist_demo = Twist2D()

# Settings for twist supercell generation
lattice_type = 'rec'
# 60 triangular = 'tri'   //   120 hexagonal = 'hex'
# square = 'squ'   //   rectangular = 'rec'

# set superlattice index m & n
m = 19
n = 10
# set interlayer distance

# rules for superlattice vector settings
# for 60 triangular lattice, layer_1: (m,n) (-n,m+n), layer_2: (n,m) (-m,n+m)
# for 90 square lattice, the same to 60 triangular case
# for 90 rectangular lattice, layer_1: (m,0) (0,n), layer_2: (m,-1) (1,n)
# for 120 hexagonal lattice, layer_1: (m,n) (-n,m-n), layer_2: (n,m) (-m,n-m)

# Initialize the different twisted layers
#  - super_a1_mult,  super_a2_mult: supercell vector a1',a2' based on a1,a2
#  - layer_dis: the layer distance of this layer to next layer, default 2A.
#  - scs_x, scs_y: supercell shift in x,y direction in angstroms, default 0A.
#  - prim_poscar: POSCAR for primitive cell of current layer, default 'POSCAR'. 

# if lattice_type == 'rec':
#--> 1st layer
super_a1_mult = [m, 0]
super_a2_mult = [0, n]
twist_demo.add_layer(super_a1_mult, super_a2_mult, layer_dis=2.981189843685226057652085314149, scs_x=1.6484998174736949, prim_poscar="POSCAR")
#--> 2nd layer
super_a1_mult = [m, -1]
super_a2_mult = [1, n]
twist_demo.add_layer(super_a1_mult, super_a2_mult, prim_poscar="POSCAR")
# #--> 3rd layer
# super_a1_mult = [n, m]
# super_a2_mult = [-m, n+m]
# twist_demo.add_layer(super_a1_mult, super_a2_mult, prim_poscar="POSCAR-BN")

# Twisting the layers
#  - start_z: The lowest atom's fractional coordinates in z, default 0.1
#  - super_a3_z: The length of the c vector in z direction, default 20A.
twist_demo.twist_layers(start_z=0.1, super_a3_z=20)

# Write results to the file
twist_demo.write_res_to_poscar()

# (Optional) Calculate the twisted angles of each layer in degree 
twisted_angles = twist_demo.calc_layers_twist_angles()
print(m, n, twisted_angles)

# PROGRAM END

#%%
# +-------------------+
# | Special condition |
# +-------------------+
#from twist2d import *
#
## If you are twisting a bilayer graphene-like system, 
##   you can write more simply like this:
#
## Twist bilayer graphene-like structures
#tbg_demo = TwistBGL()
#tbg_demo.gen_TBGL(6, 7)
##tbg_demo.gen_TBG(m=6, n=7, prim_poscar='POSCAR', poscar_out="POSCAR.T2D.vasp", start_z=0.1, super_a3_z=20.0, layer_dis=2.0, scs_x=0.0, scs_y=0.0)
#
## (Optional) Calculate the twisted angles of each layer in degree 
#twisted_angles = tbg_demo.calc_layers_twist_angles()
#print(twisted_angles)
#
##PROGRAM END

# %%
