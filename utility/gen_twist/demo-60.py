"""Twist2D Demo."""
#%%
# +------------+
# | Usage Demo |
# +------------+
from twist2d import *

# Create an object for t2d
twist_demo = Twist2D()

# Initialize the different twisted layers 
# (for different primitive cell condition)
### For the primitive cell of VSe2 in 120 degree
# m = 2
# n = 3
# The key parameters for twist-2d is the mult vectors
# --> 1st layer
# super_a1_mult = [m, n]
# super_a2_mult = [-n, m-n]
# twist_demo.add_layer(super_a1_mult, super_a2_mult, layer_dis=4.004972298593504, prim_poscar="POSCAR")
# --> 2nd layer
# super_a1_mult = [n, m]
# super_a2_mult = [-m, n-m]
# twist_demo.add_layer(super_a1_mult, super_a2_mult, prim_poscar="POSCAR")

### For the primitive cell of VSe2 in 60 degree (more common) change!!!
m = 1
n = 2
# --> 1st layer
super_a1_mult = [m, n]
super_a2_mult = [-n, m+n]
twist_demo.add_layer(super_a1_mult, super_a2_mult, layer_dis=3.01, prim_poscar="POSCAR1_60")  # change!!!
# --> 2nd layer
super_a1_mult = [n, m]
super_a2_mult = [-m, n+m]
twist_demo.add_layer(super_a1_mult, super_a2_mult, prim_poscar="POSCAR2_60")  # change!!!

# Fill the cell with the layers
twist_demo.twist_layers(start_z=0.3, super_a3_z=25)  #start_z用分数坐标 change!!!

# (Optional) Calculate the twisted angles of each layer in degree 
twisted_angles = twist_demo.calc_layers_twist_angles()
print(twisted_angles)

# Write results to the file
twist_demo.write_res_to_poscar()
