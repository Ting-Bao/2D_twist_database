import numpy as np
from scipy.interpolate import interp1d

# Example data
data = [1, 3, 5, 7, 9, 11, 13, 15, 17]

# Define the x-axis values (assumed to be evenly spaced)
x = np.arange(len(data))

# Define a function for the interpolation
interp_func = interp1d(x, data, kind='cubic')

# Choose the number of points for the smoothed data
num_points = 2*len(data)

# Compute the smoothed data using interpolation
smooth_data = interp_func(np.linspace(x.min(), x.max(), num_points))

# Print original and smooth data
print("Original data:", data)
print("Smooth data:", smooth_data)
