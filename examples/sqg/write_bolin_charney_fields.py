#!/usr/bin/env python

"""
In order to run this script you need 
	https://github.com/EPIC-model/pytools
	(we used commit fdff946a15edeeef2d87a5dd47869cee2c925c67)
and the input data
	ba.r8 - buoyancy anomaly
	ox.r8 - x-vorticity component
	oy.r8 - y-vorticity component
	oz.r8 - z-vorticity component
This script expects the four files being placed in the
sub-directory "Bolin-Charney-balance". When successful, the
NetCDF input file "bolin_charney_512x512x256.nc" is created.
"""

import tools.netcdf as nc
import numpy as np

ncf = nc.FieldWriter()

ncf.open('bolin_charney_512x512x256.nc')

# number of cells
nx = 512
ny = 512
nz = 256

def load_dataset(ff):
    N = nx*ny*(nz+1)
    in_file = open('Bolin-Charney-balance/' + ff,'r')
    raw_array = np.fromfile(in_file,dtype = np.float64)
    in_file.close()
    bb = np.empty([ny,nx,nz+1])
    bb = raw_array[1:N+1].reshape(nz+1,nx,ny)

    qq = np.empty([nz+1,ny,nx])
    for iz in range(nz+1):
        qq[iz, :, :] = bb[iz, :, :].T
    return qq


# domain origin
origin = (-np.pi, -np.pi, -0.03125 * np.pi)

# domain extent
extent = (2.0*np.pi, 2.0*np.pi, 0.03125*np.pi)

files = ['ba.r8', 'ox.r8', 'oy.r8', 'oz.r8']
names = ['buoyancy_anomaly', 'x_vorticity', 'y_vorticity', 'z_vorticity']
units = ['m/s^2', '1/s', '1/s', '1/s']

for i in range(len(files)):
    print("Write", names[i], end=' ')
    data = load_dataset(files[i])
    print("min:", data.min(), "max:", data.max())
    ncf.add_dataset(names[i], data, unit=units[i], long_name=names[i])

#ncf.add_box(origin, extent, [nx, ny, nz])

# mesh spacings
dx = extent[0] / nx
dy = extent[1] / ny
dz = extent[2] / nz

x = np.linspace(origin[0], origin[0]+extent[0], nx, endpoint=False)
y = np.linspace(origin[1], origin[1]+extent[1], ny, endpoint=False)
z = np.linspace(origin[2], origin[2]+extent[2], nz+1, endpoint=True)
ncf.add_axis('x', x)
ncf.add_axis('y', y)
ncf.add_axis('z', z)

ncf.add_axis('t', [0.0])


ncf.add_physical_quantity('l_planetary_vorticity', 'true')
ncf.add_physical_quantity('planetary_angular_velocity', 0.5)
ncf.add_physical_quantity('latitude_degrees', 90.0)
# Note: "rossby_number" is currently only used by "init_sqg.f90"
# and not part of the simulation code.
ncf.add_physical_quantity('rossby_number', 0.125)
ncf.add_physical_quantity('squared_buoyancy_frequency', 64.0)

ncf.add_parameter("ncells", [nx, ny, nz])
ncf.add_parameter("extent", extent)
ncf.add_parameter("origin", origin)
ncf.add_parameter("grid_type", "uniform")

ncf.close()
