# Pseudo-spectral code for turbulent flows in 3D
PS3D is developed to solve Euler equations for an incompressible inviscid flow in a
horizontally-periodic domain and confined between parallel free-slip (including stress) surfaces.
We especially target the solution of Beltrami flows.

## How to compile
In order to compile perform following steps
```
$ cd $PS3D_ROOT
$ ./bootstrap
$ mkdir build
$ cd build
$ ../configure --prefix=$PS3D_PREFIX --enable-openmp --enable-verbose
$ make
$ make install
```
where `$PS3D_ROOT` is the root directory of this repository and `$PS3D_PREFIX` is the location where you want to install. All executables will then be located in `$PS3D_PREFIX/bin`.

## How to run
In the directory `examples` you can find configuration scripts. You first need to generate a field input file by
running the Beltrami program providing a *.nml file, e.g.
```
$ beltrami --config beltrami32x32x32.nml
```
This will generate a field file in NetCDF format. The name of the file is specified in the *.nml file which is
`beltrami_32x32x32.nc` in this example. After that you run PS3D by
```
$ ps3d --config beltrami_32.config
```
The NetCDF file is specified in the configuration file too.
