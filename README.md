[![License](https://img.shields.io/github/license/matt-frey/ps3d)](https://github.com/matt-frey/ps3d/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6979867.svg)](https://doi.org/10.5281/zenodo.6979867)

# Pseudo-spectral code for turbulent flows in 3D
PS3D is developed to solve the Euler equations for an incompressible inviscid flow in a
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
The NetCDF file is specified in the configuration file too. In order to restart a simulation from a previous
simulation, you need to update `field_file` and add the step number input with `field_step` to the configuration
file. Please see the file `beltrami_256_restart.config` in the `examples` directory for an example.

## How to post-process
In the directory `plotting` we provide many pre-configured scripts to generate cross sections and 3D views of
fields. The class `nc_reader` is used to retrieve data from the NetCDF file. The visualisation is done with
matplotlib (version 3.5.2) and ParaView (version 5.10.1). You can install the Python environment with the provided
`requirements.txt` file found in the root directory of this repository. The basic command to accomplish this is
```
$ conda create --name <env> --file requirements.txt
```
where `<env>` denotes a user-defined environment name like `ps3d-env`. To install conda, please read the relevant
documentation.

The enstrophy production rates are calculated using the program `production`. Use the flag `--help` for further information. Note: The post-processing codes are currently switched off as we moved to distributed memory parallelism and these codes are not yet verified.


## How to cite
The source code can be cited using the latest Zenodo DOI given on the top of this readme. A corresponding paper is published and can be cited using
the BibteX below. You can also find further information at https://doi.org/10.1017/jfm.2022.1007

```
@article{dritschel_frey_2023,
  title={The stability of inviscid Beltrami flow between parallel free-slip impermeable boundaries},
  volume={954},
  DOI={10.1017/jfm.2022.1007},
  journal={Journal of Fluid Mechanics},
  publisher={Cambridge University Press},
  author={Dritschel, David and Frey, Matthias},
  year={2023},
  pages={A31}
}
```
