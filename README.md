# Pseudo-spectral code for stratified flows in 3D
In order to compile perform following steps
```
$ cd $PS3D_ROOT
$ ./bootstrap
$ mkdir build
$ cd build
$ ../configure --prefix=$PS3D_PREFIX
$ make
$ make install
```
where `$PS3D_ROOT` is the root directory of this repository and `$PS3D_PREFIX` is the location where you want to install. All executables will then be
located in `$PS3D_PREFIX/bin`.
