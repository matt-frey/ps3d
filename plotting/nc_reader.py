import netCDF4 as nc
import os
import numpy as np


class nc_reader:
    def __init__(self):
        self._ncfile = None

        self._derived_fields = [
            'vorticity_magnitude',
            'helicity'
        ]

    def open(self, fname):
        if not os.path.exists(fname):
            raise IOError("File '" + fname + "' does not exist.")
        self._ncfile = nc.Dataset(fname, "r", format="NETCDF4")

    def close(self):
        self._ncfile.close()

    def get_num_steps(self):
        return self._ncfile.dimensions['t'].size

    def get_box_extent(self):
        return self.get_global_attribute("extent")

    def get_box_ncells(self):
        return self.get_global_attribute("ncells")

    def get_box_origin(self):
        return self.get_global_attribute("origin")

    def get_axis(self, name):
        axis = self.get_all(name)
        if name == 'x' or name == 'y':
            # copy periodic grid point
            axis = np.append(axis, abs(axis[0]))
        return axis

    def get_meshgrid(self):
        x = self.get_axis('x')
        y = self.get_axis('y')
        z = self.get_axis('z')

        xg, yg, zg = np.meshgrid(x, y, z, indexing='ij')

        # 13 July 2022
        # https://stackoverflow.com/questions/1827489/numpy-meshgrid-in-3d
        assert np.all(xg[:, 0, 0] == x)
        assert np.all(yg[0, :, 0] == y)
        assert np.all(zg[0, 0, :] == z)

        return xg, yg, zg


    def get_all(self, name):
        if not name in self._ncfile.variables.keys():
            raise IOError("Dataset '" + name + "' unknown.")
        return np.array(self._ncfile.variables[name])

    # returns a dataset with axis ordering (x, y, z)
    def get_dataset(self, step, name):

        if name in self._derived_fields:
            return self._get_derived_dataset(step, name)

        if not name in self._ncfile.variables.keys():
            raise IOError("Dataset '" + name + "' unknown.")

        nsteps = self.get_num_steps()
        if step > nsteps - 1:
            raise ValueError("Dataset has only steps 0 to " + str(nsteps - 1) + ".")

        fdata = np.array(self._ncfile.variables[name][step, ...])
        fdata = self._copy_periodic_layers(fdata)

        # change ordering from (z, y, x) to (x, y, z)
        fdata = np.transpose(fdata, axes=[2, 1, 0])
        return fdata

    def _get_derived_dataset(self, step, name):
        if name == 'vorticity_magnitude':
            x_vor = self.get_dataset(step=step, name='x_vorticity')
            y_vor = self.get_dataset(step=step, name='y_vorticity')
            z_vor = self.get_dataset(step=step, name='z_vorticity')
            return np.sqrt(x_vor ** 2 + y_vor ** 2 + z_vor ** 2)
        if name == 'helicity':
            u = self.get_dataset(step=step, name='x_velocity')
            v = self.get_dataset(step=step, name='y_velocity')
            w = self.get_dataset(step=step, name='z_velocity')
            xi = self.get_dataset(step=step, name='x_vorticity')
            eta = self.get_dataset(step=step, name='y_vorticity')
            zeta = self.get_dataset(step=step, name='z_vorticity')
            return u * xi + v * eta + w * zeta

    def get_dataset_attribute(self, name, attr):
        if not name in self._ncfile.variables.keys():
            raise IOError("Dataset '" + name + "' unknown.")

        if not attr in self._ncfile.variables[name].ncattrs():
            raise IOError("Dataset attribute '" + name + "' unknown.")

        return self._ncfile.variables[name].getncattr(attr)

    def get_dataset_min_max(self, name, indices=None):
        nsteps = self.get_num_steps()
        data = self.get_dataset(0, name, indices=indices)
        vmax = data.max()
        vmin = data.min()
        for step in range(1, nsteps):
            data = self.get_dataset(step, name, indices=indices)
            vmax = max(vmax, data.max())
            vmin = min(vmin, data.min())
        return vmin, vmax

    def get_global_attribute_names(self):
        return list(self._ncfile.ncattrs())

    def get_global_attribute(self, name):
        if not name in self._ncfile.ncattrs():
            raise IOError("Global attribute '" + name + "' unknown.")
        attr = self._ncfile.getncattr(name)
        if isinstance(attr, np.bytes_):
            attr = attr.decode()
        return attr

    def get_diagnostic(self, name):
        if not name in self._ncfile.variables.keys():
            raise IOError("Dataset '" + name + "' unknown.")
        return np.array(self._ncfile.variables[name])

    def _get_step_string(self, step):
        return str(step).zfill(10)

    # 18 Feb 2022
    # https://stackoverflow.com/questions/8450472/how-to-print-a-string-at-a-fixed-width
    # 19 Feb 2022
    # https://stackoverflow.com/questions/873327/pythons-most-efficient-way-to-choose-longest-string-in-list
    def __str__(self):
        print("=" * 80)
        # print global attributes
        print("GLOBAL ATTRIBUTES:")
        l = len(max(self._ncfile.ncattrs(), key=len))
        fmt = '{0: <' + str(l) + '}'
        for key in self._ncfile.ncattrs():
            print(fmt.format(key), "\t", self._ncfile.getncattr(key))
        print("-" * 80)

        print("DIMENSIONS:")

        for dim in self._ncfile.dimensions:
            print("    ", dim, "=", self._ncfile.dimensions[dim].size)

        print("-" * 80)

        print("VARIABLES:")
        # get first variable name
        name = list(self._ncfile.variables.keys())[0]

        # get length of longest attribute string
        l = len(max(self._ncfile.variables[name].ncattrs(), key=len))
        fmt = '{0: <' + str(l) + '}'

        # print variables and their attributes
        for var in self._ncfile.variables:
            print("    ", var)
            for attr in self._ncfile.variables[var].ncattrs():
                print("\t", fmt.format(attr), "\t", self._ncfile.variables[var].getncattr(attr))
        print("=" * 80)
        return ""

    def _copy_periodic_layers(self, field):
        nz, ny, nx = field.shape
        field_copy = np.empty((nz, ny+1, nx+1))
        field_copy[:, 0:ny, 0:nx] = field.copy()
        field_copy[:, ny, :] = field_copy[:, 0, :]
        field_copy[:, :, nx] = field_copy[:, :, 0]
        return field_copy
