#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Catalogue functionality to aid LAE population analysis.

Examples:
--------
    Import the module:
        >>> import catalogue as cat

    Create a catalogue from scratch:
        >>> mock = cat.Catalogue(False, z, n, coords, mass_h, mag_uv,
        ...                      lum_lya, rew_lya, t_igm)

    Save a catalogue to HDF5:
        >>> mock.write_to_HDF5(filename)

    Load a catalogue from an HDF5 file:
        >>> mock = cat.Catalogue(True, filename)

    Find coordinates of LAEs within spatial slice:
        >>> x, y = mock.spatial_slice(zdepth, depth, axis='z')

"""

import h5py as h5
import numpy as np


class Catalogue:
    """The catalogue object describing the LAE population.

    Attributes:
    -----------
        z (float): mean redshift of LAE population
        n (int): number of LAEs in population
        coords (array_like, shape (3,n)): (x,y,z) coordinates in cMpc/h
        mass_h (array_like, shape (n)): host halo masses in Msun/h
        mag_uv (array_like, shape (n)): UV magnitude in AB system
        lum_lya (array_like, shape (n)): intrinsic Lya luminosity in erg/s
        rew_lya (array_like, shape (n)): intrinsic Lya REW in Angstroms
        t_igm (array_like, shape (n)): IGM transmission fraction

    """

    def __init__(self, read_from_file, *args):
        """Initialise catalogue object.

        Args:
        -----
            read_from_file (bool): where to load from file or arguments.
            *args: if read_from_file is True then a single str with the
                path to the HDF5 path, otherwise the arguments to pass
                to read_from_args().

        """

        if read_from_file:
            self.read_from_hdf5(args[0])
        else:
            self.read_from_args(args)

    def read_from_hdf5(self, fname):
        """Initialise catalogue from HDF5 file.

        Args:
        -----
            fname (str): path to HDF5 file

        """

        with h5.File(fname, 'r') as f:
            self.z = f['/Catalogue'].attrs['redshift']
            self.n = f['/Catalogue'].attrs['n']
            self.coords = f['/Catalogue/coords'][()]
            self.mass_h = f['/Catalogue/mass_h'][()]
            self.mag_uv = f['/Catalogue/mag_uv'][()]
            self.lum_lya = f['/Catalogue/lum_lya'][()]
            self.rew_lya = f['/Catalogue/rew_lya'][()]
            self.t_igm = f['/Catalogue/t_igm'][()]

    def read_from_args(self, *args):
        """Initialise catalogue from provided arguments.

        Args:
        -----
            z (float): mean redshift of LAE population
            n (int): number of LAEs in population
            coords (array_like, shape (3,n)): (x,y,z) coordinates in cMpc/h
            mass_h (array_like, shape (n)): host halo masses in Msun/h
            mag_uv (array_like, shape (n)): UV magnitude in AB system
            lum_lya (array_like, shape (n)): intrinsic Lya luminosity in erg/s
            rew_lya (array_like, shape (n)): intrinsic Lya REW in Angstroms
            t_igm (array_like, shape (n)): IGM transmission fraction

        """
        assert len(args) == 8

        self.z = args[0]
        self.n = args[1]

        assert args[2].shape == (3, self.n)
        assert args[3].size == self.n
        assert args[4].size == self.n
        assert args[5].size == self.n
        assert args[6].size == self.n
        assert args[7].size == self.n

        self.coords = args[2]
        self.mass_h = args[3]
        self.mag_uv = args[4]
        self.lum_lya = args[5]
        self.rew_lya = args[6]
        self.t_igm = args[7]

    def print_info(self):
        """Print useful info.

        """
        print("\nCatalogue info:")
        print("z = {:.3f}".format(self.z))
        print("n = {}".format(self.n))
        print("Halo mass [Msun/h]: min = {:.2e}, max = {:.2e}".format(
            self.mass_h.min(), self.mass_h.max()))
        print("UV magnitude [AB]: min = {:.1f}, max = {:.1f}".format(
            self.mag_uv.max(), self.mag_uv.min()))
        print("Lya luminosity [erg/s]: min = {:.2e}, max = {:.2e}".format(
            self.lum_lya.min(), self.lum_lya.max()))
        print("Lya REW [Angstroms]: min = {:.1f}, max = {:.1f}".format(
            self.rew_lya.min(), self.rew_lya.max()))
        print("IGM transmission: min = {:.3f}, max = {:.3f}".format(
            self.t_igm.min(), self.t_igm.max()))

    def write_to_HDF5(self, fname, overwrite=False):
        """Write catalogue to HDF5 file.

        Args:
        -----
            fname (str): path to output file.
            overwrite (bool): if True, overwrite existing files.

        """

        print("Writing to:", fname)
        if overwrite:
            wcode = 'w'
        else:
            wcode = 'w-'

        with h5.File(fname, wcode) as f:
            h5catalogue = f.create_group("Catalogue")

            h5catalogue.attrs['z'] = self.z
            h5catalogue.attrs['n'] = self.n
            h5catalogue.create_dataset("coords", data=self.coords)
            h5catalogue.create_dataset("mass_h", data=self.mass_h)
            h5catalogue.create_dataset("mag_uv", data=self.mag_uv)
            h5catalogue.create_dataset("lum_lya", data=self.lum_lya)
            h5catalogue.create_dataset("rew_lya", data=self.rew_lya)
            h5catalogue.create_dataset("t_igm", data=self.t_igm)

    def write_to_binary(self, fname, overwrite=False):
        """Write catalogue to custom format binary file.

        Args:
        -----
            fname (str): path to output file.
            overwrite (bool): if True, overwrite existing files.

        """

        print("Writing to:", fname)
        if overwrite:
            wcode = 'w'
        else:
            wcode = 'x'

        z = np.array(self.z).astype(np.float64)
        n = np.array(self.n).astype(np.int32)
        xcoords = self.coords[0, :].astype(np.float64)
        ycoords = self.coords[1, :].astype(np.float64)
        zcoords = self.coords[2, :].astype(np.float64)
        mass_h = self.mass_h.astype(np.float64)
        mag_uv = self.mag_uv.astype(np.float64)
        lum_lya = self.lum_lya.astype(np.float64)
        rew_lya = self.rew_lya.astype(np.float64)
        t_igm = self.t_igm.astype(np.float64)

        with open(fname, wcode) as f:
            z.tofile(f)
            n.tofile(f)
            xcoords.tofile(f)
            ycoords.tofile(f)
            zcoords.tofile(f)
            mass_h.tofile(f)
            mag_uv.tofile(f)
            lum_lya.tofile(f)
            rew_lya.tofile(f)
            t_igm.tofile(f)

    def spatial_slice(self, zdepth, depth, axis='z'):
        """Create a configuration-space slice of the catalogue volume.

        Args:
        -----
            zdepth (float): the leading edge of the slice in cMpc/h
            depth (float): the depth of the slice in cMpc/h
            axis (str): the coordinate along which to slice

        Returns:
        --------
            x, y (array_like): coordinates of LAEs in the slice.
        """

        if axis == 'z':
            x = self.coords[0, :]
            y = self.coords[1, :]
            z = self.coords[2, :]
        elif axis == 'y':
            x = self.coords[2, :]
            y = self.coords[0, :]
            z = self.coords[1, :]
        elif axis == 'x':
            x = self.coords[1, :]
            y = self.coords[2, :]
            z = self.coords[0, :]
        else:
            raise ValueError("Must slice along one of the coordinate axes!")

        zdiff = np.absolute(z - zdepth)
        mask = np.logical_and(zdiff >= 0.0, zdiff < depth)
        slice_x = x[mask]
        slice_y = y[mask]
        print("Number of LAEs in slice:", slice_x.size)

        return slice_x, slice_y
