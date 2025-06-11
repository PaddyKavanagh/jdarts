#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @authors: Patrick Kavanagh
#
"""
This module contains functions to generate quicklook plots for JWST MIRI MRS data.

It includes functions for plotting 1D spectra from text files, 
2D median images of MIRI data cubes, where the median is computed along the spectral axis, 
and velocity channel maps from a MIRI data cube, where the user specifies the line rest wavelength and velocity range.

"""
import os
import gc
import sys
import glob
import shutil
import numpy as np
import math
import warnings
#
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from regions import Regions
from astropy.visualization import ImageNormalize, PercentileInterval, LogStretch, SqrtStretch
from astropy import constants as cc

def mrs_quicklook_fullband_spectra(spectra, plot_scale='linear', save_plot=False,
                           output_filename='soectrum_plot.png'):
    """
    Generates a quicklook plot of one or more full-band spectra from text files.
    Each text file is expected to contain three columns: wavelength (in microns),
    flux (in Jy), and surface brightness.

    Parameters
    ----------
    spectra : list of str
        A list of paths to the .txt files containing spectra.
    plot_scale : str, optional
        The scaling of the plot. If set to 'linear', a linear scale is used for both axes.
        Any other value will result in a log-log scale. Defaults to 'linear'.
    save_plot : bool, optional
        If True, the plot is saved to a PNG file. Defaults to False.
    output_filename : str, optional
        The name of the PNG file to save the plot to, if save_plot is True.
        Defaults to 'spectrum_plot.png'.

    Returns
    -------
    None
    """
    # setup plot
    fig, axs = plt.subplots(1, 1, figsize=(10, 5))

    for n, f in enumerate(spectra):

        wavelength, flux, surf_bright = np.loadtxt(f, unpack=True)
        axs.step(wavelength, ap_flux, linewidth=1.5, markersize=0, alpha=0.9)

    # set other things
    axs.set_xlabel('Wavelength (micron)', fontsize=16)
    axs.set_ylabel('Flux (Jy)', fontsize=16)
    axs.set_xlim(4.8,28)

    if not plot_scale == 'linear':
        axs.set_yscale('log')
        axs.set_xscale('log')
        plt.xticks([5, 6, 7, 8, 9, 10, 12, 15, 20, 25], [5, 6, 7, 8, 9, 10, 12, 15, 20, 25])

    plt.tight_layout(h_pad=0)
    plt.show()

    if save_plot:
        fig.savefig(output_filename, format='png', dpi=300)


def mrs_quicklook_cubes(cubes, plot_scale='linear', save_plot=False,
                           output_filename='cubes_plot.png'):
    """
    Generates a quicklook plot of median images computed from a set of MIRI data cubes.
    For each cube, it calculates the median along the spectral axis and displays
    the resulting 2D spatial image with a 1 arcsecond scale bar.

    Parameters
    ----------
    cubes : list of str
        A list of paths to the MIRI datacube FITS files.
    plot_scale : str, optional
        The scaling of the plot. Defaults to 'linear'.
    save_plot : bool, optional
        If True, the plot is saved to a PNG file. Defaults to False.
    output_filename : str, optional
        The name of the PNG file to save the plot to, if save_plot is True.
        Defaults to 'cubes_plot.png'.

    Returns
    -------
    None
    """
    # setup plot
    fig, axs = plt.subplots(2, 6, figsize=(12, 4), facecolor='black')
    axs = axs.ravel()

    # iterate through cubes
    for n, c in enumerate(cubes):
        with fits.open(c) as hdu:
            data = hdu[1].data      # data arr
            scale_line = (1 / 3600) / hdu[1].header['CDELT1']       # 1 arcsec scale
            name = hdu[0].header['CHANNEL'] + hdu[0].header['BAND']     # name of band

        # use the median of the cube and setup plot normalisation
        sci_arr = np.nanmedian(data, axis=0)
        norm = ImageNormalize(sci_arr, interval=PercentileInterval(99.5), stretch=LogStretch(a=1), vmin=0)

        # plot
        axs[n].imshow(sci_arr, cmap='inferno', interpolation='nearest', origin='lower', norm=norm, alpha=1.0)
        axs[n].annotate(name, (0.05, 0.87), color='white', xycoords='axes fraction', fontweight='bold',
                        fontsize=16, rotation=0.0)
        axs[n].set_facecolor('black')

        # plot the scale line
        sl_x = data.shape[2] * 0.05
        sl_y = data.shape[1] * 0.15
        axs[n].plot([sl_x, sl_x + scale_line], [sl_y, sl_y], markersize=0, linewidth=2, color='white')

        # switch off the axes and thicken the spines
        axs[n].axis('off')
        for axis in ['top', 'bottom', 'left', 'right']:
            axs[n].spines[axis].set_linewidth(1.5)

    plt.tight_layout(pad=0.0, h_pad=0, w_pad=0)
    plt.show()

    if save_plot:
        fig.savefig(output_filename, format='png', dpi=300)


def generate_plot_grid(n):
    """
    Generates a square grid for plotting based on the number of frames.

    Parameters
    ----------
    n : int
        Number of frames.

    Returns
    -------
    tuple of int
        A tuple (rows, cols) representing the dimensions of the grid for plotting.
    
    """
    a = round(math.sqrt(n))
    b = a + 1

    if a ** 2 > n:
        return (a, a)
    else:
        return (a, b)


def find_nearest(array, value):
    """
    Finds the index of the nearest value in an array.

    Parameters
    ----------
    array : numpy.ndarray
        The input NumPy array.
    value : float or int
        The value to find the nearest element to.

    Returns
    -------
    int
        The index of the element in the array closest to value.
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def wav_to_vel(wavelength, line_wavelength):
    """
    Convert wavelength to velocity using the Doppler formula.

    Parameters
    ----------
    wavelength : astropy.units.Quantity
        Wavelength in microns.
    line_wavelength : astropy.units.Quantity
        Line rest wavelength in microns.
    Returns
    -------
    astropy.units.Quantity
        Radial velocity in km/s rounded to nearest integer.
    """
    radial_vel = ((wavelength - line_wavelength) / line_wavelength) * cc.c
    radial_vel = np.round(radial_vel.to(u.km / u.s), 0)

    return radial_vel


def plot_channel_maps(cube, line_name, line_rest_wavelength, velocity_range, radial_vel=0,
                      plot_regions=False, region_file=None, save_plot=False, plot_name='channel_maps.png'):
    """
    Plots a velocity channel map from a MIRI MRS data cube around a specified rest wavelength and velocity range. 
    Optionally plots a DS9 region on one of the maps.

    Parameters
    ----------
    cube : str
        The path to the MIRI data cube FITS file.
    line_name : str
        The name of the spectral line being studied (e.g., '[Ne II]').
    line_rest_wavelength : float
        The rest wavelength of the spectral line in microns.
    velocity_range : list or tuple of float
        A list or tuple [min_velocity, max_velocity] in km/s defining the range
        of velocities to plot channel maps for.
    radial_vel : float, optional
        The systemic radial velocity of the object in km/s. Defaults to 0.
    plot_regions : bool, optional
        If True and a region_file is provided, a DS9 region is overlaid on one
        of the channel maps. Defaults to False.
    region_file : str, optional
        The path to the DS9 region file. Required if plot_regions is True.
        Defaults to None.
    save_plot : bool, optional
        If True, the plot is saved to a PNG file. Defaults to False.
    plot_name : str, optional
        The name of the PNG file to save the plot to, if save_plot is True.
        Defaults to 'channel_maps.png'.

    Returns
    -------
    None

    Raises
    ------
    IOError
        If plot_regions is True but no region_file is supplied.
    """
    # if plot_regions is True, check a region file supplied
    regions = None
    pix_region = None

    if plot_regions:
        if region_file is None:
            raise IOError('Plot regions is True but no region file supplied. Aborting')
        else:
            regions = Regions.read(region_file, format='ds9')
            sky_region = regions[0]
            pix_region = sky_region.to_pixel(wcs)

    # read info from cube file
    with fits.open(cube) as hdu:
        data = hdu[1].data
        wcs = WCS(hdu[1].header)
        wcs = wcs.dropaxis(2)
        wavelength = (np.arange(hdu[1].header['NAXIS3']) * hdu[1].header['CDELT3'] +
                      hdu[1].header['CRVAL3']) * u.micron

    # update this
    frame_binning = 2

    # sort out some things
    line_rest_wavelength *= u.micron
    radial_vel *= u.km / u.s
    line_wavelength = line_rest_wavelength * (cc.c / (radial_vel + cc.c))
    line_ind = find_nearest(wavelength.value, line_wavelength.value)

    wav_arr_vel = wav_to_vel(wavelength, line_wavelength)
    wav_arr_vel_diff = np.diff(wav_arr_vel)

    wav_arr_vel_edges = wav_arr_vel - np.pad(wav_arr_vel_diff, (0, 1), 'constant')
    velocity_range_argmin = find_nearest(wav_arr_vel.value, velocity_range[0])
    velocity_range_argmax = find_nearest(wav_arr_vel.value, velocity_range[1])

    # generate the plot grid
    fnum = (velocity_range_argmax - velocity_range_argmin) // frame_binning
    plot_grid = generate_plot_grid(fnum - 1)

    # determine scale line in pixels
    scale_line = (1 / 3600) / hdu[1].header['CDELT1']

    # plot the apertures on the cubes
    fig, axs = plt.subplots(plot_grid[1], plot_grid[0], figsize=(plot_grid[0] * 2.1, plot_grid[1] * 2.1),
                            sharex=True, sharey=True)
    axs = axs.ravel()

    # need to select the frame for scaling, and cuts in spatial dimensions. Interval,
    # stretch and limits should also be set
    scale_frame = 10
    cuts = [0,-1,0,-1]  # y1,y2,x1,x2 - provide user with an option in future

    if pix_region is not None:
        pix_region.vertices.x = pix_region.vertices.x - cuts[2]
        pix_region.vertices.y = pix_region.vertices.y - cuts[0]

    # provide user control in future
    mrs_sci_arr = data[scale_frame, cuts[0]:cuts[1], cuts[2]:cuts[3]]
    norm = ImageNormalize(mrs_sci_arr, interval=PercentileInterval(95), stretch=SqrtStretch(), vmin=20, vmax=500)

    on_off = []
    if frame_binning == 2:
        for n in range(fnum):
            ll = n * 2
            frames = [velocity_range_argmin + ll, velocity_range_argmin + ll + 1]
            mrs_sci_arr = np.sum(data[frames[0]:frames[1], cuts[0]:cuts[1], cuts[2]:cuts[3]], axis=0)
            mrs_sci_arr = np.nan_to_num(mrs_sci_arr)
            radial_vel1 = wav_arr_vel_edges[frames[0]]
            radial_vel2 = wav_arr_vel_edges[frames[1] + 1]

            axs[n].imshow(mrs_sci_arr, cmap='inferno', interpolation='nearest', origin='lower', norm=norm, alpha=1.0)
            axs[n].annotate(line_name + ' ' + str(int(radial_vel1.value)) + r'$\rightarrow$' + str(
                int(radial_vel2.value)) + r' km s$^{-1}$', (0.05, 0.9), color='white',
                            xycoords='axes fraction', fontsize=10, rotation=0.0, fontweight='bold')
            axs[n].set_facecolor('black')

            axs[n].axis('off')
            for axis in ['top', 'bottom', 'left', 'right']:
                axs[n].spines[axis].set_linewidth(1.5)

                axs[n].axis('off')
                for axis in ['top', 'bottom', 'left', 'right']:
                    axs[n].spines[axis].set_linewidth(1.5)

                sl_x = 14
                sl_y = 1

    axs[0].annotate(r'1$^{\prime\prime}$', (0.75, 0.1), color='white',
                    xycoords='axes fraction', fontsize=12, rotation=0.0)

    axs[0].plot([sl_x, sl_x + scale_line], [sl_y, sl_y], markersize=0, linewidth=2, color='white')

    while n < (plot_grid[1] * plot_grid[0]) - 1:
        n += 1
        fig.delaxes(axs[n])

    if plot_regions and len(axs) > 5 and pix_region is not None:
        pix_region.plot(ax=axs[5], color='white', linewidth=2)

    plt.tight_layout(pad=0.0, h_pad=-9, w_pad=-3)
    plt.show()

    if save_plot:
        fig.savefig(plot_name, format='png', dpi=300)