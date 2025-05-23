#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @authors: Patrick Kavanagh
#
"""

JDarts tools for running the JWST Detector1Pipeline. 

This modules provides a wrapper around the JWST `Detector1Pipeline`, (Stage 1 Detector processing).
It processes raw detector-level integrations (`_uncal.fits`) from JWST
into corrected countrate (slope) images (`_rate.fits`), suitable for Stage 2 processing.

This tool simplifies batch processing and allows optional step overrides, including
multiprocessing configuration.

:Dependencies:
    - jwst (JWST calibration pipeline)
    - numpy
    - multiprocessing
    - CRDS (Calibration References Data System) environment for JWST reference files

:Environment:
    This script sets the following environment variables by default:
    - CRDS_SERVER_URL
    - CRDS_CONTEXT (must be configured depending on reference file versioning)

:Inputs:
    - A directory containing _uncal.fits files, or a list of raw detector files.

:Outputs:
    - Calibrated JWST Level 1b products in the specified or default output directory.

:History:

16 Mar 2020: created

"""
import os
import sys
import glob
import shutil
import numpy as np
import warnings
import multiprocessing
from jwst import datamodels
from jwst.pipeline import Detector1Pipeline

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

# set JWST CDRS environment variables/constants
os.environ["CRDS_CONTEXT"] = ""
os.environ["CRDS_SERVER_URL"] = 'https://jwst-crds.stsci.edu'

# set number of cores to be half of those on the machine
NCORES = multiprocessing.cpu_count() // 2


class JDartsDet1:
    """
    A wrapper class for the JWST Detector1Pipeline.

    This class simplifies the use of the JWST Detector1Pipeline by enabling batch
    processing of _uncal.fits files from JWST instruments. It allows step-specific
    configuration and multiprocessing for efficient pipeline execution

    Parameters
    ----------
    input: str or list, required
        Path to a directory containing `_uncal.fits` files, or a list of input files.

    output_dir: str, optional
        Directory to save output files. Defaults to location of input files.

    Methods
    -------
    run_pipeline()
        runs the JWST Detector1Pipeline on all input files. 

    Example
    -------
    >>> my_det1 = Detector1Pipeline('my_input_dir', output_dir='output_dir')
    >>> my_det1.run_pipeline()

    :Notes:

    :Todo:

    """
    def __init__(self, input, output_dir=None):

        # check, then set input to attributes
        self.input_files = self._check_input(input)

        # If no output directory is provided, use the directory of the first input file
        if output_dir is None:
            self.output_dir = os.path.split(input[0])[0]
        else:
            self.output_dir = output_dir
            if not os.path.isdir(self.output_dir):
                os.mkdir(self.output_dir)

    def _check_input(self, input):
        """
        Check if the input is a directory or a list of files. If it's a directory,
        check if it contains any `_uncal.fits` files. If it's a list, check if it
        contains any `_uncal.fits` files. If not, raise an error.

        Parameters
        ----------
        input : str or list
            Path to directory or list of `_uncal.fits` files.

        Returns
        -------
        list
            List of input files to be processed.

        Raises
        ------
        OSError
            If the input directory does not exist.
        TypeError
            If the input is not a string or a list.

        """
        if isinstance(input, str):
            if os.path.exists(input):
                input = glob.glob(os.path.join(input, '*_uncal.fits'))
            else:
                raise OSError("The supplied directory '{}' does not exist".format(input))
        elif isinstance(input, list):
            pass
        else:
            raise TypeError("Invalid input format, must be either a string pointing to file or list of files")

        return input

    def create_step_overrides(self):
        """
        Defines overrides for individual step-specific arguments of the Detector1Pipeline.

        Returns
        -------
        step_dict
            Dictionary of step overides for the Detector1Pipeline.

        """
        step_dict = {}
        step_dict['ipc'] = {'skip': 'True'}
        step_dict['jump'] = {'maximum_cores': NCORES}
        step_dict['ramp_fit'] = {'maximum_cores': NCORES}

        return step_dict

    def run_pipeline(self, step_dict=None):
        """        
        Runs the JWST Detector1Pipeline on the input files.

        Parameters
        ----------
        step_dict : dict, optional
            A dictionary of step overrides (see create_step_overrides). If None,
            default overrides will be used.

        Returns
        -------
        None

        Outputs
        -------
        Saves calibrated files (e.g., _rate.fits, _rateints.fits) in the specified output directory.
        """
        if step_dict is None:
            step_dict = self.create_step_overrides()

        # process the files to level2B
        for f in self.input_files:
            with datamodels.RampModel(f) as dm:
                _ = Detector1Pipeline.call(dm, save_results=True, output_use_model=True, output_dir=self.output_dir,
                                           steps=step_dict)

        return None
