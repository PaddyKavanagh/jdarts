#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @authors: Patrick Kavanagh
#
"""

JDarts tools for Detector1Pipeline


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

# settings
warnings.filterwarnings("ignore")

# set environment variables/constants
os.environ["CRDS_CONTEXT"] = ""
os.environ["CRDS_SERVER_URL"] = 'https://jwst-crds.stsci.edu'

# set number of cores to be half of those on the machine
NCORES = multiprocessing.cpu_count() // 2


class JDartsDet1:
    """
    Wrapper for Detector1Pipeline

    :Parameters:

    input: str or list, required
        The directory containing input files or list of input files

    output_dir: str, optional
        Location to save the output. Defaults to location of input files

    :Methods:

    run_pipeline:
        runs the pipeline

    :Examples:

    my_det1 = Detector1Pipeline('my_input_dir', output_dir='output_dir')
    my_det1.run_pipeline()

    :Notes:

    :Todo:

    """
    def __init__(self, input, output_dir=None):

        # check, then set input to attributes
        self.input_files = self._check_input(input)

        # set the output directories
        if output_dir is None:
            self.output_dir = os.path.split(input[0])[0]
        else:
            self.output_dir = output_dir
            if not os.path.isdir(self.output_dir):
                os.mkdir(self.output_dir)

    def _check_input(self, input):
        """
        Check the input file is valid
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
        Include any step-specific arguments here
        """
        step_dict = {}
        step_dict['ipc'] = {'skip': 'True'}
        step_dict['jump'] = {'maximum_cores': NCORES}
        step_dict['ramp_fit'] = {'maximum_cores': NCORES}

        return step_dict

    def run_pipeline(self, step_dict=None):
        """
        Process each file
        """
        if step_dict is None:
            step_dict = self.create_step_overrides()

        # process the files to level2B
        for f in self.input_files:
            with datamodels.RampModel(f) as dm:
                _ = Detector1Pipeline.call(dm, save_results=True, output_use_model=True, output_dir=self.output_dir,
                                           steps=step_dict)

        return None
