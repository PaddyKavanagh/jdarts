.. jdarts documentation master file, created by
   sphinx-quickstart on Mon Mar 24 18:17:49 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*******************************
JDARTS Documentation
*******************************

Introduction
============

The *jdarts* (JWST Data Analysis and Reduction Tools) package is a Python library 
for the data reduction and analysis of data from the James Webb Space Telescope (JWST).

You can install *jdarts* with `pip <http://www.pip-installer.org/en/latest/>`_::

The package is organized into several sub-packages, each with its own set of modules containing functions.
The sub-packages are:

   - **jdarts.pipeline**: Contains JWST pipeline data reduction scripts.
   - **jdarts.image**: Contains functions for image processing and analysis.
   - **jdarts.spectro**: Contains functions for spectroscopic processing and analysis.
   - **jdarts.utils**: Contains shared utility functions.

.. toctree::
   :maxdepth: 2
   :caption: Package Structure:

   jdarts/pipeline/index
   jdarts/image/index
   jdarts/spectro/index
   jdarts/utils/index
   
----

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
