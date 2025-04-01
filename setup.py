"""
Setup file for installing the jdarts software

@author: Patrick Kavanagh
"""

try:
    from setuptools import setup
except ImportError:
    raise ImportError("setuptools is required to install this package. You can install it with 'pip install setuptools'")

__version__ = '0.1.0'

setup(  
    name="jdarts",
    version=__version__,
    description="jdarts (JWST data analysis and reduction tools). A package for the analysis and reduction of JWST data.",
    packages=['jdarts'],
    author="Patrick Kavanagh",
    author_email = "patrick.kavanagh@mu.ie",
    license = "See licence file",
    url="https://github.com/PaddyKavanagh/jdarts",
    platforms = ["Linux", "Mac OS X"],
    python_requires='>=3.5',
    install_requires=['numpy', 
                      'astropy', 
                      'scipy', 
                      'matplotlib', 
                      'pytest', 
                      'photutils', 
                      'scikit-image', 
                      'jwst',  
                      'jwst_reffiles', 
                      'jwedb', 
                      'pandas', 
                      'requests', 
                      'tabulate', 
                      'astroquery>=0.4.2'],
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research/DataAnalysis",
        "Topic :: Scientific/Engineering"
    ]
)