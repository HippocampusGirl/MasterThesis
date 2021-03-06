# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
Base module variables
"""

__version__ = "0.0.1-dev"

__author__ = ""
__copyright__ = ""
__credits__ = [""]
__license__ = ""
__maintainer__ = "Lea Waller"
__email__ = "lea.waller@charite."
__status__ = ""
__url__ = ""
__packagename__ = "autaptomatic"
__description__ = ("")
__longdesc__ = ("")

DOWNLOAD_URL = ("")

SETUP_REQUIRES = [
    "setuptools>=18.0",
    "numpy",
]

REQUIRES = [
    "numpy",
    "scipy",
    "neo",
    "matplotlib",
    "pandas",
    "openpyxl",
    #"wxPython",
]

LINKS_REQUIRES = [
]

TESTS_REQUIRES = [
]

EXTRA_REQUIRES = {
}

# Enable a handle to install all extra dependencies at once
EXTRA_REQUIRES["all"] = [val for _, val in list(EXTRA_REQUIRES.items())]
CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
]