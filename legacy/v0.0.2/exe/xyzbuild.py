#!/usr/bin/env python3

"""
Get all of the xyz values for a given RGB and save them in a specified output
file. 
"""

import argparse
import numpy as np
import matplotlib.image as mpimg
from subprocess import call
import material_functions as mf
from definitions import *


# -------------------------- Command Line Arguments ---------------------------
parser = argparse.ArgumentParser(
    description="""The SeidarT software requires a
    .PNG image that is used to construct the model domain for seismic and
    electromagnetic wave propagation. Given the image file, a project file
    will be constructed which contains all the necessary parameters to be
    read in to the finite differences time domain modeling schemes."""
)

parser.add_argument(
    '-p', '--prjfile', nargs=1, type=str, required = True,
    help='the full file path for the project file', default=None
)

args = parser.parse_args()
prjfile = ''.join(args.prjfile)
# =============================================================================
domain, material, seismic, electromag = loadproject(
    prjfile,
    Domain(), 
    Material(),
    Model(),
    Model()
)


# -----------------------------------------------------------------------------
