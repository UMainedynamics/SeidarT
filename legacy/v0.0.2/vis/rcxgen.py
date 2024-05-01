#!/usr/bin/env python3

import numpy as np 
from definitions import *
import argparse 


parser = argparse.ArgumentParser(
    description=
    """
    Receiver locations can be specified from a PNG image using a unique RGB 
    value. The image must be the same dimensions as the model without absorbing
    boundaries. 
    """
)

parser.add_argument(
    '-p', '--prjfile', nargs=1, type=str, required = True,
    help='the full file path for the project file', default=None
)

parser.add_argument(
    '-i', '--image', nargs = 1, type = str, required = False,
    help = """Specify whether to run the seismic (s), or electromagnetic (e), 
    or none (default = n)""",
    default = 'n'
)

parser.add_argument(
    '-c', '--color',
    nargs = 3, type = int, required = False, default = [1],
    help = """The integer values [0,255] for the RGB color associated with the.
    the receiver locations."""
)

# Get the arguments
args = parser.parse_args()
prjfile = ''.join(args.prjfile)
imfilename = ''.join(args.image)
rgb = args.color 

domain, material, seismic, electromag = loadproject(
    prjfile,
    Domain(), 
    Material(),
    Model(),
    Model()
)
rcxgen(rgb, domain, material)