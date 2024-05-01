#!/usr/bin/env python3

# From the set of image outputs, we can build a gif. The images will be in csv
# or fortran unformatted binary


import numpy as np
from glob2 import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.io import FortranFile
from seidart.visualization.imgen import *
import os
import subprocess

# =============================================================================
def build_animation(
        prjfile, 
        channel, 
        delay, 
        numsteps, 
        alpha, 
        is_complex: bool = False,
        is_single_precision: bool = True,
        plottype: str = 'magnitude'
    ):
    '''
    This function builds a gif from the set of images output of the FDTD modeling. 
    The images can be in csv or unformatted Fortran binary, however, the program 
    runs faster to use the latter. 
    
    :param prjfile: The full file path for the project file.
    :type prjfile: str
    :param channel: The channel to be used. Available channels are Ex, Ez, Vx, and Vz 
                    for the electric field and seismic velocities, respectively.
    :type channel: str
    :param delay: The amount of delay between two frames in the GIF animation.
    :type delay: int
    :param numsteps: The time step interval between the images to be used. 
                     For example, n=20 means that every 20th image will be used.
    :type numsteps: int
    :param alpha: Change the transparency of the model plotted in the background. 
                  0 is fully transparent, 1 is fully opaque.
    :type alpha: float
    :param is_complex: Flag indicating whether the data will be complex valued. 
                       If not flagged but the data is complex, only the real part is used.
    :type is_complex: bool, optional
    :param is_single_precision: Flag indicating whether the data is in single precision.
    :type is_single_precision: bool, optional
    :param plottype: The type of plot to generate. Valid inputs are 'magnitude', 
                     'phase' (complex only), and 'energy_density' (complex only) for EM; 
                     'magnitude' and 'quiver' for seismic.
    :type plottype: str, optional
    :return: None

        EM
        --    
            'magnitude', 
            'phase' (complex only), 
            'energy_density (complex only)'
        Seismic
        -------
            'magnitude',
            'quiver'
    
    Returns
    -------
    None
    '''
    # Check if the .dat files are still around
    files = glob(channel + '*.dat')
    files.sort()

    # We'll start counting with the first frame
    n=numsteps

    for fn in files:
        if n == numsteps:
            mag = FDTDImage(
                prjfile, fn, 
                is_complex = is_complex, 
                is_single_precision = is_single_precision,
                plottype = plottype
            )
            mag.magnitudeplot(alpha = alpha)
            mag.addlabels()
            mag.plotfile = 'magnitude.' + fn[:-3] + '.png'
            plt.savefig(mag.plotfile)
            plt.close()
            n = 1
        else:
            n = n + 1


    print('Creating the GIF')
    # Use imagemagick via shell command to create the gif
    shellcommand = 'convert -delay ' + \
        str(delay) + ' -loop 0 magnitude.' + channel + '*.png ' + \
            channel + '.gif'
    subprocess.call(shellcommand, shell = True)

    # Remove the png files 
    for filepath in glob('magnitude.' + channel + '*.png'):
        os.remove(filepath)
        

# -------------------------- Command Line Arguments ---------------------------
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="""This program builds a gif from
        the set of image output of the FDTD modeling. The images can be in csv or
        unformatted Fortran binary, however, the program runs faster to use the
        latter. """ )

    parser.add_argument(
        '-p','--prjfile', 
        nargs=1, type=str, required = True,
        help='the full file path for the project file'
    )


    parser.add_argument(
        '-n', '--num_steps', 
        nargs = 1, type = int, required = True, 
        help = """The time step interval between the images that
        are going to be used. Every time step is written to file which means that
        we can take any equally spaced images to create the gif with an
        appropriate resolution, time to compute, and file size. For example,
        n=20 means that every 20 images will be used thus significantly reducing
        how long it takes to compile the gif."""
    )

    parser.add_argument(
        '-c', '--channel', 
        nargs = 1, type = str, required = True,
        help = """Specify whether a particular channel is going to be used. The
        available channels are Ex, Ez, Vx, and Vz for the electric field and
        seismic velocities, respectively."""
    )

    parser.add_argument(
        '-d', '--delay', 
        nargs = 1, type = int, required = False, default = [1], 
        help = """The amount of delay between two frames"""
    )

    parser.add_argument(
        '-a', '--alpha',
        nargs = 1, type = float, required = False, default = [0.3],
        help = """(OPTIONAL FLOAT [0,1]) Change the transparency of the model 
        plotted in the background; default = 0.3. Zeros is full transparency, 1 is 
        CIA transparency."""
    )
    
    parser.add_argument(
        '-z', '--is_complex', action = 'store_true', required = False,
        help = """Flag whether the data will be complex valued. If the data is
        not flagged but is complex, only the real data will be returned. """
    )
    
    parser.add_argument(
        '-D', '--is_double_precision', action = 'store_false', required = False,
        help = """Flag whether the data is in double precision."""
    )
    
    parser.add_argument(
        '-P', '--plot_type', required = False, default = ['magnitude'], 
        nargs = 1, type = str, 
        help = """
        Plot types are different for EM and seismic. Valid inputs are: 
        
        EM
        --    
            'magnitude', 
            'phase' (complex only), 
            'energy_density (complex only)'
        Seismic
        -------
            'magnitude',
            'quiver'
        """
    )

    # Get the arguments
    args = parser.parse_args()
    prjfile = ''.join(args.prjfile)
    channel = ''.join(args.channel)
    delay = str(args.delay[0])
    num_steps = args.num_steps[0]
    alpha = min([1, args.alpha[0]] )
    is_complex = args.is_complex
    is_single_precision = args.is_single_precision
    plot_type = args.plot_type[0]
    
    build_animation(
        prjfile, 
        channel, 
        delay, 
        numsteps, 
        alpha, 
        is_complex,
        is_single_precision,
        plot_type
    )