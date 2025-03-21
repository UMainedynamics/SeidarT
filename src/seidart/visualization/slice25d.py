#!/usr/bin/env python3

"""
Create a gif for a given plane of a 2.5D model
"""

import numpy as np 
from seidart.visualization.imgen import *
from seidart.routines.definitions import *
import argparse
from glob2 import glob
import subprocess 
import os


# !!!!! This could use some better organization of the definition

# !!!!! This needs to be moved over to FDTDImage in the imgen.py file

# ==================== Create the object and assign inputs ====================
def slicer(
            project_file, channel, indslice, num_steps, delay,
            plane = 'xz', alpha = 0.3, is_single = True
        ):
    """
    Create a gif for a given plane of a 2.5D model. The function will create a
    series of images and then use imagemagick to create a gif.

    :param project_file: str: The project file path
    :type project_file: str
    :param channel: str: The channel to query
    :type channel: str
    :param indslice: int: The index along the plane that we are slicing. *This 
        isn't corrected for the boundary layer. 
    :type indslice: int
    :param num_steps: int: The number of time steps between frames
    :type num_steps: int
    :param plane: str: Specify what plane to slice (default = 'xz'). The slice is the
        along the 3rd plane (i.e. y-ind for a 'xz' plane). Options are 'xy', 'yz', and 'xz'.
    :type plane: str
    :param alpha: float: (OPTIONAL FLOAT [0,1]) Change the transparency of the model
        plotted in the background; default = 0.3. Zeros is full transparency, 1 is
        CIA transparency
    :type alpha: float
    :param delay: int: The amount of delay between two frames
    :type delay: int
    :param is_single: bool: If the data is complex
    :type is_single: bool
    :return: None

    """
    # We don't need the material values
    domain, material, seismic, electromag = loadproject(
        project_file,
        Domain(),
        Material(),
        Model(),
        Model()
    )
    # First check to see if the inputs are indices or
    domain.dim = domain.dim
    # Adjust the object fields relative to the cpml
    domain.nx = int(domain.nx + 2*domain.cpml)
    domain.ny = int(domain.ny + 2*domain.cpml)
    domain.nz = int(domain.nz + 2*domain.cpml)

    # Get the list of files to load
    all_files = glob(channel + '*.dat')
    all_files.sort()

    m = len(all_files)

    # if channel == 'Ex':
    #     NX = domain.nz
    #     NY = domain.ny
    #     NZ = domain.nx-1
    # elif channel == 'Ey':
    #     NX = domain.nz
    #     NY = domain.ny-1
    #     NZ = domain.nx
    # elif channel == 'Ez':
    #     NX = domain.nz-1
    #     NY = domain.ny
    #     NZ = domain.nx
    # else:
    NX = domain.nz
    NY = domain.ny 
    NZ = domain.nx

    # Pre allocate depending on what dimension we are slicing
    if plane == 'xy':
        imageseries = np.zeros([NY, NZ, m])
    elif plane == 'xz':
        imageseries = np.zeros([NX, NZ, m])
    else:
        imageseries = np.zeros([NX, NY, m])


    # We'll start counting with the first frame
    n=num_steps
    axscale = np.array([domain.nz, domain.ny, domain.nx])
    axscale = axscale/axscale.max() 

    for ind, fn in enumerate(all_files, start = 0):
        if n == num_steps:
            mag = FDTDImage(project_file, fn)
            mag.getprjvals()
            npdat = read_dat(
                fn, channel, domain, single=is_single
            )
            if plane == 'xy':
                npdat = npdat[int(indslice),:,:]
            elif plane == 'yz':
                npdat = npdat[:,:,int(indslice)]
            else:
                npdat = npdat[:,int(indslice),:]
            # mag.nx = npdat.shape[1]
            # mag.nz = npdat.shape[0]
            mag.sliceplot(npdat, axscale, plane, alpha = alpha)
            mag.addlabels()
            mag.plotfile = 'magnitude.' + fn[:-3] + '.png'
            plt.savefig(mag.plotfile)
            plt.close()
            n = 1
        else:
            n = n + 1


    print('Creating the GIF')
    # Use imagemagick via shell command to create the gif
    shellcommand = 'magick -delay ' + \
        str(delay) + ' -loop 0 magnitude.' + channel + '*.png ' + \
            channel + '.gif'
    subprocess.call(shellcommand, shell = True)

    # Remove the png files 
    for filepath in glob('magnitude.' + channel + '*.png'):
        os.remove(filepath)
            

# -------------------------- Command Line Arguments ---------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This program creates a csv file of time series for each 
        receiver location listed in the the specified receiver file.""" 
    )

    parser.add_argument(
        '-p', '--project_file',
        nargs = 1, type = str, required = True,
        help = 'The project file path'
    )

    parser.add_argument(
        '-i', '--index',
        nargs=1, type=int, required = True,
        help='The index along the plane that we are slicing.'
    )

    parser.add_argument(
        '-c', '--channel',
        nargs = 1, type = str, required = True,
        help = """The channel to query. """
    )

    parser.add_argument(
        '-P', '--plane',
        nargs = 1, type = str, default = ['xz'], required = False,
        help = """Specify what plane to slice (default = 'xy'). The slice is the 
        along the 3rd plane. """
    )

    parser.add_argument(
        '-n', '--num_steps', 
        nargs = 1, type = int, required = True,
        help = """Specify the number of time steps between frames"""

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

    # Get the arguments
    args = parser.parse_args()
    project_file = ''.join(args.project_file)
    channel = ''.join(args.channel)
    indslice = args.index[0]
    num_steps = args.num_steps[0]
    plane = ''.join(args.plane)
    alpha = args.alpha[0]
    delay = str(args.delay[0])

    slice(project_file, channel, indslice, num_steps, plane, alpha, delay)