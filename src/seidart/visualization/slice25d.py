#!/usr/bin/env python3

"""
Create a gif for a given plane of a 2.5D model
"""

import numpy as np 
from seidart.visualization.imgen import *
from seidart.routines.definitions import *
import argparse


# !!!!! This could use some better organization of the definition

# !!!!! This needs to be moved over to FDTDImage in the imgen.py file

# ==================== Create the object and assign inputs ====================
def slice(prjfile, channel, indslice, num_steps, plane, alpha, delay):
    """
    
    """
    # We don't need the material values
    domain, material, seismic, electromag = loadproject(
        prjfile,
        Domain(),
        Material(),
        Model(),
        Model()
    )
    # First check to see if the inputs are indices or
    domain.dim = domain.dim
    # Adjust the object fields relative to the cpml
    domain.nx = int(domain.nx) + 2*int(domain.cpml)
    domain.ny = int(domain.ny) + 2*int(domain.cpml)
    domain.nz = int(domain.nz) + 2*int(domain.cpml)

    # Get the list of files to load
    all_files = glob.glob(channel + '*.dat')
    all_files.sort()

    m = len(all_files)

    if channel == 'Ex':
        NX = domain.nz
        NY = domain.ny
        NZ = domain.nx-1
    elif channel == 'Ey':
        NX = domain.nz
        NY = domain.ny-1
        NZ = domain.nx
    elif channel == 'Ez':
        NX = domain.nz-1
        NY = domain.ny
        NZ = domain.nx
    else:
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
            mag = FDTDImage(prjfile, fn)
            mag.getprjvals()
            npdat = read_dat(fn,channel,domain, single=True)
            if plane == 'xy':
                npdat = npdat[int(indslice),:,:]
            elif plane == 'xz':
                npdat = npdat[:,int(indslice),:]
            else:
                npdat = npdat[:,:,int(indslice)]
            mag.nx = npdat.shape[1]
            mag.nz = npdat.shape[0]
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
    shellcommand = 'convert -delay ' + \
        delay + ' -loop 0 magnitude.' + channel + '*.png ' + \
            channel + '.gif'
    call(shellcommand, shell = True)

    # Remove the png files 
    for filepath in glob.glob('magnitude.' + channel + '*.png'):
        os.remove(filepath)
            

# -------------------------- Command Line Arguments ---------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This program creates a csv file of time series for each 
        receiver location listed in the the specified receiver file.""" 
    )

    parser.add_argument(
        '-p', '--prjfile',
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
    prjfile = ''.join(args.prjfile)
    channel = ''.join(args.channel)
    indslice = args.index[0]
    num_steps = args.num_steps[0]
    plane = ''.join(args.plane)
    alpha = args.alpha[0]
    delay = str(args.delay[0])

    slice(prjfile, channel, indslice, num_steps, plane, alpha, delay)