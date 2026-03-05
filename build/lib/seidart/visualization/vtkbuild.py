#!/usr/bin/env python3

#

import numpy as np
import glob
import argparse
from scipy.io import FortranFile
from seidart.routines.definitions import *
from pyevtk.hl import imageToVTK



# ------------------------------ Run the program ------------------------------
def build():
    domain, material, seismic, electromag = loadproject(
        project_file,
        Domain(),
        Material(),
        Model(),
        Model()
    )

    # Define some plotting inputs
    domain.cpml = int(domain.cpml)
    domain.nx = int(domain.nx) + 2*domain.cpml
    domain.ny = int(domain.ny) + 2*domain.cpml
    domain.nz = int(domain.nz) + 2*domain.cpml
    ncells = domain.nx*domain.ny*domain.nz

    X = np.linspace(1, domain.nx, num = domain.nx)*domain.dx
    Y = np.linspace(1, domain.nx, num = domain.ny)*domain.dy
    Z = np.linspace(domain.nz, 1, num = domain.nz)*domain.dz
    x = np.zeros([domain.nx, domain.ny, domain.nz])
    y = np.zeros([domain.nx, domain.ny, domain.nz])
    z = np.zeros([domain.nx, domain.ny, domain.nz])

    for i in range(0, domain.nx):
        for j in range(0, domain.ny):
            for k in range(0, domain.nz):
                x[i,j,k] = X[i]
                y[i,j,k] = Y[j]
                z[i,j,k] = Z[k]

    # Add the source location to plot
    if channel == 'Ex' or channel == 'Ey' or channel == 'Ez':
        electromag.x = float(electromag.x)
        electromag.y = float(electromag.y)
        electromag.z = float(electromag.z)
        ex = electromag.x/domain.dx + domain.cpml+1
        ey = electromag.y/domain.dy + domain.cpml+1
        ez = electromag.z/domain.dz + domain.cpml+1
        source_location = np.array([ex, ey, ez])
        dt = float(electromag.dt)
    else:
        seismic.x = float(seismic.x)
        seismic.y = float(seismic.y)
        seismic.z = float(seismic.z)
        sx = seismic.x/domain.dx + domain.cpml+1
        sy = seismic.y/domain.dy + domain.cpml+1
        sz = seismic.z/domain.dz + domain.cpml+1
        source_location = np.array([sx, sy, sz])
        dt = float(seismic.dt)

    # Proceed accordingly to the channel flag

    # Check if the .dat files are still around
    files = glob.glob(channel + '*.dat')
    ind = 0
    files.sort()

    # We'll start counting with the first frame
    n=num_steps
    for fn in files:
        if n == num_steps:
            dat = read_dat(fn,channel,domain, single=True)
            
            # Zero out any values below our given threshold
            duration = dt*ind

            # Reset the countern = 1
            ind = ind + 1
            n = 1

            vtkfilename = "./image" + channel + "." + str(ind)
            imageToVTK(vtkfilename, cellData = {"Displacement" : dat} )

        else:
            ind = ind + 1
            n = n + 1


def main(
        project_file: str,
        channel: str,
        numsteps: int,
        
    ):
    """
    This main function ...
    
    :param project_file: Path to the project file.
    :type project_file: str
    :param channel: Data channel to be analyzed.
    :type channel: str
    :param numsteps: The time step interval between the images to be used. 
                     For example, n=20 means that every 20th image will be used.
    :type numsteps: int
    :return: None
    """
    # -------------------------- Command Line Arguments ---------------------------
    parser = argparse.ArgumentParser(description="""This program builds .VTI
        (Visualization Toolkit Image) files from the 3d array outputs of the FDTD
        modeling. These files can be displayed using Paraview.""" )

    parser.add_argument(
        '-p', '--project_file',
        nargs=1, type=str, required = True,
        help='the full file path for the project file'
    )

    parser.add_argument(
        '-c', '--channel',
        nargs = 1, type = str, required = True,
        help = """Specify whether a particular channel is going to be used. The
        available channels are Ex, Ez, Vx, and Vz for the electric field and
        seismic velocities, respectively."""
    )

    parser.add_argument(
        '-n', '--numsteps',
        nargs = 1, type = int, required = True,
        help = """The time step interval between the images that
        are going to be used. Every time step is written to file which means that
        we can take any equally spaced images to create the gif with an
        appropriate resolution, time to compute, and file size. For example,
        n=20 means that every 20 images will be used thus significantly reducing
        how long it takes to compile."""
    )

    #-- Get the arguments
    args = parser.parse_args()
    project_file = ''.join(args.project_file)
    channel = ''.join(args.channel)
    numsteps = args.num_steps[0]