#!/usr/bin/env python3

"""
Functions for plotting a 2D vector/quiver plot over the model image
"""

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl

import argparse
import gzip
from pathlib import Path
from seidart.routines.definitions import *
from seidart.routines.classes import Domain, Material, Model
import matplotlib.image as mpimg

from glob2 import glob 



def is_em25_block_ref(filename: str) -> bool:
    return isinstance(filename, str) and '|EM25|' in filename


def make_em25_block_ref(filename: str, channel: str, step: int) -> str:
    return f'{filename}|EM25|{channel}|{int(step)}'


def parse_em25_block_ref(filename: str):
    block_file, _, channel, step = filename.rsplit('|', 3)
    return block_file, channel, int(step)


_EM25_BLOCK_CACHE = {}


def _open_em25_header(filename: str):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rb')
    return open(filename, 'rb')


def _read_em25_block_metadata(filename: str):
    with _open_em25_header(filename) as f:
        magic = f.read(16).decode('ascii').strip()
        if magic != 'SEIDART_EM25B1':
            raise ValueError(f'{filename} is not a SeidarT EM25 block file')
        version = int(np.frombuffer(f.read(4), dtype='<i4')[0])
        nx, ny, nz = np.frombuffer(f.read(12), dtype='<i4')
        first_step, step_count, component_count = np.frombuffer(f.read(12), dtype='<i4')
    return {
        'version': version,
        'nx': int(nx),
        'ny': int(ny),
        'nz': int(nz),
        'first_step': int(first_step),
        'step_count': int(step_count),
        'component_count': int(component_count),
    }


def em25_block_frames(channel: str, path: str = '.', pattern: str = 'EM25.*.blk*'):
    frames = []
    for block_file in sorted(glob(f'{path}/{pattern}')):
        metadata = _read_em25_block_metadata(block_file)
        first_step = metadata['first_step']
        step_count = metadata['step_count']
        for step in range(first_step, first_step + step_count):
            frames.append(make_em25_block_ref(block_file, channel, step))
    return frames


def select_frame_interval(files, frame_interval: int):
    frame_interval = int(frame_interval)
    if frame_interval < 1:
        raise ValueError('frame_interval must be greater than or equal to 1')
    return list(files)[::frame_interval]


def fdtd_output_frames(channel: str, frame_interval: int = 1, path: str = '.'):
    # Keep legacy behavior first. Only use EM25 blocks when the per-step files
    # are absent, then apply the same interval to logical timesteps.
    files = glob(f'{path}/{channel}*.dat')
    files.sort()
    if not files and channel in ['Ex', 'Ey', 'Ez']:
        files = em25_block_frames(channel, path=path)
    return select_frame_interval(files, frame_interval)


def read_fdtd_image(filename: str, channel: str, domain, single: bool = True) -> np.ndarray:
    if is_em25_block_ref(filename):
        block_file, block_channel, step = parse_em25_block_ref(filename)
        cache_key = (block_file, block_channel)
        if cache_key not in _EM25_BLOCK_CACHE:
            _EM25_BLOCK_CACHE[cache_key] = read_em25_block(block_file, channels=[block_channel])
        block = _EM25_BLOCK_CACHE[cache_key]
        time_index = step - block['metadata']['first_step']
        # EM25 blocks are stored as (time, nx, ny, nz). The legacy visualization
        # path expects (nz, ny, nx), matching read_dat for 2.5D fields.
        return np.transpose(block['fields'][block_channel][time_index], (2, 1, 0))
    return read_dat(filename, channel, domain, single=single)


def frame_label(filename: str) -> str:
    if is_em25_block_ref(filename):
        block_file, channel, step = parse_em25_block_ref(filename)
        return f'{channel}.{step:06d}.{Path(block_file).stem}'
    return filename[:-3]

# ============================ Create the objects =============================
class FDTDImage:
    """
    A class to represent an FDTD image for visualization.

    This class provides methods to build and manipulate images generated from FDTD modeling.

    :param project_file: The full file path for the project file.
    :type project_file: str
    :param inputfile: The name of the input file containing simulation data.
    :type inputfile: str
    :param is_single_precision: Indicates if the data is in single precision format.
    :type is_single_precision: bool, optional
    :param plottype: The type of plot to generate ('magnitude', 'quiver', 'phase', 'energy_density').
    :type plottype: str, optional
    """
    def __init__(
            self, 
            project_file, inputfile,
            is_single_precision: bool = True,
            plottype: str = 'magnitude',
            numerical_method: str = 'fdtd'
        ):
        self.project_file = project_file
        self.x = None
        self.z = None
        self.inputfile = inputfile
        self.srcx = None
        self.srcz = None
        self.channel = None
        # Plot values
        self.extent = None
        self.background = None
        self.ax = None 
        self.fig = None
        self.dx = None 
        self.dz = None
        self.nx = None 
        self.nz = None
        self.numerical_method = numerical_method
        # Flags
        self.is_single_precision = is_single_precision
        
        # Check the inputs that need to be checked
        if plottype not in [
            'magnitude', 'quiver', 'phase', 'energy_density'
        ]:
            raise ValueError(
                "plottype must be one of \
                    'magnitude', 'quiver', 'phase', 'energy_density'"
            )
        
        self.plottype = plottype
        self.build()
        self.getprjvals()
         
    def build(self):
        """
        Builds the domain, material, seismic, and electromagnetic models from the project file.
        """
        if self.numerical_method == 'dg':
            model_seismic = Biot()
        else:
            model_seismic = Model()
        
        self.domain, self.material, self.seismic, self.electromag = loadproject(
            self.project_file,
            Domain(), 
            Material(),
            model_seismic,
            Model()
        )
        
        if self.numerical_method == 'dg':
            self.domain.cpml = 0 
        
        # Define the channel given the input file name
        if is_em25_block_ref(self.inputfile):
            block_file, self.channel, step = parse_em25_block_ref(self.inputfile)
        else:
            block_file = None
            step = None
            self.channel = self.inputfile[0:2]
        if self.plottype == 'quiver':
            if 'E' in self.channel:
                if block_file is None:
                    self.xfile = 'Ex' + self.inputfile[2:]
                    self.zfile = 'Ez' + self.inputfile[2:]
                else:
                    self.xfile = make_em25_block_ref(block_file, 'Ex', step)
                    self.zfile = make_em25_block_ref(block_file, 'Ez', step)
            else:
                self.xfile = 'Vx' + self.inputfile[2:]
                self.zfile = 'Vz' + self.inputfile[2:]
            
        self.plotfile = self.plottype + '.' + frame_label(self.inputfile) + '.png'
    
    def getprjvals(self):
        """
        Retrieves project values and initializes domain parameters.
        """
        # Let's initiate the domain
        self.domain.cpml = int(self.domain.cpml)
        self.domain.nx = int(self.domain.nx) + 2*self.domain.cpml
        self.domain.nz = int(self.domain.nz) + 2*self.domain.cpml
        self.dx = float(self.domain.dx)
        self.dz = float(self.domain.dz)
        # # The EM code uses a slightly different mesh
        # if self.channel == 'Ex':
        #     self.nx = self.domain.nx
        #     self.nz = self.domain.nz
        # elif self.channel == 'Ez':
        #     self.nz = self.domain.nz - 1
        #     self.nx = self.domain.nx
        # else:
        self.nz = self.domain.nz 
        self.nx = self.domain.nx
        
        self.extent = (
            -self.domain.cpml, 
            (self.nx-self.domain.cpml), 
            (self.nz-self.domain.cpml), 
            -self.domain.cpml
        )

        if self.channel == 'Ex' or self.channel == 'Ez':
            self.srcx = float(self.electromag.x)/self.dx + self.domain.cpml + 1
            self.srcz = float(self.electromag.z)/self.dz + self.domain.cpml + 1
        else:    
            self.srcx = float(self.seismic.x)/self.dx + self.domain.cpml + 1
            self.srcz = float(self.seismic.z)/self.dz + self.domain.cpml + 1
        
        # Define tick locations for plotting
        self.xticklocs = np.array(
            [
                0, 
                (self.nx-2*self.domain.cpml)/4, 
                int((self.nx-2*self.domain.cpml)/2), 
                3*(self.nx-2*self.domain.cpml)/4, 
                self.nx - 2*self.domain.cpml
            ]
        )
        self.yticklocs = np.array(
            [
                0, 
                (self.nz-2*self.domain.cpml)/4, 
                int((self.nz-2*self.domain.cpml)/2), 
                3*(self.nz-2*self.domain.cpml)/4, 
                self.nz - 2*self.domain.cpml
            ]
        )
        # Load the model image and assign variables
        self.background = mpimg.imread(self.domain.imfile)
        
    # -------------------------------------------------------------------------
    def quiverplot(self, papercolumnwidth: float = 7.2) -> None:
        """
        Generates a quiver plot for seismic outputs.

        :param papercolumnwidth: The width of the paper column for the plot.
        :type papercolumnwidth: float, optional
        """
        # buildmesh
        # Get axes values
        x = (np.arange(-self.domain.cpml, self.nx-self.domain.cpml))
        z = (np.arange(-self.domain.cpml, self.nz-self.domain.cpml))
        # Append the cpml values 
        x, z = np.meshgrid(x,z)

        u = read_fdtd_image(
            self.xfile,
            self.channel[0] + 'x',
            self.domain,
            single = self.is_single_precision
        )
        v = read_fdtd_image(
            self.zfile,
            self.channel[0] + 'z',
            self.domain,
            single = self.is_single_precision
        )        
        # Set the figure size to be for a full two column width
        
        # Create the figure and axes objects
        self.fig = plt.figure(
            figsize = [
                papercolumnwidth, 
                self.nz*papercolumnwidth/self.nx
            ]
        )
        self.ax = plt.gca()
        
        # add the model 
        self.ax.imshow(
            self.background,
            alpha = 0.7, 
            extent=[
                0, self.nx-2*self.domain.cpml, 
                self.nz-2*self.domain.cpml, 0
            ]
        )
        
        # add quiver
        q = self.ax.quiver(
            x, z, u, v, 
            headwidth = 0.5, 
            headlength = 1,
            headaxislength = 1, 
            scale = 16, 
            minlength = 0.1
        )
    
    # -------------------------------------------------------------------------
    def magnitudeplot(
            self, alpha: float = 0.3, papercolumnwidth: float = 7.2
        ) -> None:
        """
        Plots the real data as a magnitude plot.

        :param alpha: Transparency of the model plotted in the background.
        :type alpha: float, optional
        :param papercolumnwidth: The width of the paper column for the plot.
        :type papercolumnwidth: float, optional
        """
        dat = read_fdtd_image(
            self.inputfile,
            self.channel,
            self.domain,
            single = self.is_single_precision
        )
        
        self.fig = plt.figure(
            figsize = [
                papercolumnwidth, 
                self.nz*papercolumnwidth/self.nx
            ]
        )
        self.ax = plt.gca()
        self.ax.imshow(
            dat, 
            cmap = 'seismic',
            extent = self.extent,
            norm = mpl.colors.CenteredNorm()
        )
        self.ax.imshow(
            self.background,
            alpha = alpha,
            extent = [
                0, self.nx-2*self.domain.cpml, 
                self.nz-2*self.domain.cpml, 0
            ]
        )
        
    # -------------------------------------------------------------------------
    def sliceplot(
            self, 
            dat: np.ndarray, 
            axscale: tuple, 
            sliceaxes: str, 
            alpha: float = 0.3, 
            papercolumnwidth: float = 7.2
        ) -> None:
        """
        Generates a slice plot from the given data.

        :param dat: The data array to plot.
        :type dat: np.ndarray
        :param axscale: The scale for the axes.
        :type axscale: tuple
        :param sliceaxes: The axes to slice ('xy', 'xz', 'yz').
        :type sliceaxes: str
        :param alpha: Transparency of the model plotted in the background.
        :type alpha: float, optional
        :param papercolumnwidth: The width of the paper column for the plot.
        :type papercolumnwidth: float, optional
        """
        # sliceaxes is a string that specifies 'xy', 'xz', or 'yz'
        # All of the slices need to be scaled equally and normalized
        # axscale = (nz, ny, nx)/max(nz, ny, nx)
        if sliceaxes == 'xy': # dat dims will be y then x
            figuredims = [
                axscale[1]*papercolumnwidth, axscale[2]*papercolumnwidth
            ]
            exaggeration = axscale[2]/axscale[1]
        elif sliceaxes == 'xz': # Dims will be z then x
            figuredims = [
                axscale[0]*papercolumnwidth, axscale[2]*papercolumnwidth
            ]
            exaggeration = axscale[2]/axscale[0]
        else:
            figuredims = [
                axscale[0]*papercolumnwidth, axscale[1]*papercolumnwidth
            ]
            exaggeration = axscale[1]/axscale[0]
        
        self.fig = plt.figure(
            figsize = figuredims
        )
        self.ax = plt.gca()
        self.ax.imshow(
            dat,
            cmap = 'seismic', 
            extent = self.extent,
            norm = mpl.colors.CenteredNorm()
        )
        self.ax.set_aspect(exaggeration)
        self.ax.imshow(
            self.background,
            alpha = alpha,
            extent = [
                0, self.nx-2*self.domain.cpml, 
                self.nz-2*self.domain.cpml, 0
            ]
        )
    
    # -------------------------------------------------------------------------
    def imvector(self) -> None:
        """
        Generates a vector image from the project values and saves it.
        """
        self.getprjvals()
        self.quiverplot()
        self.addlabels()
        plt.savefig(self.plotfile)
        plt.close()
    
    # This is setup for seismic
    def vectoranim(self, files: list = glob('Vx*.dat')) -> None:
        """
        Creates an animation of vector images from a list of files.

        :param files: A list of file names to create the animation from.
        :type files: list, optional
        """
        # Project file needs to already be assigned in the object
        if self.project_file:
            files.sort()

            print('Creating PNG snapshots')
            # We want to do the first file; initial condition
            n=num_steps
            for fn in files:
                if n == num_steps:
                    self.inputfile = fn
                    self.imvector()
                    n = 1
                else:
                    n = n + 1
            # Create the gif using Imagemagick
            print('Creating the GIF')
            call('convert -delay 20 -loop 0 vector.*.png vector.gif', shell = True)
        else:
            print('No project file has been assigned')
    
    # -------------------------------------------------------------------------
    
    # -------------------------------------------------------------------------
    
    # -------------------------------------------------------------------------
    def addlabels(self) -> None:
        """
        Adds labels to the plot and updates the axis and figure objects.
        """
        # Set axes labels
        xticklabels = (self.xticklocs)*self.dx
        yticklabels = (self.yticklocs)*self.dz
        self.ax.set_xlabel(r'X (m)')
        self.ax.xaxis.tick_top()
        self.ax.xaxis.set_label_position('top')
        self.ax.set_xticks(self.xticklocs)
        self.ax.set_xticklabels(xticklabels)
        self.ax.set_ylabel(r'Z (m)')
        self.ax.set_yticks(self.yticklocs)
        self.ax.set_yticklabels(yticklabels)
        
    def addsource(self) -> None:
        """
        Adds the source location to the plot.
        """
        # Source location
        self.ax.scatter(
            self.srcx, 
            self.srcz,
            marker = '*', 
            s = 30, 
            linewidths = 1,
            edgecolor = (0.2, 0.2, 0.2, 1 ) 
        )

    def addrcx(self):
        """
        Placeholder for adding receiver markers to the current plot.

        :return: None
        :rtype: None
        """
        pass 
    
