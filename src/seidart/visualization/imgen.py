#!/usr/bin/env python3

"""
Functions for plotting a 2D vector/quiver plot over the model image
"""

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl

import argparse
from seidart.routines.definitions import *
import matplotlib.image as mpimg

from glob2 import glob 

# ============================ Create the objects =============================
class FDTDImage:
    """
    A class to represent an FDTD image for visualization.

    This class provides methods to build and manipulate images generated from FDTD modeling.

    :param prjfile: The full file path for the project file.
    :type prjfile: str
    :param inputfile: The name of the input file containing simulation data.
    :type inputfile: str
    :param is_complex: Indicates if the data is complex valued.
    :type is_complex: bool
    :param is_single_precision: Indicates if the data is in single precision format.
    :type is_single_precision: bool, optional
    :param plottype: The type of plot to generate ('magnitude', 'quiver', 'phase', 'energy_density').
    :type plottype: str, optional
    """
    def __init__(
            self, 
            prjfile, inputfile,
            is_complex: bool, 
            is_single_precision: bool = True,
            plottype: str = 'magnitude'
        ):
        self.prjfile = prjfile
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
        # Flags
        self.is_complex = is_complex
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
        self.domain, self.material, self.seismic, self.electromag = loadproject(
            self.prjfile,
            Domain(), 
            Material(),
            Model(),
            Model()
        )
        
        if self.is_complex:
            self.imag_part = None
        
        # Define the channel given the input file name
        self.channel = self.inputfile[0:2]
        if self.plottype == 'quiver':
            if 'E' in self.channel:
                self.xfile = 'Ex' + self.inputfile[2:]
                self.zfile = 'Ez' + self.inputfile[2:]
            else:
                self.xfile = 'Vx' + self.inputfile[2:]
                self.zfile = 'Vz' + self.inputfile[2:]
            
        self.plotfile = self.plottype + '.' + self.inputfile[2:-3] + '.png'
    
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
        # The EM code uses a slightly different mesh
        if self.channel == 'Ex':
            self.nx = self.domain.nx - 1
            self.nz = self.domain.nz
        if self.channel == 'Ez':
            self.nz = self.domain.nz - 1
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
            self.srcx = float(self.seismic.x)/self.dx + cpml + 1
            self.srcz = float(self.seismic.z)/self.dz + cpml + 1
        
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

        u = read_dat(
            self.xfile, 
            self.channel[0] + 'x', 
            self.domain, 
            is_complex = False, 
            single = self.is_single_precision
        )
        v = read_dat(
            self.zfile, 
            self.channel[0] + 'z', 
            self.domain, 
            is_complex = False, 
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
        dat = read_dat(
            self.inputfile, 
            self.channel, 
            self.domain, 
            is_complex = self.is_complex, 
            single = self.is_single_precision
        )
        if self.plottype == 'energy_density':
            # Convert the complex electric field to the 
            dat = dat.conj() * dat 
            dat = dat.real
        elif self.plottype == 'phase':
            dat = np.angle(dat)
        else:
            dat = dat.real
        
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
        if self.prjfile:
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
        pass 
    
