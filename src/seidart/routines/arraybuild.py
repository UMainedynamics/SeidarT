import numpy as np
import pandas as pd
from glob2 import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from typing import Tuple
from seidart.routines.definitions import *

# =============================================================================
class Array:
    def __init__(
        self, 
        channel: str,
        prjfile: str, 
        receiver_file: str,
        receiver_indices: bool = False, 
        single_precision: bool = True,
        is_complex: bool = False
    ):
        """
        Initializes the Array object with project settings and receiver details.
        
        :param channel: The data channel to be analyzed.
        :type channel: str
        :param prjfile: Path to the project file.
        :type prjfile: str
        :param receiver_file: Path to the receiver locations file.
        :type receiver_file: str
        :param receiver_indices: Flag if receiver file contains indices (True) 
                                 or Cartesian coordinates (False).
        :type receiver_indices: bool
        :param single_precision: Use single precision for numerical data 
            if True.
        :type single_precision: bool
        :param is_complex: Treat the data as complex if True.
        :type is_complex: bool
        """
        self.prjfile = prjfile 
        self.channel = channel
        self.receiver_file = receiver_file
        self.receiver_indices = receiver_indices
        self.single_precision = single_precision
        self.is_complex = is_complex
        self.stream = None
        self.gain = None
        self.exaggeration = 0.5
        self.build()
    
    # -------------------------------------------------------------------------
    def build(self):
        """
        Constructs domain and models based on project file and sets up source 
        and receiver configurations.
        """
        self.domain, self.material, self.seismic, self.electromag = loadproject(
            self.prjfile,
            Domain(), 
            Material(),
            Model(),
            Model()
        )
        if self.channel in ['Vx','Vy','Vz']:
            self.is_seismic = True
            self.dt = self.seismic.dt 
            self.gain = self.seismic.time_steps
        else:
            self.is_seismic = False
            self.dt = self.electromag.dt
            self.gain = self.electromag.time_steps
        
         
        # __, self.receivers_xyz = self.loadxyz()
        self.loadxyz()
        if self.channel == 'Vx' or self.channel == 'Vy' or self.channel == 'Vz':
            if self.domain.dim == '2.0':
                self.source = np.array(
                    [self.seismic.x, self.seismic.z]
                )
            else:
                self.source = np.array(
                    [
                        self.seismic.x, 
                        self.seismic.y, 
                        self.seismic.z
                    ]
                )
        else:
            if self.domain.dim == '2.0':
                self.source = np.array(
                    [self.electromag.x, self.electromag.z]
                )
            else:
                self.source = np.array(
                    [
                        self.electromag.x, 
                        self.electromag.y, 
                        self.electromag.z
                    ]
                )
        
        # Load the time series for all receivers
        self.getrcx()
    
    # -------------------------------------------------------------------------
    def loadxyz(self):
        """
        Loads and sorts receiver locations from a file and adjusts them 
        according to the domain and CPML layer.  If the source_file flag is 
        True, the source file is loaded instead of the receiver file and saved
        to the object under source_xyz.

        :param source_file: Flag to indicate if the file is a source file.
        :type source_file: bool 

        """
        xyz = pd.read_csv(self.receiver_file)
        
        # We need to make sure the recievers are ordered correctly and the 
        # absorbing boundary is corrected for
        # First check to see if the inputs are indices or
        cpml = int(self.domain.cpml)
        # Adjust the object fields relative to the cpml. The y-direction will be
        # adjusted when we are 2.5/3D modeling
        self.domain.nx = self.domain.nx + 2*cpml
        self.domain.nz = self.domain.nz + 2*cpml
        if self.domain.dim == 2.5:
            self.domain.ny = self.domain.ny + 2*cpml
        
        xyz = xyz.to_numpy() 
        
        if xyz.shape[1] == 1:
            xyz = xyz.T
        
        # We want to make sure the shape of the self.receiver_xyz array is in 
        # the correct shape. We won't be able to differentiate if the array is 
        # in correct shape if it is 3x3
        if xyz.shape[0] == 3 and np.prod(xyz.shape) > 9:
            xyz = xyz.T
        if xyz.shape[0] == 3 and np.prod(xyz.shape) == 6:
            xyz = xyz.T
        
        # If the receiver file contains values that are relative to the 
        # cartesian space of the domain, we want to change them to the indices
        # of the 
        if not self.receiver_indices:
            xyz = xyz / \
                np.array(
                    [
                        float(self.domain.dx), 
                        float(self.domain.dy), 
                        float(self.domain.dz) 
                    ]
                )
            xyz.astype(int)
        
        self.receiver_xyz = xyz + cpml
    
    # -------------------------------------------------------------------------
    def getrcx(self):
        """
        Loads the time series data for all receivers and handles complex data
        conditions.
        """
        # input rcx as an n-by-2 array integer values for their indices.
        src_ind = (
            self.source / \
            np.array([self.domain.dx, self.domain.dy, self.domain.dz])
        ).astype(int)
        if self.domain.dim == 2.5:
            all_files = glob(
                self.channel + '*.' + '.'.join(src_ind.astype(str)) + '.dat'
            )
        else: 
            all_files = glob(
                self.channel + '*.' + '.'.join(src_ind[np.array([0,2])].astype(str)) + '..dat'
            )
        
        all_files.sort()
        m = len(all_files)
        if len(self.receiver_xyz.shape) == 1:
            n = 1
            timeseries = np.zeros([m]) 
        else:
            n = len(self.receiver_xyz[:,0])
            timeseries = np.zeros([m,n])   
            
        timeseries_complex = timeseries.copy()
        
        if self.domain.dim == 2.0:
            self.domain.ny = None
        
        if self.domain.dim == 2.5:
            for i in range(m):
                npdat = read_dat(
                    all_files[i], 
                    self.channel, 
                    self.domain, 
                    self.is_complex, 
                    single = self.single_precision
                )
                if n == 1:
                    timeseries[:,i] = npdat.real[
                        int(self.receiver_xyz[2]), 
                        int(self.receiver_xyz[1]), 
                        int(self.receiver_xyz[0])
                    ]
                    timeseries_complex[:,i] = npdat.imag[
                        int(self.receiver_xyz[2]), 
                        int(self.receiver_xyz[1]), 
                        int(self.receiver_xyz[0])
                    ]
                else:
                    for j in range(0, n):
                        # Don't forget x is columns and z is rows
                        timeseries[i,j] = npdat.real[
                            int(self.receiver_xyz[j,2]),
                            int(self.receiver_xyz[j,1]),
                            int(self.receiver_xyz[j,0])
                        ]
                        timeseries_complex[i,j] = npdat.imag[
                            int(self.receiver_xyz[j,2]),
                            int(self.receiver_xyz[j,1]),
                            int(self.receiver_xyz[j,0])
                        ]
        else:
            for i in range(m):
                npdat = read_dat(
                    all_files[i], 
                    self.channel, 
                    self.domain, 
                    self.is_complex,
                    single = self.single_precision
                )
                if n == 1:
                    timeseries[i] = npdat.real[
                        int(self.receiver_xyz[2]), int(self.receiver_xyz[0])
                    ]
                    timeseries_complex[i] = npdat.imag[
                        int(self.receiver_xyz[2]), 
                        int(self.receiver_xyz[0])
                    ]
                else:
                    for j in range(n):
                        # Don't forget x is columns and z is rows
                        timeseries[i,j] = npdat.real[
                            int(self.receiver_xyz[j,2]),
                            int(self.receiver_xyz[j,0])
                        ]
                        timeseries_complex[i,j] = npdat.imag[
                            int(self.receiver_xyz[j,2]),
                            int(self.receiver_xyz[j,0])
                        ]
        
        # Store both of the time series in the array object
        self.timeseries_complex = timeseries_complex
        self.timeseries = timeseries
    
    # -------------------------------------------------------------------------
    def sectionplot(
            self, 
            plot_complex: bool = False
        ):
        """
        Creates a grayscale section plot of the time series data.
        
        :param plot_complex: Plot the complex part of the solution if True.
        :type plot_complex: bool
        """
        if plot_complex:
            # Use complex values 
            dat = self.timeseries_complex 
        else:
            dat = self.timeseries
        
        m,n = dat.shape
        
        if self.is_seismic:
            mult = 1e2
        else:
            mult = 1e6
        
        timelocs = np.arange(0, m, int(m/10) ) # 10 tick marks along y-axis
        rcxlocs = np.arange(0, np.max([n, 5]), int(np.max([n, 5])/5) ) # 5 tick marks along x-axis
        
        if self.is_seismic:
            timevals = np.round(timelocs*float(self.dt) * mult, 2)
        else:
            timevals = np.round(timelocs*float(self.dt) * mult, 2)
        
        if self.gain == 0:
            self.gain = 1
        
        if self.gain < m:
            for j in range(0, n):
                # Subtract the mean value
                # dat[:,j] = dat[:,j] - np.mean(dat[:,j])
                dat[:,j] = agc(dat[:,j], self.gain, "mean")
        
        self.fig = plt.figure()#figsize =(n/2,m/2) )
        self.ax = plt.gca()
        
        self.ax.imshow(dat, cmap = 'Greys', aspect = 'auto')
        self.ax.set_xlabel(r'Receiver #')
        self.ax.xaxis.tick_top()
        self.ax.xaxis.set_label_position('top')
        self.ax.set_xticks(rcxlocs)
        self.ax.set_xticklabels(rcxlocs)
        self.ax.set_ylabel(r'Two-way Travel Time (s)')
        self.ax.set_yticks(timelocs)
        self.ax.set_yticklabels(timevals)
        
        # Other figure handle operations
        self.ax.set_aspect(aspect = self.exaggeration)
        
        if self.is_seismic:
            self.ax.text(0, m + 0.03*m, 'x $10^{-2}$')
        else:
            self.ax.text(0, m + 0.03*m, 'x $10^{-6}$')
        
        self.ax.update_datalim( ((0,0),(m, n)))
        plt.show()
    
    # -------------------------------------------------------------------------
    def wiggleplot(
            self, 
            receiver_number: int, 
            plot_complex: bool = False, 
            plot_background: str = 'none',
            figure_size: Tuple[float, float] = (8, 5),
            **kwargs
        ):
        '''
        Plot an individual receiver from the list of the receiver files. Use 
        indices in the Python reference frame with the first component being 0. 
        
        :param receiver_number: Specify the indice of the reciever from the 
            receiver file to plot.  
        :type receiver_number: int
        :param plot_complex: Plot the complex part of the solution if True.
        :type plot_complex: bool
        :param plot_background: Specify the plot background color. Default is 
            transparent.
        :type plot_background: str
        :param figure_size: Specify the figure dimensions in inches. Default is
            (8,5) width and height, respectively.
        :type figure_size: Tuple[float, float]
        :param **kwargs: Additional plotting parameters as defined in 
            matplotlib.pyplot.plot. 
        :type **kwargs: dict
        '''
        default_plotspec = {
            'linewidth': 2,
            'linestyle': '-',
            'color': 'k',
            'alpha': 1
        }
        
        plot_params = {**default_plotspec, **kwargs}
        
        if plot_complex:
            # Use complex values 
            dat = self.timeseries_complex[:,n]
        else:
            dat = self.timeseries[:,n]
            
        timevector = np.arange(0, len(dat) ) * self.dt 
        
        self.wigglefig = plt.figure(facecolor='none', figsize = figure_size )
        self.wiggleax = plt.gca()
        
        self.wiggleax.plot(timevector, dat, **plot_params)
        self.wiggleax.set_facecolor(plot_background)
        ax.set_xlabel('Two-way travel time (s)')
        if self.is_seismic:
            self.wiggleax.set_ylabel('Velocity (m/s)')
        else:
            self.wiggleax.set_ylabel('Electric Field (V/m)')
        
        plt.tight_layout()
        
    
    # -------------------------------------------------------------------------
    def save(self, as_csv=False):
        """
        Save the object as a pickle formatted file or the numpy array of 
        receiver time series to a CSV.
        
        :param as_csv: Flag if you would like to export the time series to CSV. 
            This creates a filename formatted as 
            <project_name>.<channel>.<src_x>.<src_y>.<src_z>.csv. DEFAULT is False. 
        :type as_csv: bool
        """
        filename = '.'.join( 
            [
                self.prjfile.split('.')[:-1][0],
                self.channel, 
                '.'.join(self.source.astype(str))
            ]
        )
        if as_csv:
            filename = filename + '.csv'
            df = pd.DataFrame(self.timeseries)
            df.to_csv(filename, header = False, index = False)
        else:
            # Pickle the object and save to file
            filename = filename + '.pkl'
            with open(filename, 'wb') as file:
                pickle.dump(self, filename)
            
        
        


# =============================================================================
# ============================== Main Function ================================
def main(
        prjfile: str, 
        receiver_file: str, 
        channel: str, 
        rind: bool, 
        is_complex: bool, 
        single_precision: bool, 
        exaggeration: float, 
        gain: int,
        plot_complex: bool = False,
        plot: bool = False
    ) -> None:
    """
    The main function creating an Array object and triggering plotting or 
    saving actions based on user inputs.
    
    :param prjfile: Path to the project file.
    :type prjfile: str
    :param receiver_file: Path to the receiver file.
    :type receiver_file: str
    :param channel: Data channel to be analyzed.
    :type channel: str
    :param rind: Flag for coordinate indices.
    :type rind: bool
    :param is_complex: Flag for complex data handling.
    :type is_complex: bool
    :param single_precision: Flag for single precision data.
    :type single_precision: bool
    :param exaggeration: Aspect ratio between the x and y axes for plotting.
    :type exaggeration: float
    :param gain: Smoothing length for the data.
    :type gain: int
    :param plot_complex: Flag to plot complex part of the solution.
    :type plot_complex: bool
    :param plot: Flag to enable plotting.
    :type plot: bool
    """
    array = Array(
        channel,
        prjfile, 
        receiver_file,
        rind, 
        single_precision = single_precision,
        is_complex = is_complex
    )
    
    if gain:
        array.gain = gain
    if exaggeration:
        array.exaggeration = exaggeration
    if save:
        array.save()
    if plot:
        array.section_plot(plot_complex = plot_complex)



# =============================================================================
# ========================== Command Line Arguments ===========================
# =============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This program creates an Array class that can be saved in 
        pickle format. Receiver locations listed in the the specified receiver 
        file are extracted to create the time series. The minimum number of 
        plots is 1 which will return a wiggle plot.""" 
    )
    
    parser.add_argument(
        '-p', '--prjfile',
        nargs = 1, type = str, required = True,
        help = 'The project file path'
    )
    
    parser.add_argument(
        '-r', '--rcxfile',
        nargs=1, type=str, required = True,
        help='the file path for the text file of receiver locations'
    )
    
    parser.add_argument(
        '-i', '--index',
        action = 'store_true', required = False,
        help = """Indicate whether the receiver file contains coordinate indices 
        or if these are the locations in meters. Default (0 - meters)"""
    )
    
    parser.add_argument(
        '-c', '--channel',
        nargs = 1, type = str, required = True,
        help = """The channel to query. """
    )
    
    parser.add_argument(
        '-z', '--is_complex', action = 'store_true', required = False,
        help = """Flag whether the data will be complex valued. If the data is
        not flagged but is complex, only the real data will be returned. """
    )
    
    parser.add_argument(
        '-d', 'double_precision', action = 'store_false', required = False,
        help = '''Flag whether the model outputs are in double precision. 
        Default is single precision. '''
    )
    
    parser.add_argument(
        '-g', '--gain',
        nargs = 1, type = float, required = False, default = [None],
        help = "The smoothing length"
    )
    
    parser.add_argument(
        '-e', '--exaggeration',
        nargs=1, type = float, required = False, default = [None],
        help = """Set the aspect ratio between the x and y axes for
        plotting."""
    )
    
    parser.add_argument(
        '-s', '--save', action = 'store_true', required = False, 
        help = """Flag to save the array object as a .pkl file."""
    )
    
    parser.add_argument(
        '-P', '--plot', action = 'store_true', required = False,
        help = """Flag whether you want to return a plot to the console."""
    )
    
    parser.add_argument(
        '-C', '--plot_complex', action = 'store_true', required = False,
        help = """
        If flagged, plot the complex part of the solution otherwise default to 
        the real valued solution. 
        """
    )
    
    # Get the arguments
    args = parser.parse_args()
    prjfile = ''.join(args.prjfile)
    receiver_file = ''.join(args.rcxfile)
    channel = ''.join(args.channel)
    rind = args.index
    is_complex = args.is_complex
    single_precision = args.double_precision 
    exaggeration = args.exaggeration[0] 
    gain = args.gain[0]
    plot_complex = args.plot_complex
    plot = args.plot 
    
    # ==================== Create the object and assign inputs ====================    
    main(
        prjfile, 
        receiver_file, 
        channel, 
        rind, 
        is_complex, 
        single_precision, 
        exaggeration, 
        gain, 
        plot_complex = plot_complex,
        plot = plot
    )
        
