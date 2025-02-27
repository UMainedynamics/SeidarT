import numpy as np
import pandas as pd
from glob2 import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from typing import Tuple
from seidart.routines.definitions import *
import dill as pickle
from scipy.fft import fft2, fftshift
from scipy.signal import spectrogram, windows
from matplotlib.colors import TwoSlopeNorm
from seidart.routines.classes import Domain, Material, Model

# =============================================================================
class Array:
    def __init__(
            self, 
            channel: str,
            project_file: str, 
            receiver_file: str,
            single_precision: bool = True,
            csvfile: str = None, 
        ):
        """
        Initializes the Array object with project settings and receiver details.
        This will build a set of time series for each specified receiver from 
        the receiver file, or alternatively, a previously built CSV file can be 
        loaded. The save function produces both a CSV and Pickle (.pkl) file. 
        Loading the pickle file will restore the object and all of the variables
        within it. 
        
        :param channel: The data channel to be analyzed.
        :type channel: str
        :param project_file: Path to the project file.
        :type project_file: str
        :param receiver_file: Path to the receiver locations file.
        :type receiver_file: str
        :param receiver_indices: Flag if receiver file contains indices (True) 
                                 or Cartesian coordinates (False).
        :type receiver_indices: bool
        :param single_precision: Use single precision for numerical data 
            if True.
        :type single_precision: bool
        :param csvfile: The CSV file containing the set of receiver time series.
        :type csvfile: str
        """
        self.project_file = project_file 
        self.channel = channel
        self.receiver_file = receiver_file
        self.receiver_indices = None
        self.receiver_locations = None # The locations in meters
        self.single_precision = single_precision
        self.source_xyz = None 
        self.stream = None 
        self.agc_gain_window = None
        self.exaggeration = 0.5
        self.geometric_correction_exponent = 1.0 # Geometric spreading correction parameter
        self.alpha_exponent = 1.0 # Exponential gain parameter
        self.distances = None
        self.csvfile = csvfile
        self.csvfilename = None
        self.pklfilename = None
        self.simulation_type = 'CSG'
        self.build()
    
    # -------------------------------------------------------------------------
    def build(self):
        """
        Constructs domain and models based on project file and sets up source 
        and receiver configurations.
        """
        self.domain, self.material, self.seismic, self.electromag = loadproject(
            self.project_file,
            Domain(), 
            Material(),
            Model(),
            Model()
        )
        if self.channel in ['Vx','Vy','Vz']:
            self.is_seismic = True
            self.dt = self.seismic.dt 
            self.agc_gain_window = int(self.seismic.time_steps)
        else:
            self.is_seismic = False
            self.dt = self.electromag.dt
            self.agc_gain_window = int(self.electromag.time_steps)
        
        
        # Load the receiver file as a numpy array
        self.loadxyz()
        # Add the source locations into their respective variables. 
        if self.channel == 'Vx' or self.channel == 'Vy' or self.channel == 'Vz':
            if self.domain.dim == '2.0':
                self.source_xyz = np.array(
                    [self.seismic.x, self.seismic.z]
                )
            else:
                self.source_xyz = np.array(
                    [
                        self.seismic.x, 
                        self.seismic.y, 
                        self.seismic.z
                    ]
                )
        else:
            if self.domain.dim == '2.0':
                self.source_xyz = np.array(
                    [self.electromag.x, self.electromag.z]
                )
            else:
                self.source_xyz = np.array(
                    [
                        self.electromag.x, 
                        self.electromag.y, 
                        self.electromag.z
                    ]
                )
        
        if self.csvfile:
            self.timeseries = pd.read_csv(
                self.csvfile, header = None
            ).to_numpy()
        else:
            # Load the time series for all receivers
            self.getrcx()
        
    # -------------------------------------------------------------------------
    def loadxyz(self):
        """
        Loads and sorts receiver locations from a file and adjusts them 
        according to the domain and CPML layer.  If the source_file flag is 
        True, the source file is loaded instead of the receiver file and saved
        to the object under source.

        :param source_file: Flag to indicate if the file is a source file.
        :type source_file: bool 

        """
        
        
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
        
        # Load the receiver file
        xyz = pd.read_csv(self.receiver_file)
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
        
        self.receiver_indices = (xyz / \
            np.array(
                [
                    float(self.domain.dx), 
                    float(self.domain.dy), 
                    float(self.domain.dz) 
                ]
            )).astype(int) + cpml
        self.receiver_xyz = xyz.copy()
        # The following line is to make this compatible with the CommonOffset 
        # object setup. 
        self.receiver_xyz_all = self.receiver_xyz.copy()
    
    # -------------------------------------------------------------------------
    def getrcx(self):
        """
        Loads the time series data for all receivers
        """
        # input rcx as an n-by-2 array integer values for their indices.
        src_ind = (
            self.source_xyz / \
            np.array([self.domain.dx, self.domain.dy, self.domain.dz])
        ).astype(int)
        if self.domain.dim == 2.5:
            all_files = glob(
                self.channel + '*.' + '.'.join(src_ind.astype(str)) + '.dat'
            )
        else: 
            all_files = glob(
                self.channel + '*.' + '.'.join(src_ind[np.array([0,2])].astype(str)) + '.dat'
            )
        all_files.sort()
        m = len(all_files)
        if len(self.receiver_indices) == 1:
            n = 1
            timeseries = np.zeros([m])
            # Sometimes the 1 row array needs to be
            if self.receiver_indices.ndim == 2:
                self.receiver_indices = self.receiver_indices[0]
        else:
            n = len(self.receiver_indices[:,0])
            timeseries = np.zeros([m,n])   
                    
        if self.domain.dim == 2.0:
            self.domain.ny = None
        
        if self.domain.dim == 2.5:
            for i in range(m):
                npdat = read_dat(
                    all_files[i], 
                    self.channel, 
                    self.domain, 
                    single = self.single_precision
                )
                if n == 1:
                    timeseries[:,i] = npdat.real[
                        int(self.receiver_indices[2]), 
                        int(self.receiver_indices[1]), 
                        int(self.receiver_indices[0])
                    ]
                else:
                    for j in range(0, n):
                        # Don't forget x is columns and z is rows
                        timeseries[i,j] = npdat.real[
                            int(self.receiver_indices[j,2]),
                            int(self.receiver_indices[j,1]),
                            int(self.receiver_indices[j,0])
                        ]
        else:
            for i in range(m):
                npdat = read_dat(
                    all_files[i], 
                    self.channel, 
                    self.domain, 
                    single = self.single_precision
                )
                if n == 1:
                    timeseries[i] = npdat.real[
                        int(self.receiver_indices[2]), int(self.receiver_indices[0])
                    ]
                else:
                    for j in range(n):
                        # Don't forget x is columns and z is rows
                        timeseries[i,j] = npdat.real[
                            int(self.receiver_indices[j,2]),
                            int(self.receiver_indices[j,0])
                        ]
        
        # Store both of the time series in the array object
        self.timeseries = timeseries
        
    # -------------------------------------------------------------------------
    def srcrcx_distance(self):
        """
        Compute the source and receiver(s) distance 
        """
        if self.source_xyz.ndim == 1:
            source_xyz = self.source_xyz.reshape(1,3)
        else:
            source_xyz = self.source_xyz.copy()
        
        if self.domain.dim == '2.0':
            self.distances = (source_xyz[:,(0,2)]  - self.receiver_xyz[:,(0,2)]  )**2
        else: 
            self.distances = (source_xyz - self.receiver_xyz)**2
        
        self.distances = np.sqrt(np.sum(self.distances,axis = 1))
    
    # -------------------------------------------------------------------------
    def amplitude_correction(self, dat: np.ndarray, correction_type: str = None):
        """
        
        :param amplitude_correction_type: Correct the amplitudes for geometric spreading ('GS')or 
            using an auto gain control function ('AGC'); Default is None
        :type amplitude_correction_type: str
        """
        # ------------------------------
        m, n = dat.shape
        # Apply any amplitude corrections if specified
        if correction_type == 'AGC' and self.agc_gain_window < m:
            if self.agc_gain_window == 0:
                self.agc_gain_window = 1
            
            for j in range(0, n):
                # Subtract the mean value
                # dat[:,j] = dat[:,j] - np.mean(dat[:,j])
                dat[:,j] = agc(dat[:,j], self.agc_gain_window, "mean")
            
            self.agc_corrected_timeseries = dat.copy() 
        # ------------------------------
        if correction_type == 'GS':
            dat = correct_geometric_spreading(
                dat, self.distances, self.geometric_correction_exponent
            )
            self.gs_corrected_timeseries = dat.copy()

        if correction_type == 'EXP':
            for j in range(0, n):
                dat[:,j] = exponential_gain(dat[:,j], self.alpha_exponent)
            self.exp_corrected_timeseries = dat.copy() 
        
        return dat

    # -------------------------------------------------------------------------
    def alpha_attenuation(self, phase_velocity, direction = 'z', material_indice = 0, fc = None):
            # Get the attenuation coefficient 
            gamma = self.seismic.attenuation_coefficients[f'gamma_{direction}'][material_indice]
            
            if not fc:
                fc = self.seis.f0
                
            # Compute the phase velocity for the given direction
            if Q:
                alpha = 2*np.pi * fc / (2*Q*phase_velocity)
            else:
                alpha = np.pi * fc * gamma / phase_velocity 
            
            return alpha

    # -------------------------------------------------------------------------
    def amplitude_calculation(self, distance, alpha, A0 = 1, spreading = "spherical"):
            dA = np.exp(-alpha * distance)
            
            if spreading == "spherical":
                dG = 1 / distance if distance > 0 else 1
            elif spreading == "cylindrical":
                dG = 1 / np.sqrt(ditance) if distance > 0 else 1 
            else:
                dG = 1 
            
            At = A0 * dA *dG
            return At, dA, dG 

    # -------------------------------------------------------------------------
    def sectionplot(
            self, 
            amplitude_correction_type: str = None,
            colormap: str = 'Greys',
            figure_size: Tuple[float, float] = (8,8),
        ):
        """
        Creates a grayscale section plot of the time series data.
        
        :param amplitude_correction_type: Correct the amplitudes for geometric 
            spreading ('GS') or using an auto gain control function ('AGC'); 
            Default is None.
        :type amplitude_correction_type: str
        :param colormap: Specify the colormap to be used in the plot. Default is grey.
        :type colormap: str 
        :param figure_size: Specify the figure size.
        :type figure_size: Tuple 
        """

        dat = self.timeseries.copy()
        
        m,n = dat.shape
        time_max = self.dt * m 
        
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
        
        # Amplitude correction
        if amplitude_correction_type:
            dat = self.amplitude_correction(
                dat, correction_type = amplitude_correction_type
            )
        
        norm = TwoSlopeNorm(vmin=-np.max(dat), vmax=np.max(dat), vcenter=0)

        self.fig_section = plt.figure(figsize = figure_size )
        self.ax_section = plt.gca()
        
        self.ax_section.imshow(
            dat, cmap = colormap, aspect = 'auto', norm = norm,
            extent = [0, n, m, 0]
        )
        self.ax_section.set_xlabel(r'Receiver #')
        self.ax_section.xaxis.tick_top()
        self.ax_section.xaxis.set_label_position('top')
        self.ax_section.set_xticks(rcxlocs)
        self.ax_section.set_xticklabels(rcxlocs)
        self.ax_section.set_ylabel(r'Two-way Travel Time (s)')
        self.ax_section.set_yticks(timelocs)
        self.ax_section.set_yticklabels(timevals)
        
        # Other figure handle operations
        self.ax_section.set_aspect(aspect = self.exaggeration)
        
        if self.is_seismic:
            self.ax_section.text(0, m + 0.03*m, 'x $10^{-2}$')
        else:
            self.ax_section.text(0, m + 0.03*m, 'x $10^{-6}$')
        
        # self.ax_section.update_datalim( ((0,0),(m, n)))
        plt.show()
    
    # -------------------------------------------------------------------------
    def wiggleplot(
            self, 
            receiver_number: int, 
            plot_background: str = 'none',
            plot_vertical = True,
            positive_fill_color = None,
            negative_fill_color = None,
            figure_size: Tuple[float, float] = (8, 5),
            **kwargs
        ):
        """
        Plot an individual receiver from the list of the receiver files. Use 
        indices in the Python reference frame with the first component being 0. 
        
        :param receiver_number: Specify the indice of the reciever from the 
            receiver file to plot.  
        :type receiver_number: int
        :param plot_background: Specify the plot background color. Default is 
            transparent.
        :type plot_background: str
        :param figure_size: Specify the figure dimensions in inches. Default is
            (8,5) width and height, respectively.
        :type figure_size: Tuple[float, float]
        :param **kwargs: Additional plotting parameters as defined in 
            matplotlib.pyplot.plot. 
        :type **kwargs: dict
        
        """
        
        default_plotspec = {
            'linewidth': 2,
            'linestyle': '-',
            'color': 'k',
            'alpha': 1
        }
        
        plot_params = {**default_plotspec, **kwargs}
        
        dat = self.timeseries[:,receiver_number]
            
        timevector = np.arange(0, len(dat) ) * self.dt 
        
        self.fig_wiggle = plt.figure(facecolor='none', figsize = figure_size )
        self.ax_wiggle = plt.gca()
        
        if plot_vertical:
            self.ax_wiggle.plot(dat, timevector, **plot_params)
            self.ax_wiggle.set_ylabel('Two-way travel time (s)')
            if self.is_seismic:
                self.ax_wiggle.set_xlabel('Velocity (m/s)')
            else:
                self.ax_wiggle.set_xlabel('Electric Field (V/m)')
            # self.ax_wiggle.set_ylim(self.ax_wiggle.get_ylim()[::-1])
            self.ax_wiggle.invert_yaxis() 
        else:
            self.ax_wiggle.plot(timevector, dat, **plot_params)
            self.ax_wiggle.set_xlabel('Two-way travel time (s)')
            if self.is_seismic:
                self.ax_wiggle.set_ylabel('Velocity (m/s)')
            else:
                self.ax_wiggle.set_ylabel('Electric Field (V/m)')
            
        self.ax_wiggle.set_facecolor(plot_background)
            
        plt.tight_layout()
        plt.show()
    
    # -------------------------------------------------------------------------
    def sectionwiggleplot(
            self, 
            receiver_indices: np.ndarray,
            scaling_factor: float = 1.0,
            receiver_distance: np.ndarray = None,
            plot_vertical = False,
            positive_fill_color = None,
            negative_fill_color = None,
            figure_size: Tuple[float, float] = (5,8),
            amplitude_correction_type = None,
            yaxis_label: str = 'Source-Receiver Distance (m)',
            gridded = True
        ):
        """
        Plots a seismic or electromagnetic record section using wiggle traces. 
        This function generates a wiggle plot for the selected receiver traces, 
        with options for vertical plotting, scaling, and amplitude corrections. 
        The positive and negative amplitudes can be filled with custom colors.

        :param receiver_indices: An array of indices specifying which receiver 
            traces to plot.
        :type receiver_indices: np.ndarray
        :param scaling_factor: A scaling factor for adjusting the amplitude of 
            the traces (default is 1.0).
        :type scaling_factor: float, optional
        :param receiver_distance: An array specifying the horizontal 
            (or vertical) positions of the receiver traces. If None, a 
            sequential index is used (default is None).
        :type receiver_distance: np.ndarray, optional
        :param plot_vertical: If True, plots the traces vertically rather than 
            horizontally (default is False).
        :type plot_vertical: bool, optional
        :param positive_fill_color: Color to fill the positive amplitudes of the 
            traces. Use any valid matplotlib color (default is None).
        :type positive_fill_color: str, optional
        :param negative_fill_color: Color to fill the negative amplitudes of the 
            traces. Use any valid matplotlib color (default is None).
        :type negative_fill_color: str, optional
        :param figure_size: Size of the figure in inches (default is (5, 8)).
        :type figure_size: Tuple[float, float], optional
        :param amplitude_correction_type: Type of amplitude correction to apply 
            before plotting. This should be a valid correction type recognized 
            by the `amplitude_correction` method (default is None).
        :type amplitude_correction_type: str, optional
        :param yaxis_label: The label corresponding to the y-axis if not using 
            source-receiver distance.
        :type yaxis_label: str
        :param gridded: Specify if a gridded plot should be created. Grids will 
            be drawn in vertically for a horizontally oriented graph (default) 
            and drawn horizontally for a vertically oriented graph. 

        :returns: None
        :rtype: None

        :notes:
            - Vertical plotting hasn't yet been implemented 
            - Receiver distances is generically named and is merely a time shift
              in the time series so that they are layed out as a section plot. 
              The traces can be sorted and plotted according to other inputs.  
        """
        dat = self.timeseries.copy()
        m,n = dat.shape
        n_traces_to_plot = len(receiver_indices) 
        time = np.arange(m)*self.dt 

        if amplitude_correction_type:
            dat = self.amplitude_correction(
                dat, correction_type = amplitude_correction_type
            )

        fig_wiggles, ax_wiggles = plt.subplots(figsize = figure_size)

        if plot_vertical:
            for i in range(n_traces_to_plot):
                indice = receiver_indices[i]
                ax_wiggles.plot(
                    dat[:, indice] * scaling_factor + receiver_distance[indice], 
                    time,
                    'k', linewidth=0.8
                )
                if positive_fill_color:
                    ax_wiggles.fill_betweenx(
                        time, receiver_distance[indice],
                        dat[:, indice] * scaling_factor + receiver_distance[indice],
                        where=(dat[:, indice] > 0), color=positive_fill_color
                    )
                if negative_fill_color:
                    ax_wiggles.fill_betweenx(
                        time, receiver_distance[indice],
                        dat[:, indice] * scaling_factor + receiver_distance[indice],
                        where=(dat[:, indice] < 0), color=negative_fill_color
                    )

            ax_wiggles.invert_yaxis()
            ax_wiggles.set_ylabel('Time (s)')
            ax_wiggles.set_xlabel(yaxis_label)
        else:
            for i in range(n_traces_to_plot):
                indice = receiver_indices[i]
                ax_wiggles.plot(
                    time, dat[:,indice]*scaling_factor + receiver_distance[indice],
                    'k', linewidth = 0.8
                )
                if positive_fill_color:
                    ax_wiggles.fill_between(
                        time, receiver_distance[indice],
                        dat[:,indice] * scaling_factor + receiver_distance[indice],
                        where=(dat[:,indice] > 0), color = positive_fill_color
                    )
                if negative_fill_color:
                    ax_wiggles.fill_between(
                        time, receiver_distance[indice],
                        dat[:,indice] * scaling_factor + receiver_distance[indice],
                        where=(dat[:,indice] < 0), color = positive_fill_color
                    )

            ax_wiggles.invert_yaxis()
            ax_wiggles.set_xlabel('Time (s)')
            ax_wiggles.set_ylabel(yaxis_label)
            if gridded:
                ax_wiggles.grid(
                    True, which = 'both', axis = 'x', 
                    linestyle = '--', color = '#c7c7c7', linewidth = 0.7
                )
        
        plt.show() 

    # -------------------------------------------------------------------------
    def fk_analysis(
            self, 
            d_rcx: float, 
            figure_size: Tuple[float,float] = (10,8), 
            colormap: str = 'viridis',
            contour_levels: int = 100,
            frequency_limits: Tuple[float,float] = None,
            wavenumber_limits: Tuple[float,float] = None,
        ):
        """
        Compute and plot the frequency-wavenumber analysis for an array with 
        regularly spaced intervals.
        
        :param d_rcx: The receiver interval
        :type d_rcx: float 
        
        """
        m, n = self.timeseries.shape 
        fk_spectrum = fftshift(fft2(self.timeseries))
        
        # Create the frequency and wavenumber vectors then we need to shift the 
        # frequencies and wavenumbers for plotting 
        freqs = fftshift(np.fft.fftfreq(m, self.dt))
        wavenumbers = -fftshift(np.fft.fftfreq(n, d_rcx)) # positive down coordinates means we have to flip the f-k 
        
        
        self.fig_fk, self.ax_fk = plt.subplots(figsize = figure_size)
        contour = self.ax_fk.contourf(
            wavenumbers, freqs, np.abs(fk_spectrum), 
            cmap = colormap, levels = contour_levels,
        )
        cbar = self.fig_fk.colorbar(contour, ax = self.ax_fk, label = 'Amplitude')
        self.ax_fk.set_xlabel(r'Wavenumber (m$^{-1}$)')
        self.ax_fk.set_ylabel('Frequency (Hz)')
        if frequency_limits:
            self.ax_fk.set_ylim(frequency_limits)
        
        if wavenumber_limits:
            self.ax_fk.set_xlim(wavenumber_limits)
        
        plt.show()
    
    # -------------------------------------------------------------------------
    def dispersion_analysis(
            self,
            d_rcx: float,
            figure_size: Tuple[float, float] = (10, 8),
            colormap: str = 'viridis',
            contour_levels: int = 100,
            frequency_limits: Tuple[float, float] = None,
            velocity_limits: Tuple[float, float] = None,
            time_bandwidth: float = 4.0,  # Time-bandwidth product for DPSS
            n_tapers: int = 5,  # Number of tapers
        ):
        """
        Compute and plot the dispersion image using multitaper spectral analysis
        with Scipy's DPSS.

        :param d_rcx: The receiver interval (spacing between receivers).
        :type d_rcx: float
        :param time_bandwidth: Time-bandwidth product for DPSS tapers.
        :type time_bandwidth: float
        :param n_tapers: Number of DPSS tapers to use.
        :type n_tapers: int
        """
        # Step 1: Preprocess data (multitaper applied per time-series)
        m, n = self.timeseries.shape
        fs = 1.0 / self.dt  # Sampling frequency
        multitaper_spectrum = np.zeros((m, n), dtype=np.complex128)

        for i in range(n):
            # Generate DPSS tapers
            tapers = windows.dpss(m, NW=time_bandwidth, Kmax=n_tapers)
            spectra = []
            for taper in tapers:
                # Compute the FFT of the tapered data
                tapered_data = self.timeseries[:, i] * taper
                spectrum = np.fft.fft(tapered_data)
                spectra.append(spectrum)
            # Average the multitaper spectra
            multitaper_spectrum[:, i] = np.mean(spectra, axis=0)

        # Step 2: Perform 2D FFT for f-k spectrum
        fk_spectrum = fftshift(fft2(multitaper_spectrum))
        freqs = fftshift(np.fft.fftfreq(m, self.dt))
        wavenumbers = -fftshift(np.fft.fftfreq(n, d_rcx))  # Flip for positive down

        # Step 3: Compute phase velocity (c = f / k)
        fk_mesh_freq, fk_mesh_k = np.meshgrid(freqs, wavenumbers, indexing='ij')
        with np.errstate(divide='ignore', invalid='ignore'):
            phase_velocity = np.abs(fk_mesh_freq / fk_mesh_k)

        # Step 4: Bin amplitude into f-c space
        if velocity_limits is None:
            velocity_limits = (0, np.max(phase_velocity[np.isfinite(phase_velocity)]))
        velocity_bins = np.linspace(*velocity_limits, 500)

        amplitude_fc = np.zeros((len(freqs), len(velocity_bins)))
        for i, f_row in enumerate(np.abs(fk_spectrum)):
            c_indices = np.digitize(phase_velocity[i], velocity_bins)
            for j, c_idx in enumerate(c_indices):
                if 0 <= c_idx < len(velocity_bins):
                    amplitude_fc[i, c_idx] += f_row[j]

        # Step 5: Plot the dispersion image
        self.fig_dispersion, self.ax_dispersion = plt.subplots(figsize=figure_size)
        contour = self.ax_dispersion.contourf(
            velocity_bins, freqs, amplitude_fc, cmap=colormap, levels=contour_levels
        )
        cbar = self.fig_dispersion.colorbar(contour, ax=self.ax_dispersion, label='Amplitude')
        self.ax_dispersion.set_xlabel('Phase Velocity (m/s)')
        self.ax_dispersion.set_ylabel('Frequency (Hz)')
        if frequency_limits:
            self.ax_dispersion.set_ylim(frequency_limits)
        if velocity_limits:
            self.ax_dispersion.set_xlim(velocity_limits)

        plt.show()

    # -------------------------------------------------------------------------
    def save(self, save_object = True, output_basefile = None):
        """
        Save the object as a pickle formatted file and the numpy array of 
        receiver time series to a CSV.
        
        """
        if output_basefile:
            filename = output_basefile
        else:
            filename = '-'.join( 
                [
                    self.project_file.split('.')[:-1][0],
                    self.channel, 
                    '-'.join(self.source_xyz.astype(str))
                ]
            )
        
        self.csvfilename = filename + '.csv'
        self.pklfilename = filename + '.pkl'
        
        df = pd.DataFrame(self.timeseries)
        df.to_csv(self.csvfilename, header = False, index = False)
        
        # Pickle the object and save to file
        with open(self.pklfilename, 'wb') as file:
            pickle.dump(self, file)
        
# =============================================================================
# ============================== Main Function ================================
def main(
        project_file: str, 
        receiver_file: str, 
        channel: str, 
        rind: bool, 
        single_precision: bool, 
        exaggeration: float, 
        gain: int,
        plot: bool = False
    ) -> None:
    """
    The main function creating an Array object and triggering plotting or 
    saving actions based on user inputs.
    
    :param project_file: Path to the project file.
    :type project_file: str
    :param receiver_file: Path to the receiver file.
    :type receiver_file: str
    :param channel: Data channel to be analyzed.
    :type channel: str
    :param rind: Flag for coordinate indices.
    :type rind: bool
    :param single_precision: Flag for single precision data.
    :type single_precision: bool
    :param exaggeration: Aspect ratio between the x and y axes for plotting.
    :type exaggeration: float
    :param gain: Smoothing length for the data.
    :type gain: int
    :param plot: Flag to enable plotting.
    :type plot: bool
    """
    parser = argparse.ArgumentParser(
        description="""This program creates an Array class that can be saved in 
        pickle format. Receiver locations listed in the the specified receiver 
        file are extracted to create the time series. The minimum number of 
        plots is 1 which will return a wiggle plot.""" 
    )
    
    parser.add_argument(
        '-p', '--project_file',
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
        '-w', '--wiggleplot', required = False, nargs = 1, type = int, 
        help = '''Provide the index of the receiver in the receiver file that 
        you would like to plot as a single 1D time series. The channel provided 
        in the channel argument is the dimension that will be plotted.'''
    )
    
    # Get the arguments
    args = parser.parse_args()
    project_file = ''.join(args.project_file)
    receiver_file = ''.join(args.rcxfile)
    channel = ''.join(args.channel)
    rind = args.index
    single_precision = args.double_precision 
    exaggeration = args.exaggeration[0] 
    gain = args.gain[0]
    plot = args.plot 
    
    array = Array(
        channel,
        project_file, 
        receiver_file,
        rind, 
        single_precision = single_precision,
    )
    
    if gain:
        array.agc_gain_window = gain
    if exaggeration:
        array.exaggeration = exaggeration
    if save:
        array.save()
    if plot:
        array.section_plot()



# =============================================================================
# ========================== Command Line Arguments ===========================
# =============================================================================
if __name__ == "__main__":
    main()
    
    
    
