import numpy as np
import pandas as pd
from glob2 import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from typing import Tuple, Optional, Dict, List, Union
from seidart.routines.definitions import *
import dill as pickle
from scipy.ndimage import gaussian_filter
from numpy.fft import fft2, fftshift, rfft, rfftfreq
from scipy.fft import fft2, fftshift
from scipy.signal import spectrogram, windows, butter, filtfilt, sosfiltfilt
from scipy.signal.windows import tukey
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
                dat[:,j] = exponential_gain(dat[:,j], self.dt, self.alpha_exponent)
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
            use_filtered: bool = False
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

        if use_filtered:
            dat = self.timeseries_filtered.copy()
        else:
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
            use_filtered: bool = False,
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
        
        if use_filtered:
            dat = self.timeseries_filtered.copy()
        else:
            dat = self.timeseries.copy()
            
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
            plot_vertical = False,
            positive_fill_color = None,
            negative_fill_color = None,
            figure_size: Tuple[float, float] = (5,8),
            amplitude_correction_type = None,
            yaxis_label: str = 'Source-Receiver Distance (m)',
            gridded: bool = True, 
            use_filtered: bool = False
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
        
        if self.distances is None:
            self.srcrcx_distance() 
        
        if use_filtered:
            dat = self.timeseries_filtered.copy()
        else:
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
                    dat[:, indice] * scaling_factor + self.distances[indice], 
                    time,
                    'k', linewidth=0.8
                )
                if positive_fill_color:
                    ax_wiggles.fill_betweenx(
                        time, self.distances[indice],
                        dat[:, indice] * scaling_factor + self.distances[indice],
                        where=(dat[:, indice] > 0), color=positive_fill_color
                    )
                if negative_fill_color:
                    ax_wiggles.fill_betweenx(
                        time, self.distances[indice],
                        dat[:, indice] * scaling_factor + self.distances[indice],
                        where=(dat[:, indice] < 0), color=negative_fill_color
                    )

            ax_wiggles.invert_yaxis()
            ax_wiggles.set_ylabel('Time (s)')
            ax_wiggles.set_xlabel(yaxis_label)
        else:
            for i in range(n_traces_to_plot):
                indice = receiver_indices[i]
                ax_wiggles.plot(
                    time, dat[:,indice]*scaling_factor + self.distances[indice],
                    'k', linewidth = 0.8
                )
                if positive_fill_color:
                    ax_wiggles.fill_between(
                        time, self.distances[indice],
                        dat[:,indice] * scaling_factor + self.distances[indice],
                        where=(dat[:,indice] > 0), color = positive_fill_color
                    )
                if negative_fill_color:
                    ax_wiggles.fill_between(
                        time, self.distances[indice],
                        dat[:,indice] * scaling_factor + self.distances[indice],
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

    # --------------------------------------------------------------------------
    def butterworth_filter(self,
            filter_type: str,
            lowcut: float = None,
            highcut: float = None,
            order: int = 4,
            taper_duration: float = 0.0,
            pad_samples: int = 100,
        ):
        """
        Unified seismic filter with tapering and padding.
        
        Parameters
        ----------
        filter_type : str
            'lowpass', 'highpass', or 'bandpass'
        lowcut : float, optional
            Lower cutoff frequency in Hz (used for 'highpass' and 'bandpass')
        highcut : float, optional
            Upper cutoff frequency in Hz (used for 'lowpass' and 'bandpass')
        order : int
            Butterworth filter order
        taper_duration : float
            Duration in seconds for a pre-arrival taper (default = 0.0, no taper)
        pad_samples : int
            Number of samples to pad symmetrically at top and bottom
        
        Result
        ------
        Sets self.timeseries_filtered (same shape as self.timeseries)
        """
        
        data = self.timeseries
        fs = 1 / self.dt
        nyq = 0.5 * fs
        
        # Validate filter parameters
        if filter_type == 'lowpass':
            if not highcut or not (0 < highcut < nyq):
                raise ValueError(f"Invalid highcut: {highcut} Hz")
            Wn = highcut / nyq
            sos = butter(order, Wn, btype='low', output='sos')
        
        elif filter_type == 'highpass':
            if not lowcut or not (0 < lowcut < nyq):
                raise ValueError(f"Invalid lowcut: {lowcut} Hz")
            Wn = lowcut / nyq
            sos = butter(order, Wn, btype='high', output='sos')
        
        elif filter_type == 'bandpass':
            if not lowcut or not highcut or not (0 < lowcut < highcut < nyq):
                raise ValueError(f"Invalid bandpass range: {lowcut}–{highcut} Hz")
            Wn = [lowcut / nyq, highcut / nyq]
            sos = butter(order, Wn, btype='band', output='sos')
        
        else:
            raise ValueError(f"Unsupported filter type: '{filter_type}'")
        
        # Apply pre-arrival taper if requested
        if taper_duration > 0.0:
            nt = data.shape[0]
            taper_samples = int(taper_duration / self.dt)
            taper = np.ones(nt)
            ramp = 0.5 * (1 - np.cos(np.pi * np.arange(taper_samples) / taper_samples))
            taper[:taper_samples] = ramp
            data = data * taper[:, None]
        
        # Pad with zeros before filtering
        pad = np.zeros((pad_samples, data.shape[1]))
        padded = np.vstack([pad, data, pad])
        filtered = sosfiltfilt(sos, padded, axis=0)
        unpadded = filtered[pad_samples:-pad_samples, :]
        
        self.timeseries_filtered = unpadded
    
    # -------------------------------------------------------------------------
    def fk_analysis(
            self,
            d_rcx: float,
            ntfft: Optional[int] = None,
            nxfft: Optional[int] = None,
            taper: Union[str, bool] = 'time',
            figure_size: Tuple[float, float] = (10, 8),
            colormap: str = 'viridis',
            contour_levels: int = 100,
            frequency_limits: Tuple[float, float] = None,
            wavenumber_limits: Tuple[float, float] = None,
            wavenumber_units: str = 'cycles',
            use_filtered: bool = True,
            to_db: Union[str, bool] = 'normalized',
            mode_lines: Optional[dict] = None,
            mask_db: float = -40
        ):
        """
        Compute and plot the frequency–wavenumber (f-k) power spectrum of the time series
        for a regularly spaced seismic receiver array.

        Applies a 2D FFT to the array, with optional tapering, zero-padding, and
        amplitude-to-dB conversion. Overlays theoretical mode dispersion curves if provided.

        Parameters
        ----------
        d_rcx : float
            Receiver spacing in meters.
        ntfft : int, optional
            Number of FFT points in time. If None, defaults to the number of time samples.
        nxfft : int, optional
            Number of FFT points in space. If None, defaults to the number of receivers.
        taper : {'time', 'space', 'both', 'tukey_time', 'tukey_space', 'tukey_both', False}, default='time'
            Apply a Hanning taper in the specified direction(s) to suppress spectral leakage.
        figure_size : tuple of float, default=(10, 8)
            Size of the matplotlib figure in inches.
        colormap : str, default='viridis'
            Colormap used for plotting the contour.
        contour_levels : int, default=100
            Number of levels in the contour plot.
        frequency_limits : tuple of float, optional
            y-axis limits for frequency (min, max) in Hz.
        wavenumber_limits : tuple of float, optional
            x-axis limits for wavenumber (min, max) in cycles/m or rad/m.
        wavenumber_units : {'cycles', 'radians'}, default='cycles'
            Units of the wavenumber axis. Affects velocity overlays.
        use_filtered : bool, default=True
            Whether to use the filtered time series (`self.timeseries_filtered`) or raw (`self.timeseries`).
        to_db : {'normalized', 'raw', True, False}, default='normalized'
            If 'normalized', converts to decibel scale with 0 dB peak normalization.
            If 'raw', converts to decibels without normalization.
            If False, keeps the amplitude linear.
        mode_lines : dict, optional
            A dictionary of dispersion curves to overlay. Each key is a label (e.g., "Rayleigh"),
            and each value is either a velocity or a (velocity, color) tuple.

        Sets
        ----
        self.fk_freqs : ndarray
            Frequency axis in Hz.
        self.fk_wavenumbers : ndarray
            Wavenumber axis in cycles/m or rad/m.
        self.fk_spectrum : ndarray
            The 2D FFT power spectrum (in dB or linear units depending on `to_db`).
        self.fk_velocities : ndarray
            Phase velocity grid (freqs / wavenumbers).
        self.fk_power_linear : ndarray
            Linear (non-dB) power spectrum.
        self.fig_fk : matplotlib.figure.Figure
            The matplotlib figure object for the plot.
        self.ax_fk : matplotlib.axes.Axes
            The matplotlib axes object for the plot.

        Raises
        ------
        ValueError
            If FFT sizes are smaller than the input array dimensions, or if cutoff frequency limits are invalid.
        """
        if use_filtered: 
            dat = self.timeseries_filtered.copy() 
        else:
            dat = self.timeseries.copy() 
        
        m, n = dat.shape
        ntfft = ntfft or m
        nxfft = nxfft or n
        
        if ntfft < m or nxfft < n:
            raise ValueError("ntfft and nxfft must be >= data dimensions for zero-padding.")
        
        if taper:
            if 'tukey' in taper:
                from scipy.signal import tukey
                time_win = tukey(m, alpha=tukey_alpha)
                space_win = tukey(n, alpha=tukey_alpha)
            else:
                time_win = np.hanning(m)
                space_win = np.hanning(n)
            
            if taper in ['both', 'tukey_both']:
                dat *= time_win[:, None] * space_win[None, :]
            elif taper in ['time', 'tukey_time']:
                dat *= time_win[:, None]
            elif taper in ['space', 'tukey_space']:
                dat *= space_win[None, :]
        
        fk_spectrum = fftshift(fft2(dat, s = (ntfft, nxfft)) )
        
        epsilon = 1e-12
        if to_db == 'normalized':
            fk_spectrum = 10 * np.log10((np.abs(fk_spectrum)**2 + epsilon) / \
                np.max(np.abs(fk_spectrum)**2))
        elif to_db == 'raw':
            fk_spectrum = 10 * np.log10(np.abs(fk_spectrum)**2 + epsilon)
        
        # Create the frequency and wavenumber vectors then we need to shift the 
        # frequencies and wavenumbers for plotting 
        freqs = fftshift(np.fft.fftfreq(ntfft, self.dt))
        wavenumbers = -fftshift(np.fft.fftfreq(nxfft, d_rcx)) # positive down coordinates means we have to flip the f-k 
        velocities = freqs[:, None] / wavenumbers[None, :]  # cycles/m
        
        # Save the values 
        self.fk_freqs = freqs 
        self.fk_wavenumbers = wavenumbers 
        self.fk_spectrum = fk_spectrum 
        self.fk_velocities = velocities # multiply by 2*np.pi to get rad/m
        self.fk_power_linear = np.abs(fk_spectrum)**2

        mask = fk_spectrum < mask_db
        fk_spectrum_masked = np.ma.masked_where(mask, fk_spectrum)
        
        self.fig_fk, self.ax_fk = plt.subplots(figsize = figure_size)
        contour = self.ax_fk.contourf(
            wavenumbers, freqs, fk_spectrum_masked, 
            cmap = colormap, levels = contour_levels,
        )
        cbar = self.fig_fk.colorbar(contour, ax = self.ax_fk, label = 'Amplitude')
        self.ax_fk.set_xlabel(r'Wavenumber (m$^{-1}$)')
        self.ax_fk.set_ylabel('Frequency (Hz)')
        if frequency_limits:
            self.ax_fk.set_ylim(frequency_limits)
        
        if wavenumber_limits:
            self.ax_fk.set_xlim(wavenumber_limits)
        
        # overlay mode lines
        if mode_lines:
            for label, spec in mode_lines.items():
                if isinstance(spec, (tuple, list)):
                    v, col = spec
                else:
                    v, col = spec, 'white'
                # compute k_line in correct units
                if wavenumber_units == 'radians':
                    kline = 2 * np.pi * freqs / v
                else:
                    kline = freqs / v
                self.ax_fk.plot(kline, freqs, '--', color=col, lw=1, label=label)
                self.ax_fk.plot(-kline, freqs, '--', color=col, lw=1)
        
        plt.show()

    # -------------------------------------------------------------------------
    def multichannel_analysis(
            self,
            d_rcx: float,
            figure_size: Tuple[float, float] = (10, 8),
            colormap: str = 'viridis',
            contour_levels: int = 100,
            frequency_limits: Optional[Tuple[float, float]] = None,
            velocity_limits: Optional[Tuple[float, float]] = None,
            time_bandwidth: float = 4.0,
            n_tapers: int = 5,
            ntfft: Optional[int] = None,
            nxfft: Optional[int] = None,
            vmin = -40,
            vmax = None
        ):
        """
        Compute and plot the dispersion image using multitaper spectral analysis
        with DPSS tapers and f-k to f-c transformation.

        Parameters
        ----------
        d_rcx : float
            Receiver spacing (in meters).
        figure_size : tuple of float, optional
            Size of the output figure in inches.
        colormap : str, optional
            Matplotlib colormap name.
        contour_levels : int, optional
            Number of contour levels in the image.
        frequency_limits : tuple of float, optional
            Frequency axis limits (Hz).
        velocity_limits : tuple of float, optional
            Velocity axis limits (m/s).
        time_bandwidth : float, optional
            Time–bandwidth product for DPSS tapers.
        n_tapers : int, optional
            Number of DPSS tapers.
        ntfft : int, optional
            FFT size in time (zero-padded if larger than trace length).
        nxfft : int, optional
            FFT size in space (zero-padded if larger than number of traces).
        """
        m, n = self.timeseries.shape
        fs = 1.0 / self.dt

        ntfft = ntfft or m
        nxfft = nxfft or n

        # Generate DPSS tapers once
        tapers = windows.dpss(m, NW=time_bandwidth, Kmax=n_tapers)

        # Multitaper spectral estimation
        multitaper_spectrum = np.zeros((m // 2 + 1, n), dtype=np.complex128)
        for i in range(n):
            spectra = []
            for taper in tapers:
                tapered_data = self.timeseries[:, i] * taper
                spectrum = rfft(tapered_data)
                spectra.append(spectrum)
            multitaper_spectrum[:, i] = np.mean(spectra, axis=0)

        # f-k spectrum
        fk_spectrum = fftshift(fft2(multitaper_spectrum, s=(ntfft // 2 + 1, nxfft)))
        freqs = rfftfreq(ntfft, self.dt)
        wavenumbers = -fftshift(np.fft.fftfreq(nxfft, d_rcx))

        self.masw_spectrum = fk_spectrum
        self.masw_freqs = freqs

        # f-k to f-c
        fk_mesh_freq, fk_mesh_k = np.meshgrid(freqs, wavenumbers, indexing='ij')
        with np.errstate(divide='ignore', invalid='ignore'):
            phase_velocity = np.abs(fk_mesh_freq / fk_mesh_k)

        if velocity_limits is None:
            vmax = np.nanmax(phase_velocity[np.isfinite(phase_velocity)])
            velocity_limits = (0.0, vmax)
        velocity_bins = np.linspace(*velocity_limits, 500)
        self.masw_velocity_bins = velocity_bins

        # Bin f-k amplitude into f-c space
        amplitude_fc = np.zeros((len(freqs), len(velocity_bins)))
        for i, f_row in enumerate(np.abs(fk_spectrum)):
            c_indices = np.digitize(phase_velocity[i], velocity_bins)
            c_indices = np.clip(c_indices, 0, len(velocity_bins) - 1)
            np.add.at(amplitude_fc[i], c_indices, f_row)

        amplitude_fc /= n_tapers
        amplitude_fc = gaussian_filter(amplitude_fc, sigma=0.75)

        # Convert to dB and limit contrast
        epsilon = 1e-12
        amplitude_db = 10 * np.log10(amplitude_fc + epsilon)
        # amplitude_db = np.clip(amplitude_db, -40, 0)

        # Apply frequency limits to the plotted range
        if frequency_limits:
            fmin, fmax = frequency_limits
            freq_mask = (freqs >= fmin) & (freqs <= fmax)
            freqs = freqs[freq_mask]
            amplitude_db = amplitude_db[freq_mask, :]
        
        mask = amplitude_db < -20
        amplitude_db_masked = np.ma.masked_where(mask, amplitude_db)
        # Plot frequency on x and velocity on y
        self.fig_dispersion, self.ax_dispersion = plt.subplots(figsize=figure_size)
        # Clip everything below -40 dB for visualization purposes
        if not vmax:
            vmax = np.max(amplitude_db)

        # Create colormap that clips under vmin
        cmap = plt.get_cmap(colormap).copy()
        cmap.set_under(cmap(0))  # Optional: choose a different color if desired

        # Plot using imshow (note: extent uses [x0, x1, y0, y1])
        im = self.ax_dispersion.imshow(
            amplitude_db_masked.T,
            extent=[freqs[0], freqs[-1], velocity_bins[0], velocity_bins[-1]],
            origin='lower',
            aspect='auto',
            cmap=cmap,
            vmin=vmin,
            vmax=vmax
        )

        # Add colorbar
        cbar = self.fig_dispersion.colorbar(im, ax=self.ax_dispersion, label='Amplitude (dB)')

        # Axis labels
        self.ax_dispersion.set_xlabel('Frequency (Hz)')
        self.ax_dispersion.set_ylabel('Phase Velocity (m/s)')

        # Optional axis limits
        if frequency_limits:
            self.ax_dispersion.set_xlim(frequency_limits)
        if velocity_limits:
            self.ax_dispersion.set_ylim(velocity_limits)
        
        plt.tight_layout()
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
    
    
    
