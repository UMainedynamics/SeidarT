
import argparse
import numpy as np 
from seidart.routines.definitions import *
from scipy import signal
import matplotlib.pyplot as plt 
from scipy.io import FortranFile 
import scipy.signal
from typing import Tuple 

# ================================ Definitions ================================
def wavelet(timevec: np.ndarray, f: float, stype: str) -> np.ndarray:
    """
    Generates a wavelet based on the specified type and parameters.
    
    :param timevec: The time vector for the wavelet.
    :param f: The center frequency of the wavelet.
    :param stype: The type of wavelet to generate. Options include 'gaus0', 
                    'gaus1', 'gaus2', 'chirp', and 'chirplet'.
    :type timevec: np.ndarray
    :type f: float
    :type stype: str
    :return: The generated wavelet.
    :rtype: np.ndarray
    """
    
    # Create the wavelet given the parameters
    a = np.pi**2 * f**2
    to = 1/f
    if stype == 'gaus0':
        x = np.exp(-a*(timevec - to)**2)
    if stype == "gaus1":
        x = - 2.0 * a * (timevec - to) * np.exp(-a * ((timevec - to)**2))    
    if stype == "gaus2":
        x = 2.0 * a * np.exp(-a * (timevec - to)**2) * (2.0 * a * (timevec - to)**2 - 1)
    if stype == "chirp":
        f
        x = signal.chirp(timevec, 10*f, to, f, phi = -90)    
    if stype == "chirplet":
        x = signal.chirp(timevec, f, to, 20*f, phi = -90)
        g = np.exp(-(a/4)*(timevec - to)**2)
        x = x * g        
    x = x/x.max()
    return(x)

def multimodesrc(timevec: np.ndarray, f: float, stype: str) -> np.ndarray:
    """
    Creates a multi-mode source by linearly combining wavelets at different 
    frequencies.
    
    :param timevec: The time vector for the source.
    :param f: The center frequency around which the modes are generated.
    :param stype: The type of wavelet to use for each mode.
    :type timevec: np.ndarray
    :type f: float
    :type stype: str
    :return: The combined multi-mode source.
    :rtype: np.ndarray
    """
    # Create a double octave sweep centered at f0 from the addition of multiple
    # sources. 
    # The change will be linear in 1/8 octave steps
    fmin = f0/4 
    fmax = f0
    df = (f0 - f0/2)/8
    f = np.arange(fmin, fmax, df)
    stf = np.zeros([len(timevec)])
    for freq in f:
        stf = stf + wavelet(timevec, freq, stype)
    return(stf)

def plotsource(t: np.ndarray, x: np.ndarray) -> Tuple[plt.Figure, np.ndarray]:
    """
    Plots the source function and its power spectrum.
    
    :param t: The time vector corresponding to the source function.
    :param x: The source function.
    :type t: np.ndarray
    :type x: np.ndarray
    :return: A tuple containing the matplotlib figure and axes array.
    :rtype: Tuple[plt.Figure, np.ndarray]
    """
    fs = 1/np.mean(np.diff(t) )
    f, pxx = signal.welch(x, fs = fs)
    db = 10*np.log10(pxx)
    fig, ax = plt.subplots(2)
    ax[0].plot(t, x, '-b')
    ax[0].set( xlabel = 'Time (s)', ylabel= 'Amplitude')
    ax[1].plot(f, db, '-b')
    ax[1].set(xlabel = 'Frequency (Hz)', ylabel = 'Power (dB)')
    ax[1].set_xlim([f.min(), np.min([20*f0, f.max()])] )
    return(fig,ax)

def writesrc(fn: str, srcarray: np.ndarray) -> None:
    """
    Writes a source array to a Fortran-formatted binary file.
    
    :param fn: The filename for the output file.
    :param srcarray: The source array to write.
    :type fn: str
    :type srcarray: np.ndarray
    """
    f = FortranFile(fn, 'w')
    f.write_record(srcarray)
    f.close()

def sourcefunction(
        modelclass: Model, 
        factor: float, 
        source_type: str, 
        model_type: str, 
        multimodal: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Generates a source function based on the provided model class and parameters.
    
    :param modelclass: The model class containing time steps, frequency, and 
        orientation.
    :param factor: Amplitude scaling factor for the source.
    :param source_type: The type of source to generate.
    :param model_type: Specifies whether the source is for a seismic ('s') or 
        electromagnetic ('e') model.
    :param multimodal: Whether to generate a multimodal source. Defaults to False.
    :type modelclass: Model
    :type factor: float
    :type source_type: str
    :type model_type: str
    :type multimodal: bool
    :return: A tuple containing the time vector and the x, y, z components of 
        the force, along with the source function.
    :rtype: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    """
    # Create the source 
    N = int(modelclass.time_steps)
    timevec = np.linspace(1, N, num = N ) * \
        float(modelclass.dt)
    f0 = float(modelclass.f0)
    
    # Create the source function
    if multimodal:
        srcfn = factor * multimodesrc(timevec, f0, source_type)
    else:
        srcfn = factor * wavelet(timevec, f0, source_type)
    # rotate 
    theta = np.pi * modelclass.theta * 180
    phi = np.pi * modelclass.phi * 180
    forcex = np.sin( theta ) * np.cos( phi ) * srcfn
    forcey = np.sin( theta ) * np.sin( phi ) * srcfn
    forcez = np.cos( theta ) * srcfn
    if model_type == 's' or model_type == 'seismic':
        writesrc("seismicsourcex.dat", forcex)
        writesrc("seismicsourcey.dat", forcey)
        writesrc("seismicsourcez.dat", forcez)
    if model_type == 'e' or model_type == 'electromag':
        writesrc("electromagneticsourcex.dat", forcex)
        writesrc("electromagneticsourcey.dat", forcey)
        writesrc("electromagneticsourcez.dat", forcez)
    return(timevec, forcex, forcey, forcez, srcfn)

def main(
        prjfile: str, 
        source_type: str, 
        factor: float, 
        model_type: str,
        multimodal: bool, 
        make_plot: bool
    ) -> None:
    """
    Main function to generate and optionally plot a source function based on 
    input parameters.
    
    :param prjfile: Path to the project file.
    :param source_type: Type of the source function to generate.
    :param factor: Amplitude scaling factor for the source.
    :param model_type: Specifies the model type ('s' for seismic, 'e' for 
        electromagnetic).
    :param multimodal: Indicates if a multimodal source should be generated. 
        Defaults to False.
    :param make_plot: Whether to plot the generated source function. Defaults 
        to False.
    :type prjfile: str
    :type source_type: str
    :type factor: float
    :type model_type: str
    :type multimodal: bool
    :type make_plot: bool
    """
    # Load the project file 
    # Let's initiate the domain
    domain, material, seismic, electromag = loadproject(
        prjfile, 
        Domain(), 
        Material(), 
        Model(), 
        Model()
    )
    if model_type == 's' or model_type == 'seismic':
        timevec, fx, fy, fz, srcfn = sourcefunction(
            seismic, 
            factor, 
            source_type, 
            model_type, 
            multimodal = multimodal
        )
    if model_type == 'e' or model_type == 'electromagnetic':
        timevec, fx, fy, fz, srcfn = sourcefunction(
            electromag, 
            factor, 
            source_type, 
            model_type, 
            multimodal = multimodal
        )
    if make_plot:
        plotsource(timevec, srcfn)
        plt.show()

# --------------------------- Command Line Arguments --------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = """We support building a few different source time 
        functions and writing them to a text file. From a specified project 
        file we can create the time series for a source function and save it 
        in a Fortran formatted .dat file. The default is an impulse function 
        with center frequency f0, but we can also generate other wavelets and 
        chirplets.
        """
    )
    
    parser.add_argument(
        '-p', '--projectfile', nargs = 1, type = str, required = True,
        help = """The path to the project file"""
    )
    
    parser.add_argument(
        '-S', '--sourcetype', nargs = 1, type = str, required = False, 
        default = "gaus0",
        help = """Specify the source type. Available wavelets are: 
        gaus0, gaus1, gaus2 (gaussian n-th derivative), chirp, chirplet, 
        multimodal. (Default = gaus0)"""
    )
    
    parser.add_argument(
        '-m', '--modeltype', nargs = 1, type = str, required = False, 
        default = 's',
        help = """Specify whether to construct the source for an em or seismic
        model. s-seismic, e-electromagnetic, b-both"""
    )
    
    parser.add_argument(
        '-a', '--amplitude', nargs = 1, type = float, required = False,
        default = 1.0,
        help = """Input the scalar factor for source amplification. 
        (Default = 1.0)"""
    )
    
    parser.add_argument(
        '-M', '--multimodal', action='store_true',
        help = """Multimodal source is computed across 2 octaves at 1/8 steps 
        and centered at f0 """
    )
    
    parser.add_argument(
        '-P', '--plot', action='store_true',
        help = """Plot the source and spectrum"""
    )
    
    args = parser.parse_args()
    prjfile = ''.join(args.projectfile)
    source_type = ''.join(args.sourcetype)
    factor = args.amplitude[0]
    model_type = ''.join(args.modeltype)
    multimodal = args.multimodal 
    make_plot = args.plot
    #
    
    main(prjfile, source_type, factor, model_type, multimodal, make_plot)
    







