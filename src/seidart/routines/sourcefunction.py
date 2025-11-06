import argparse
import numpy as np
from seidart.routines.definitions import *
from seidart.routines.materials import * 
from scipy import signal
import numpy.fft as fft
from scipy.signal import butter, filtfilt, firwin, minimum_phase
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from typing import Tuple

__all__ = [
    'wavelet',
    'multimodesrc',
    'broadbandsrc',
    'plotsource',
    'writesrc',
    'pointsource',
    'find_zero_crossing'
]

# ================================ Definitions ================================

def wavelet(timevec: np.ndarray, f: float, stype: str, omni=False) -> np.ndarray:
    """
    Generates a wavelet based on the specified type and parameters.
    """
    a = (2.0 * np.pi * f) ** 2
    to = 1.0 / f
    if stype == 'gaus0':
        bw = 0.2 / (np.pi * fc) 
        a = (timevec - to) / bw
        x = np.exp(-a*a)
    elif stype == 'gaus1':
        a = pi * fc * (timevec - t0) 
        e = exp(-a * a) 
        x = (1.0 - 2.0*a*a) * e 
    elif stype == 'gaus2':
        dt = (timevec - to) 
        k = np.pi * np.pi * fc * fc 
        e = exp(-k*dt*dt)
        x = 2.0 * k*dt * (2.0*k*dt*dt - 3.0) * e
    elif stype == 'chirp':
        x = signal.chirp(timevec, 10 * f, to, f, phi=-90)
    elif stype == 'chirplet':
        x = signal.chirp(timevec, f, to, 20 * f, phi=-90)
        g = np.exp(-(a / 4) * (timevec - to) ** 2)
        x = x * g
    else:
        # fallback: normalize input series if known source
        try:
            x = np.array(stype)
        except:
            raise ValueError(f"Unknown wavelet type: {stype}")
    
    return x / np.max(np.abs(x))

# ----------------------------------------------------------------------------
def wavelet_center0(timevec: np.ndarray, f: float, stype: str) -> np.ndarray:
    """
    Same as wavelet but centered at t=0 for simultaneous onset.
    Only supports gaus0, gaus1, gaus2; others fall back to center-shifted wavelet.
    """
    a = (2.0 * np.pi * f) ** 2
    t = timevec
    if stype == 'gaus0':
        x = np.exp(-a * t**2)
    elif stype == 'gaus1':
        x = -2.0 * a * t * np.exp(-a * t**2)
    elif stype == 'gaus2':
        x = 2.0 * a * np.exp(-a * t**2) * (2.0 * a * t**2 - 1)
    else:
        # fallback: shift original wavelet
        x = wavelet(timevec + 1.0/f, f, stype)
    return x / np.max(np.abs(x))

# ----------------------------------------------------------------------------

def multimodesrc(
        timevec: np.ndarray,
        f: float,
        stype: str,
        center: bool = False, 
        num_octaves = 3
    ) -> np.ndarray:
    """
    Creates a multi-mode source by combining octave-spaced wavelets around f.
    If center=True, each wavelet is shifted to t=0 for simultaneous onset.
    """
    freqs = f * 2 ** np.linspace(-num_octaves, 0, num_octaves * 4 - 1)
    stf = np.zeros_like(timevec)
    for freq in freqs:
        if center:
            stf += wavelet_center0(timevec, freq, stype)
        else:
            stf += wavelet(timevec, freq, stype)
    return stf / np.max(np.abs(stf))

# ----------------------------------------------------------------------------

def bandlimited_impulse_fir(
        timevec: np.ndarray,
        fmin: float,
        fmax: float,
        numtaps: int = 201,
        window: str = 'hann'
    ) -> np.ndarray:
    """
    Minimum-phase band-limited impulse via FIR design:
    - flat response in [fmin, fmax]
    - minimal pre-ringing (causal)

    :param timevec: sample times
    :param fmin: low cutoff (Hz)
    :param fmax: high cutoff (Hz)
    :param numtaps: FIR filter length (odd recommended)
    :param window: window for firwin
    :return: impulse response of length len(timevec)
    """
    N = len(timevec)
    fs = 1.0 / np.mean(np.diff(timevec))
    # design linear-phase bandpass FIR
    coeffs = firwin(numtaps, [fmin, fmax], pass_zero=False, fs=fs, window=window)
    # convert to minimum-phase to remove pre-ringing
    mp = minimum_phase(coeffs)
    # causal impulse: pad mp at start
    h = np.zeros(N)
    L = len(mp)
    h[:L] = mp
    return h / np.max(np.abs(h))

# ----------------------------------------------------------------------------

def broadbandsrc(timevec: np.ndarray, fmin: float, fmax: float, **kwargs) -> np.ndarray:
    """
    Broadband source: all frequencies in [fmin, fmax] onset simultaneously
    using a minimum-phase FIR impulse response.
    """
    return bandlimited_impulse_fir(timevec, fmin, fmax,
                                   numtaps=kwargs.get('numtaps', 201),
                                   window=kwargs.get('window', 'hann'))

# ----------------------------------------------------------------------------

def plotsource(t: np.ndarray, x: np.ndarray) -> Tuple[plt.Figure, np.ndarray]:
    """
    Plots the source function and its power spectrum.
    """
    fs = 1 / np.mean(np.diff(t))
    f, pxx = signal.welch(x, fs=fs)
    db = 10 * np.log10(pxx)
    fig, ax = plt.subplots(2)
    ax[0].plot(t, x, '-b')
    ax[0].set(xlabel='Time (s)', ylabel='Amplitude')
    ax[1].plot(f, db, '-b')
    ax[1].set(xlabel='Frequency (Hz)', ylabel='Power (dB)')
    ax[1].set_xlim([f.min(), np.min([20 * f0, f.max()])])
    return fig, ax

# ----------------------------------------------------------------------------

def writesrc(fn: str, srcarray: np.ndarray) -> None:
    """
    Writes a source array to a Fortran-formatted binary file.
    """
    f = FortranFile(fn, 'w')
    f.write_record(srcarray)
    f.close()

# ----------------------------------------------------------------------------
def pointsource(
        modelclass,
        multimodal: bool = False,
        broadband: bool = False,
        fmin: float = None,
        fmax: float = None,
        center: bool = False,
        num_octaves: int = 2,
        stress_output: bool = True,
        a_iso=1.0,                     # weight for ISO (diagonal) part → omni P
        b_dc=0.7,  
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build and write point-source time series.
    """
    N = int(modelclass.time_steps)
    dt = float(modelclass.dt)
    timevec = np.arange(N) * dt
    f0 = float(modelclass.f0)

    if broadband:
        srcfn = modelclass.source_amplitude * broadbandsrc(timevec, fmin, fmax)
    elif multimodal:
        srcfn = modelclass.source_amplitude * multimodesrc(
            timevec, f0, modelclass.source_type, center=center, num_octaves = num_octaves
        )
    else:
        if center:
            srcf n = modelclass.source_amplitude * wavelet_center0(timevec, f0, modelclass.source_type)
        else:
            srcfn = modelclass.source_amplitude * wavelet(timevec, f0, modelclass.source_type)

    theta = np.deg2rad(modelclass.theta)
    phi = np.deg2rad(modelclass.phi)
    psi = np.deg2rad(modelclass.psi) 
    
    R = rotator_zxz(np.array([theta, phi, psi) )

    forcez = np.sin(theta) * np.cos(phi) * srcfn
    forcey = np.sin(theta) * np.sin(phi) * srcfn
    forcex = np.cos(theta) * srcfn

    if modelclass.is_seismic:
        if stress_output:
            
            s_iso = a_iso * srcfn 
            s_dc = b_dc * srcfn 
            
            sigma_xx = s_iso.copy() 
            sigma_yy = s_iso.copy() 
            sigma_zz = s_iso.copy() 
            sigma_xy = s_dc.copy() 
            sigma_xz = s_dc.copy() 
            simga_yz = s_dc.copy() 
            
            
            writesrc("seismicsourcexx.dat", sigma_xx)
            writesrc("seismicsourcexy.dat", sigma_xy)
            writesrc("seismicsourcexz.dat", sigma_xz)
            writesrc("seismicsourceyy.dat", sigma_yy)
            writesrc("seismicsourceyz.dat", sigma_yz)
            writesrc("seismicsourcezz.dat", sigma_zz)
            
            modelclass.source_sigma_xx = sigma_xx
            modelclass.source_sigma_xy = sigma_xy 
            modelclass.source_sigma_xz = sigma_xz 
            modelclass.source_sigma_yy = sigma_yy 
            modelclass.source_sigma_yz = sigma_yz 
            modelclass.source_sigma_zz = sigma_zz 
        else:
            writesrc("seismicsourcex.dat", forcex)
            writesrc("seismicsourcey.dat", forcey)
            writesrc("seismicsourcez.dat", forcez)
    else:
        writesrc("electromagneticsourcex.dat", forcex)
        writesrc("electromagneticsourcey.dat", forcey)
        writesrc("electromagneticsourcez.dat", forcez)

    modelclass.sourcefunction_x = forcex
    modelclass.sourcefunction_y = forcey
    modelclass.sourcefunction_z = forcez
    return timevec, forcex, forcey, forcez, srcfn

# ----------------------------------------------------------------------------
# def omnisource(
#         modelclass,
#         multimodal: bool = False,
#         broadband: bool = False,
#         fmin: float = None,
#         fmax: float = None,
#         center: bool = False,
#         num_octaves: int = 2
#     ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
#     """
#     Build and write point-source time series.
#     """
#     N = int(modelclass.time_steps)
#     dt = float(modelclass.dt)
#     timevec = np.arange(N) * dt
#     f0 = float(modelclass.f0)

#     if broadband:
#         srcfn = modelclass.source_amplitude * broadbandsrc(timevec, fmin, fmax)
#     elif multimodal:
#         srcfn = modelclass.source_amplitude * multimodesrc(
#             timevec, f0, modelclass.source_type, center=center, num_octaves = num_octaves
#         )
#     else:
#         if center:
#             srcf n = modelclass.source_amplitude * wavelet_center0(timevec, f0, modelclass.source_type)
#         else:
#             srcfn = modelclass.source_amplitude * wavelet(timevec, f0, modelclass.source_type)

#     theta = np.deg2rad(modelclass.theta)
#     phi = np.deg2rad(modelclass.phi)
#     forcez = np.sin(theta) * np.cos(phi) * srcfn
#     forcey = np.sin(theta) * np.sin(phi) * srcfn
#     forcex = np.cos(theta) * srcfn

#     if modelclass.is_seismic:
#         if stress_output:
#             writesrc("seismicsourcexx.dat", forcexx)
#             writesrc("seismicsourcexy.dat", forcexy)
#             writesrc("seismicsourcexz.dat", forcexz)
#             writesrc("seismicsourceyy.dat", forceyy)
#             writesrc("seismicsourceyz.dat", forceyz)
#             writesrc("seismicsourcezz.dat", forcezz)    
#         else:
#             writesrc("seismicsourcex.dat", forcex)
#             writesrc("seismicsourcey.dat", forcey)
#             writesrc("seismicsourcez.dat", forcez)
#     else:
#         writesrc("electromagneticsourcex.dat", forcex)
#         writesrc("electromagneticsourcey.dat", forcey)
#         writesrc("electromagneticsourcez.dat", forcez)

#     modelclass.sourcefunction_x = forcex
#     modelclass.sourcefunction_y = forcey
#     modelclass.sourcefunction_z = forcez
#     return timevec, forcex, forcey, forcez, srcfn




# ----------------------------------------------------------------------------

def find_zero_crossing(timevec, wavelet_vals, noise_region=slice(0, 100), k=3):
    sigma = np.std(wavelet_vals[noise_region])
    tol = k * sigma
    filtered_vals = np.where(np.abs(wavelet_vals) < tol, 0, wavelet_vals)
    sign_changes = np.where(np.diff(np.sign(filtered_vals)) != 0)[0]
    if len(sign_changes) == 0:
        return None
    i0 = sign_changes[0]
    t0_approx = timevec[i0] + (timevec[i0+1] - timevec[i0]) * (
        -filtered_vals[i0] / (filtered_vals[i0+1] - filtered_vals[i0])
    )
    return t0_approx

# ----------------------------------------------------------------------------
def main(prjfile: str, source_type: str, factor: float, multimodal: bool, make_plot: bool) -> None:
    parser = argparse.ArgumentParser(
        description="""Generate source time functions in Fortran .dat format."""
    )
    parser.add_argument('-p', '--projectfile', nargs=1, type=str, required=True)
    parser.add_argument('-S', '--sourcetype', nargs=1, type=str, default="gaus0")
    parser.add_argument('-m', '--modeltype', nargs=1, type=str, default='s')
    parser.add_argument('-a', '--amplitude', nargs=1, type=float, default=[1.0])
    parser.add_argument('-M', '--multimodal', action='store_true')
    parser.add_argument('-B', '--broadband', action='store_true',
                        help="Enable broadband source between fmin and fmax")
    parser.add_argument('--fmin', type=float, help="Minimum frequency for broadband")
    parser.add_argument('--fmax', type=float, help="Maximum frequency for broadband")
    parser.add_argument('-P', '--plot', action='store_true')
    args = parser.parse_args()

    prjfile = args.projectfile[0]
    source_type = args.sourcetype[0]
    factor = args.amplitude[0]
    model_type = args.modeltype[0]
    multimodal = args.multimodal
    broadband_flag = args.broadband
    fmin = args.fmin
    fmax = args.fmax
    make_plot = args.plot

    domain, material, seismic, electromag = loadproject(
        prjfile, Domain(), Material(), Model(), Model()
    )
    if model_type.lower().startswith('s'):
        timevec, fx, fy, fz, srcfn = pointsource(
            seismic, multimodal, broadband_flag, fmin, fmax
        )
    else:
        timevec, fx, fy, fz, srcfn = pointsource(
            electromag, multimodal, broadband_flag, fmin, fmax
        )
    if make_plot:
        plotsource(timevec, srcfn)
        plt.show()

if __name__ == "__main__":
    main()
