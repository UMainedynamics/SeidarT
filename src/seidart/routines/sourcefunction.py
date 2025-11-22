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
    'find_zero_crossing',
    'double_couple_tensor'
]

# ================================ Definitions ================================

def wavelet(timevec: np.ndarray, f: float, stype: str, omni=False) -> np.ndarray:
    """
    Generates a wavelet based on the specified type and parameters.
    """
    a = (2.0 * np.pi * f) ** 2
    to = 1.0 / f
    if stype == 'gaus0':
        bw = 0.2 / (np.pi * f) 
        a = (timevec - to) / bw
        x = np.exp(-a*a)
    elif stype == 'gaus1':
        a = np.pi * f * (timevec - to) 
        e = np.exp(-a * a) 
        x = (1.0 - 2.0*a*a) * e 
    elif stype == 'gaus2':
        dt = (timevec - to) 
        k = np.pi * np.pi * f * f 
        e = np.exp(-k*dt*dt)
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
def double_couple_tensor(strike_deg: float, dip_deg: float, rake_deg: float,
                         M0: float = 1.0) -> np.ndarray:
    """
    Return (Mxx, Myy, Mzz, Mxy, Mxz, Myz) for a unit double-couple (scaled by M0).
    """
    strike = np.deg2rad(strike_deg)
    dip    = np.deg2rad(dip_deg)
    rake   = np.deg2rad(rake_deg)
    
    # Fault normal n and slip direction d (global coordinates)
    # (Aki & Richards 2002; n points out of the fault plane, d is in-plane slip)
    s_hat = np.array([-np.sin(strike), np.cos(strike), 0.0])
    d_hat  = np.array([ np.cos(strike)*np.cos(dip), np.sin(strike)*np.cos(dip), -np.sin(dip)])  # down-dip
    s_hat /= np.linalg.norm(s_hat)
    d_hat /= np.linalg.norm(d_hat)
    n_hat  = np.cross(s_hat, d_hat)
    n_hat /= np.linalg.norm(n_hat)
    
    # slip direction from rake (in-plane rotation)
    slip = np.cos(rake)*s_hat + np.sin(rake)*d_hat
    slip /= np.linalg.norm(slip)
    
    # DC tensor
    M = M0*(np.outer(slip, n_hat) + np.outer(n_hat, slip))
    M = 0.5 * (M + M.T) 
    M -= np.trace(M)/3.0 * np.eye(3) 
    
    return M

def isotropic_tensor(M0: float = 1.0) -> np.ndarray:
    """
    Isotropic moment tensor 
    M0 > 0 - explosion; M0 < 0 - implosion
    """
    return M0 * np.eye(3) 

# ----------------------------------------------------------------------------
def clvd_tensor(axis_vec: np.ndarray, M0: float = 1.0, mode: str = "tensile") -> np.ndarray:
    """
    Build a CLVD moment tensor with principal axis along `axis_vec`.
    mode: "tensile" -> eigvals (2, -1, -1)  (axis is T)
          "compressional" -> eigvals (-2, 1, 1) (axis is P)
    M0 scales the tensor (N·m or arbitrary).
    """
    a = np.asarray(axis_vec, dtype=float)
    a /= (np.linalg.norm(a) + 1e-15)  # e3
    # pick e1 not parallel to a
    tmp = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = tmp - np.dot(tmp, a) * a
    e1 /= (np.linalg.norm(e1) + 1e-15)
    e2 = np.cross(a, e1)
    R = np.column_stack((e1, e2, a))  # columns are eigenvectors
    
    if mode == "tensile":
        lam = np.diag([1.0, -0.5, -0.5]) * 2.0  # -> (2, -1, -1)
    elif mode == "compressional":
        lam = np.diag([-1.0, 0.5, 0.5]) * 2.0   # -> (-2, 1, 1)
    else:
        raise ValueError("mode must be 'tensile' or 'compressional'")
    
    M = (R @ lam @ R.T) * M0
    # hygiene: enforce symmetry & zero trace numerically
    M = 0.5*(M + M.T)
    M -= np.trace(M)/3.0 * np.eye(3)
    return M

# ----------------------------------------------------------------------------
def pointsource(
        modelclass,
        multimodal: bool = False,
        broadband: bool = False,
        fmin: float = None,
        fmax: float = None,
        center: bool = False,
        num_octaves: int = 2,
        source: str = "ac",
        mt: np.ndarray = None,
        a_iso=1.0,                     # weight for ISO (diagonal) part → omni P
        b_dc=0.7,  
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build and write point-source time series.
    
    source can be:
        "ac" - accelerated weight drop
        "dc" - double couple 
        "tnt" - explosive 
        "omni_ps" - omnidirectional p and s waves; 3 direction dc plus explosive source (for testing)
    """
    N = int(modelclass.time_steps)
    dt = float(modelclass.dt)
    timevec = np.arange(N) * dt
    f0 = float(modelclass.f0)
    amp = float(modelclass.source_amplitude)

    if broadband:
        srcfn = amp * broadbandsrc(timevec, fmin, fmax)
    elif multimodal:
        srcfn = amp * multimodesrc(
            timevec, f0, modelclass.source_wavelet, center=center, num_octaves = num_octaves
        )
    else:
        if center:
            srcfn = amp * wavelet_center0(timevec, f0, modelclass.source_wavelet)
        else:
            srcfn = amp * wavelet(timevec, f0, modelclass.source_wavelet)
    
    # Pre-allocate output time series 
    forcex = np.zeros(N)
    forcey = np.zeros(N) 
    forcez = np.zeros(N) 
    sigma_xx = np.zeros(N)
    sigma_yy = np.zeros(N)
    sigma_zz = np.zeros(N)
    sigma_xy = np.zeros(N)
    sigma_xz = np.zeros(N)
    sigma_yz = np.zeros(N)
    
    if modelclass.is_seismic:
        if modelclass.source_type=='ac':  
            R = rotate_sdr(modelclass.phi, modelclass.theta, modelclass.psi)       
            forcez = R[0,2] * srcfn
            forcey = R[1,2] * srcfn
            forcex = R[2,2] * srcfn      
        elif modelclass.source_type=='dc':
            M = double_couple_tensor(modelclass.phi, modelclass.theta, modelclass.psi, M0=1.0)
            modelclass.moment_tensor = M.copy() 
            sigma_xx = srcfn * M[0,0]
            sigma_yy = srcfn * M[1,1] 
            sigma_zz = srcfn * M[2,2]
            sigma_xy = srcfn * M[0,1] 
            sigma_xz = srcfn * M[0,2] 
            sigma_yz = srcfn * M[1,2]
        elif modelclass.source_type == 'tnt':
            # Isotropic stress-rate: σxx=σyy=σzz = srcfn, shears = 0 
            sigma_xx = srcfn.copy() 
            sigma_yy = srcfn.copy() 
            sigma_zz = srcfn.copy() 
            modelclass.moment_tensor = isotropic_tensor(M0=1.0)
    
        elif modelclass.source_type == 'omni_ps':    
            R = rotate_sdr(modelclass.phi, modelclass.theta, modelclass.psi)    
            # First omni-P/explosion
            s_iso = a_iso * srcfn 
            sigma_xx += s_iso
            sigma_yy += s_iso 
            sigma_zz += s_iso 
            # Now omin-S/3 double couples 
            s_dc = b_dc * srcfn 
            e1, e2, e3 = R[:,0], R[:,1], R[:,2]
            DCs = [
                np.outer(e1, e2) + np.outer(e2, e1),
                np.outer(e1, e3) + np.outer(e3, e1),
                np.outer(e2, e3) + np.outer(e3, e2),
            ]
            Mavg = sum(DCs) / 3.0
            sigma_xx += s_dc * Mavg[0,0]
            sigma_yy += s_dc * Mavg[1,1]
            sigma_zz += s_dc * Mavg[2,2]
            sigma_xy += s_dc * Mavg[0,1]
            sigma_xz += s_dc * Mavg[0,2]
            sigma_yz += s_dc * Mavg[1,2]
        elif modelclass.source_type == 'omni_s':    
            R = rotate_sdr(modelclass.phi, modelclass.theta, modelclass.psi)    
            # omin-S/3 double couples 
            s_dc = b_dc * srcfn 
            e1, e2, e3 = R[:,0], R[:,1], R[:,2]
            DCs = [
                np.outer(e1, e2) + np.outer(e2, e1),
                np.outer(e1, e3) + np.outer(e3, e1),
                np.outer(e2, e3) + np.outer(e3, e2),
            ]
            Mavg = sum(DCs) / 3.0
            sigma_xx += s_dc * Mavg[0,0]
            sigma_yy += s_dc * Mavg[1,1]
            sigma_zz += s_dc * Mavg[2,2]
            sigma_xy += s_dc * Mavg[0,1]
            sigma_xz += s_dc * Mavg[0,2]
            sigma_yz += s_dc * Mavg[1,2]
        else: 
            raise ValueError(f"Unkown source type '{source}'")
        
        writesrc("seismicsourcex.dat", forcex)
        writesrc("seismicsourcey.dat", forcey)
        writesrc("seismicsourcez.dat", forcez)
        
        writesrc("seismicsourcexx.dat", sigma_xx)
        writesrc("seismicsourcexy.dat", sigma_xy)
        writesrc("seismicsourcexz.dat", sigma_xz)
        writesrc("seismicsourceyy.dat", sigma_yy)
        writesrc("seismicsourceyz.dat", sigma_yz)
        writesrc("seismicsourcezz.dat", sigma_zz)
        
    else:
        R = rotator_sdr(modelclass.phi. modelclass.theta, modelclass.psi)       
        forcez = R[0,2] * srcfn
        forcey = R[1,2] * srcfn
        forcex = R[2,2] * srcfn  
        writesrc("electromagneticsourcex.dat", forcex)
        writesrc("electromagneticsourcey.dat", forcey)
        writesrc("electromagneticsourcez.dat", forcez)
    
    modelclass.source_sigma_xx = sigma_xx
    modelclass.source_sigma_xy = sigma_xy 
    modelclass.source_sigma_xz = sigma_xz 
    modelclass.source_sigma_yy = sigma_yy 
    modelclass.source_sigma_yz = sigma_yz 
    modelclass.source_sigma_zz = sigma_zz 
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
#             timevec, f0, modelclass.source_wavelet, center=center, num_octaves = num_octaves
#         )
#     else:
#         if center:
#             srcf n = modelclass.source_amplitude * wavelet_center0(timevec, f0, modelclass.source_wavelet)
#         else:
#             srcfn = modelclass.source_amplitude * wavelet(timevec, f0, modelclass.source_wavelet)

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
def main(prjfile: str, source_wavelet: str, factor: float, multimodal: bool, make_plot: bool) -> None:
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
    source_wavelet = args.sourcetype[0]
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
