#!/usr/bin/env python3 

"""
Compute the cpml matrix and write it to a Fortran unformatted binary
"""

import numpy as np 
from scipy.io import FortranFile
from definitions import *
import argparse

parser = argparse.ArgumentParser(
    description = "Compute the cpml for all boundaries"
)

parser.add_argument(
    '-p', '--prjfile', 
    nargs = 1, type = str, required = True,
    help = """Path to the project file."""
)

parser.add_argument(
    '-H', '--half', action='store_true'
)
parser.add_argument(
    '-d', '--direction',
    nargs = 1, type = str, required = True,
    help = """Specify the direction (i.e. x, y, or z) to compute the cpml 
    coefficients."""
)

args = parser.parse_args()
prjfile = ''.join(args.prjfile)
directory = ''.join(args.direction)
half = args.half

# -----------------------------------------------------------------------------

# prjfile = "twolayer.prj"
# directory = 'x'
# half = False

# Load the project file
domain, material, seismic, electromag = loadproject(
    prjfile,
    Domain(),
    Material(),
    Model(),
    Model()
)


domain.nx = domain.nx + 2*domain.cpml
domain.ny = domain.ny + 2*domain.cpml
domain.nz = domain.nz + 2*domain.cpml
dim = domain.dim
domain.dim = 2

# Load sigma_xx, sigma_yy, or sigma_zz since these are going to be the maximum
# values.
if directory == 'x':
    N = domain.nx
    dx = domain.dx
elif directory == 'y':
    N = domain.ny 
    dx = domain.dy 
else:
    N = domain.nz
    dx = domain.dz 

# -----------------------------------------------------------------------------
# Define constants
NP = 2 
NPA = 2 
k_max = 1.1e1 # This is the value determined from the litterature. 
eps0 = 8.85418782e-12
mu0 = 4.0*np.pi*1.0e-7

# Compute the distance along the absorbing boundary relative to the end of the 
# original model space. 
dist = dx * np.arange(0, domain.cpml)
if half:
    dist = dist + dx/2 

dist = dx*domain.cpml - dist
dist = dist/(dx*domain.cpml)

# Allocate space
alpha_max = 2*np.pi*eps0*electromag.f0/10
sig_max = 0.7 * (NP+1) / (dx * np.sqrt(mu0/eps0) )

kappa = np.ones([N])
alpha = np.zeros([N])
sigma = np.zeros([N])
acoeff = np.zeros([N])
bcoeff = np.zeros([N])

# Compute in the x, and z directions
for ind in range(0, domain.cpml):
    # From 0
    sigma[ind] = sig_max*dist[ind]**NP
    kappa[ind] = 1 + (k_max - 1) * dist[ind]**NP
    alpha[ind] = alpha_max * (1 - dist[ind])**NPA
    sigma[-(ind+1)] = sig_max*dist[ind]**NP
    kappa[-(ind+1)] = 1 + (k_max - 1) * dist[ind]**NP
    alpha[-(ind+1)] = alpha_max * (1 - dist[ind])**NPA
    bcoeff[ind] = np.exp(- (sigma[ind] / kappa[ind] + alpha[ind]) * electromag.dt / eps0 )
    bcoeff[-(ind+1)] = np.exp(- (sigma[-(ind+1)] / kappa[-(ind+1)] + alpha[-(ind+1)]) * electromag.dt / eps0)

# Compute the a-coefficients 
alpha[np.where(alpha < 0.0)] = 0.0
indices = np.where(sigma > 1.0e-6)
acoeff[indices] = sigma[indices] * (bcoeff[indices] - 1) / \
        ( ( kappa[indices] * sigma[indices] + kappa[indices] * alpha[indices] ) )

if half:
    sigma.tofile('sigma' + directory + '_half_cpml.dat')
    kappa.tofile('kappa' + directory + '_half_cpml.dat')
    alpha.tofile('alpha' + directory + '_half_cpml.dat')
    acoeff.tofile('acoef' + directory + '_half_cpml.dat')
    bcoeff.tofile('bcoef' + directory + '_half_cpml.dat')
else:
    sigma.tofile('sigma' + directory + '_cpml.dat')
    kappa.tofile('kappa' + directory + '_cpml.dat')
    alpha.tofile('alpha' + directory + '_cpml.dat')
    acoeff.tofile('acoef' + directory + '_cpml.dat')
    bcoeff.tofile('bcoef' + directory + '_cpml.dat')
