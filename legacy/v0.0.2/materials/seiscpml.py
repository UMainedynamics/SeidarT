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
direction = ''.join(args.direction)
half = args.half

# -----------------------------------------------------------------------------

# prjfile = "isotropic.prj"
# direction = 'x'
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

# Allocate space
if direction == 'x':
    N = domain.nx
    dx = domain.dx
elif direction == 'y':
    N = domain.ny 
    dx = domain.dy 
else:
    N = domain.nz
    dx = domain.dz


# -----------------------------------------------------------------------------
# Define constants
NP = 2 
NPA = 2 
k_max = 1.25e1 # This is the value determined from the litterature. 
Rcoef = 0.0010

# Compute the distance along the absorbing boundary relative to the end of the 
# original model space. 
dist = dx * np.arange(0, domain.cpml)
if half:
    dist = dist + dx/2 

dist = dx*domain.cpml - dist
dist = dist/(dx*domain.cpml)

quasi_cp_max = 0.7* np.min([domain.dx, domain.dy, domain.dz]) / (2.0 * seismic.dt)
alpha_max = np.pi*seismic.f0
sig_max = - np.log(Rcoef) * (NP+1) * quasi_cp_max / (2.0 * domain.cpml )

kappa = np.ones([N])
alpha = np.zeros([N])
sigma = np.zeros([N])
acoeff = np.zeros([N])
bcoeff = np.zeros([N])

# Compute in the x, and z directions
for ind in range(0, domain.cpml):
    # From 0
    sigma[ind] = sig_max*dist[ind]**NP
    kappa[ind] = 1.0 + (k_max - 1.0) * dist[ind]**NP
    alpha[ind] = alpha_max * (1 - dist[ind])**NPA
    sigma[-(ind+1)] = sig_max*dist[ind]**NP
    kappa[-(ind+1)] = 1 + (k_max - 1) * dist[ind]**NP
    alpha[-(ind+1)] = alpha_max * (1 - dist[ind])**NPA
    bcoeff[-(ind+1)] = np.exp(- (sigma[-(ind+1)] / kappa[-(ind+1)] + alpha[-(ind+1)]) * seismic.dt)
    bcoeff[ind] = np.exp( - (sigma[ind] / kappa[ind] + alpha[ind]) * seismic.dt)

# Compute the a-coefficients 
alpha[np.where(alpha < 0.0)] = 0.0
indices = np.where(np.abs(sigma) > 1.0e-6)
acoeff[indices] = sigma[indices] * (bcoeff[indices] - 1) / \
        (kappa[indices] * sigma[indices] + kappa[indices] * alpha[indices] )

# Save the results to a fortran binary
if half:
    sigma.tofile('sigma' + direction + '_half_cpml.dat')
    kappa.tofile('kappa' + direction + '_half_cpml.dat')
    alpha.tofile('alpha' + direction + '_half_cpml.dat')
    acoeff.tofile('acoef' + direction + '_half_cpml.dat')
    bcoeff.tofile('bcoef' + direction + '_half_cpml.dat')
else:
    sigma.tofile('sigma' + direction + '_cpml.dat')
    kappa.tofile('kappa' + direction + '_cpml.dat')
    alpha.tofile('alpha' + direction + '_cpml.dat')
    acoeff.tofile('acoef' + direction + '_cpml.dat')
    bcoeff.tofile('bcoef' + direction + '_cpml.dat')
