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
    eps = read_dat("eps11.dat", 'Ey', domain)
    nx = domain.nx + 2*domain.cpml
    dx = domain.dx
elif directory == 'y':
    eps = read_dat("eps22.dat", 'Ey', domain)
    nx = domain.ny + 2*domain.cpml 
    dx = domain.dy 
else:
    eps = read_dat("eps33.dat", 'Ey', domain) 
    nx = domain.nz + 2*domain.cpml 
    dx = domain.dz 

# -----------------------------------------------------------------------------
# Define constants
NP = 2 
NPA = 2 
k_max = 1.1e1 # This is the value determined from the litterature. 
eps0 = 8.85418782e-12

# Change relative permittivity to absolute permittivity 
eps[:,:] = eps*eps0 

# Compute the distance along the absorbing boundary relative to the end of the 
# original model space. 
dist = dx * np.arange(0, domain.cpml)
if half:
    dist = dist + dx/2 

dist = dx*domain.cpml - dist
dist = dist/(dx*domain.cpml)

# Allocate space
alpha_max = 2*np.pi*eps*electromag.f0/10
sig_max = 0.7 * (NP+1) / (dx * np.sqrt(eps) )

kappa = np.ones([domain.nz, domain.ny, domain.nx])
alpha = np.zeros([domain.nz, domain.ny, domain.nx])
sigma = np.zeros([domain.nz, domain.ny, domain.nx])
acoeff = np.zeros([domain.nz, domain.ny, domain.nx])
bcoeff = np.zeros([domain.nz, domain.ny, domain.nx])

# Protrude the values for epsilon in the y-direction
eps = np.repeat(eps[:,np.newaxis,:], domain.ny, axis = 1 )

# Compute in the x, and z directions
for ind in range(0, domain.cpml):
    # From 0
    sigma[ind,:,:] = np.array([sig_max[ind,:]*dist[-ind]**NP]*domain.ny)
    kappa[ind,:,:] = 1 + (k_max - 1) * dist[-ind]**NP
    alpha[ind,:,:] = np.array([alpha_max[ind,:] * (1 - dist[-ind])**NPA]*domain.ny)
    sigma[:,:,ind] = np.array([sig_max[:,ind]*dist[-ind]**NP]*domain.ny).T
    kappa[:,:,ind] = 1 + (k_max - 1) * dist[-ind]**NP 
    alpha[:,:,ind] = np.array([alpha_max[:,ind] * (1 - dist[-ind])**NPA]*domain.ny).T
    sigma[:,ind,:] = sig_max*dist[-ind]**NP
    kappa[:,ind,:] = 1 + (k_max - 1) * dist[-ind]**NP
    alpha[:,ind,:] = alpha_max * (1 - dist[-ind])**NPA  
    bcoeff[ind,:,:] = np.exp(- (sigma[ind,:,:] / kappa[ind,:,:] + alpha[ind,:,:]) * electromag.dt / eps[ind,:,:] )
    bcoeff[:,ind,:] = np.exp(- (sigma[:,ind,:] / kappa[:,ind,:] + alpha[:,ind,:]) * electromag.dt / eps[:,ind,:] )
    bcoeff[:,:,ind] = np.exp(- (sigma[:,:,ind] / kappa[:,:,ind] + alpha[:,:,ind]) * electromag.dt / eps[:,:,ind] )
    #From nx 
    sigma[-ind,:,:] = np.array([sig_max[-ind,:]*dist[-ind]**NP]*domain.ny)
    kappa[-ind,:,:] = 1 + (k_max - 1) * dist[-ind]**NP
    alpha[-ind,:,:] = np.array([alpha_max[-ind,:] * (1 - dist[-ind])**NPA]*domain.ny)
    sigma[:,:,-ind] = np.array([sig_max[:,-ind]*dist[-ind]**NP]*domain.ny).T
    kappa[:,:,-ind] = 1 + (k_max - 1) * dist[-ind]**NP 
    alpha[:,:,-ind] = np.array([alpha_max[:,-ind] * (1 - dist[-ind])**NPA]*domain.ny).T
    sigma[:,-ind,:] = sig_max*dist[-ind]**NP
    kappa[:,-ind,:] = 1 + (k_max - 1) * dist[-ind]**NP
    alpha[:,-ind,:] = alpha_max * (1 - dist[-ind])**NPA
    bcoeff[-ind,:,:] = np.exp(- (sigma[-ind,:,:] / kappa[-ind,:,:] + alpha[-ind,:,:]) * electromag.dt / eps[-ind,:,:])
    bcoeff[:,-ind,:] = np.exp(- (sigma[:,-ind,:] / kappa[:,-ind,:] + alpha[:,-ind,:]) * electromag.dt / eps[:,-ind,:])
    bcoeff[:,:,-ind] = np.exp(- (sigma[:,:,-ind] / kappa[:,:,-ind] + alpha[:,:,-ind]) * electromag.dt / eps[:,:,-ind])

    
# Compute the a-coefficients 
alpha[np.where(alpha < 0.0)] = 0.0
indices = np.where(sigma > 1.0e-6)
acoeff[indices] = sigma[indices] * (bcoeff[indices] - 1) / ( (sigma[indices] + kappa[indices] * alpha[indices] ) )/ kappa[indices]

if dim != 2.5:
    sigma = sigma[:,int(domain.ny/2),:]
    kappa = kappa[:,int(domain.ny/2),:]
    alpha = alpha[:,int(domain.ny/2),:]
    acoeff = acoeff[:,int(domain.ny/2),:]
    bcoeff = bcoeff[:,int(domain.ny/2),:] 

# Save the results to a fortran binary
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
