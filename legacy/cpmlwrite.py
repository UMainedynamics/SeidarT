#!/usr/bin/env python3


"""
Compute either the cpml for either the seismic or the electromagnetic model 
"""

import numpy as np 
from scipy.io import FortranFile
from definitions import *


def cpmlcomp(domain, modobj, direction, half, seismic):
    # Domain is the object output from loadproject()
    # seismic is boolean
    # Account for the boundary layer thickness
    domain.nx = domain.nx + 2*domain.cpml
    domain.ny = domain.ny + 2*domain.cpml
    domain.nz = domain.nz + 2*domain.cpml
    dim = domain.dim
    domain.dim = 2
    # dim = domain.dim
    # domain.dim = 2
    # Define constants
    NP = 2 
    NPA = 2 
    k_max = 1.5e1

    # Load sigma_xx, sigma_yy, or sigma_zz since these are going to be the 
    # maximum values. 
    if seismic:
        if direction == 'x':
            eps = read_dat("c11.dat", 'Vx', domain)
            # nx = domain.nx + 2*domain.cpml
            dx = domain.dx
        elif dir == 'y':
            eps = read_dat("c22.dat", 'Vx', domain)
            # nx = domain.ny + 2*domain.cpml 
            dx = domain.dy 
        else:
            eps = read_dat("c33.dat", 'Vz', domain) 
            # nx = domain.nz + 2*domain.cpml 
            dx = domain.dz 
        
        quasi_cp_max = np.min([domain.dx, domain.dy, domain.dz]) / (2.0 * modobj.dt)
        alpha_max = np.pi*modobj.f0
        sig_max = np.log(Rcoef) * (NP+1) * quasi_cp_max / (2.0 * domain.cpml )
    else:
        eps0 = 8.85418782e-12
        if direction == 'x':
            eps = read_dat("sig11.dat", 'Ey', domain)
            # nx = domain.nx + 2*domain.cpml
            dx = domain.dx
        elif direction == 'y':
            eps = read_dat("sig22.dat", 'Ey', domain)
            # nx = domain.ny + 2*domain.cpml 
            dx = domain.dy 
        else:
            eps = read_dat("sig33.dat", 'Ey', domain) 
            # nx = domain.nz + 2*domain.cpml 
            dx = domain.dz 
        
        eps0 = 8.85418782e-12
        # Change relative permittivity to absolute permittivity 
        eps[:,:] = eps*eps0 
        alpha_max = 2*np.pi*eps*modobj.f0/10
        sig_max = 0.7 * (NP+1) / (dx * np.sqrt(eps) )
        eps = np.repeat(eps[:,np.newaxis,:], domain.ny, axis = 1 )

    # Compute the distance along the absorbing boundary relative to the end of the 
    # original model space. 
    dist = dx * np.arange(0, domain.cpml)
    if half:
        dist = dist + dx/2 

    dist = dx*domain.cpml - dist
    dist_norm = dist/(dx*domain.cpml)

    # Allocate space 
    kappa = np.ones([domain.nz, domain.ny, domain.nx])
    alpha = np.zeros([domain.nz, domain.ny, domain.nx])
    sigma = np.zeros([domain.nz, domain.ny, domain.nx])
    acoef = np.zeros([domain.nz, domain.ny, domain.nx])
    bcoef = np.zeros([domain.nz, domain.ny, domain.nx])

    # Compute the boundary
    for ind in range(0, domain.cpml):
        # From 0
        sigma[ind,:,:] = np.array([sig_max[ind,:]*dist[ind]**NP]*domain.ny)
        kappa[ind,:,:] = 1.0 + (k_max - 1.0) * dist[ind]**NP
        alpha[ind,:,:] = np.array([alpha_max[ind,:] * (1 - dist[ind])**NPA]*domain.ny)
        sigma[:,:,ind] = np.array([sig_max[:,ind]*dist[ind]**NP]*domain.ny).T
        kappa[:,:,ind] = 1 + (k_max - 1) * dist[ind]**NP 
        alpha[:,:,ind] = np.array([alpha_max[:,ind] * (1 - dist[ind])**NPA]*domain.ny).T
        sigma[:,ind,:] = sig_max*dist[ind]**NP
        kappa[:,ind,:] = 1 + (k_max - 1) * dist[ind]**NP
        alpha[:,ind,:] = alpha_max * (1 - dist[ind])**NPA
        
        #From nx 
        sigma[-ind,:,:] = np.array([sig_max[ind,:]*dist[ind]**NP]*domain.ny)
        kappa[-ind,:,:] = 1 + (k_max - 1) * dist[ind]**NP
        alpha[-ind,:,:] = np.array([alpha_max[ind,:] * (1 - dist[ind])**NPA]*domain.ny)
        sigma[:,:,-ind] = np.array([sig_max[:,-ind]*dist[-ind]**NP]*domain.ny).T
        kappa[:,:,-ind] = 1 + (k_max - 1) * dist[-ind]**NP 
        alpha[:,:,-ind] = np.array([alpha_max[:,-ind] * (1 - dist[-ind])**NPA]*domain.ny).T
        sigma[:,-ind,:] = sig_max*dist[-ind]**NP
        kappa[:,-ind,:] = 1 + (k_max - 1) * dist[-ind]**NP
        alpha[:,-ind,:] = alpha_max * (1 - dist[-ind])**NPA
    
    # Compute the b-coefficients 
    if seismic:
        for ind in range(0, domain.cpml):
            bcoef[ind,:,:] = np.exp(
                -sigma[ind,:,:] / kappa[ind,:,:] + alpha[ind,:,:]
            ) * modobj.dt 
            bcoef[:,ind,:] = np.exp(
                -sigma[:,ind,:] / kappa[:,ind,:] + alpha[:,ind,:]
            ) * modobj.dt 
            bcoef[:,:,ind] = np.exp(
                -sigma[:,:,ind] / kappa[:,:,ind] + alpha[:,:,ind]
            ) * modobj.dt 
    else:
        for ind in range(0, domain.cpml):
            bcoef[ind,:,:] = np.exp(
                -sigma[ind,:,:] / kappa[ind,:,:] + alpha[ind,:,:]
            ) * modobj.dt / eps[ind,:,:]
            bcoef[:,ind,:] = np.exp(
                -sigma[:,ind,:] / kappa[:,ind,:] + alpha[:,ind,:]
            ) * modobj.dt / eps[:,ind,:]
            bcoef[:,:,ind] = np.exp(
                -sigma[:,:,ind] / kappa[:,:,ind] + alpha[:,:,ind]
            ) * modobj.dt / eps[:,:,ind]
    
    alpha[np.where(alpha < 0.0)] = 0.0
    indices = np.where(sigma > 1.0e-6)
    acoef[indices] = sigma[indices] * (bcoef[indices] - 1) / \
        ( (sigma[indices] + kappa[indices] * alpha[indices] ) )/ kappa[indices]

    return(sigma, kappa, alpha, bcoef, acoef)
    

def cpmlwrite(domain, modobj, seismic):
    direction = ['x', 'y', 'z']
    nx = domain.nx 
    ny = domain.ny 
    nz = domain.nz
    for d in direction:
        for half in [True, False]:
            sigma, kappa, alpha, bcoef, acoef = cpmlcomp(
                domain, modobj, d, half, seismic
            )
            domain.nx = nx 
            domain.ny = ny 
            domain.nz = nz
            if domain.dim != 2.5:
                sigma = sigma[:,int(domain.ny/2),:]
                kappa = kappa[:,int(domain.ny/2),:]
                alpha = alpha[:,int(domain.ny/2),:]
                acoef = acoef[:,int(domain.ny/2),:]
                bcoef = bcoef[:,int(domain.ny/2),:] 
        
            if half:
                sigma.tofile('sigma' + d + '_half_cpml.dat')
                kappa.tofile('kappa' + d + '_half_cpml.dat')
                alpha.tofile('alpha' + d + '_half_cpml.dat')
                acoef.tofile('acoef' + d + '_half_cpml.dat')
                bcoef.tofile('bcoef' + d + '_half_cpml.dat')
            else:
                sigma.tofile('sigma' + d + '_cpml.dat')
                kappa.tofile('kappa' + d + '_cpml.dat')
                alpha.tofile('alpha' + d + '_cpml.dat')
                acoef.tofile('acoef' + d + '_cpml.dat')
                bcoef.tofile('bcoef' + d + '_cpml.dat')

