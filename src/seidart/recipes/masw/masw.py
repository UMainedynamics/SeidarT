import numpy as np 
from seidart.routines import prjbuild, prjrun, sourcefunction
from seidart.routines.arraybuild import Array
from seidart.visualization.im2anim import build_animation


prjfile = 'masw.prj' 
rcxfile = 'receivers.xyz'

## Initiate the model and domain objects
dom, mat, seis, em = prjrun.domain_initialization(prjfile)

## Compute the stiffness coefficients
# prjrun.status_check(
#     em, mat, dom, prjfile, seismic = True, append_to_prjfile = True
# )
prjrun.status_check(
    em, mat, dom, prjfile, seismic = True, append_to_prjfile = False
)

timevec, fx, fy, fz, srcfn = sourcefunction(seis, 1e7, 'gaus1', 's')

prjrun.runseismic(seis, mat, dom)


# Since this is MASW we are only going to use the vertical component
array_vz = Array('Vz', prjfile, rcxfile)
array_vz.gain = int(seis.time_steps/3)
array_vz.exaggeration = 0.1
array_vz.sectionplot(
    plot_complex = False
)


n_rcx = array_vz.receiver_xyz.shape[0]
time_array = np.tile(
    np.arange(seis.time_steps) * seis.dt, ( n_rcx , 1)
).T

# Calculate the distance of the receiver from the source. 
receiver_distances = 
receiver_distances = np.sqrt(
    np.sum( 
        (array_vz.receiver_xyz * dom.dz - array_vz.source)**2, axis = 1
    )
)
# We want distance values to the left of the source to be negative
receiver_distances = receiver_distances * np.sign( (array_vz.receiver_xyz * dom.dz - array_vz.source)[:,0])

receiver_distances = receiver_distances.repeat(
    int(seis.time_steps)
).reshape(n_rcx, int(seis.time_steps)).T


plot_phase_picks(array_vz.timeseries, time, receiver_distances, labels, np.unique(labels), gain = 301, exaggeration = 600)


build_animation(
        prjfile, 
        'Vz', 10, 20, 0.3, 
        is_single_precision = True,
)
