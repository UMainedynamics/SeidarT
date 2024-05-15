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
build_animation(
        prjfile, 
        'Vz', 10, 20, 0.3, 
        is_single_precision = True,
)
