import numpy as np 
from seidart.routines import prjbuild, prjrun, sourcefunction
from seidart.routines.arraybuild import Array
from seidart.visualization.im2anim import build_animation


prjfile = 'single_source.prj' 
rcxfile = 'receivers.xyz'

## Initiate the model and domain objects
dom, mat, seis, em = prjrun.domain_initialization(prjfile)

## Compute the permittivity coefficients
# prjrun.status_check(em, mat, dom, prjfile, seismic = False)
prjrun.status_check(
    em, mat, dom, prjfile, seismic = False, append_to_prjfile = True
)

timevec, fx, fy, fz, srcfn = sourcefunction(em, 10, 'gaus1', 'e')

complex_values = False
prjrun.runelectromag(em, mat, dom, use_complex_equations = complex_values)

array_ex = Array('Ex', prjfile, rcxfile, is_complex = complex_values)
array_ex.gain = int(em.time_steps/3)
array_ex.exaggeration = 0.1 
array_ex.sectionplot(
    plot_complex = False
)
build_animation(
        prjfile, 
        'Ex', 10, 10, 0.3, 
        is_complex = complex_values, 
        is_single_precision = True
)


array_ez = Array('Ez', prjfile, rcxfile, is_complex = complex_values)
array_ez.gain = int(em.time_steps/3)
array_ez.exaggeration = 0.1
array_ez.sectionplot(
    plot_complex = False
)
build_animation(
        prjfile, 
        'Ex', 10, 10, 0.3, 
        is_complex = complex_values, 
        is_single_precision = True,
        plottype = 'energy_density'
)