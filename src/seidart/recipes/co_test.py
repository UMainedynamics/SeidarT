import numpy as np 
from seidart.routines import prjbuild, prjrun, sourcefunction
from seidart.simulations.common_offset import CommonOffset
from seidart.visualization.im2anim import build_animation


prjfile = 'glacial_water1.prj' 
rcxfile = 'common_offset_receivers.xyz' #'../../src/seidart/recipes/receivers.xyz'
srcfile = 'sources.xyz'


complex = False

co = CommonOffset(
    srcfile, 
    'Ex', 
    prjfile, 
    rcxfile, 
    receiver_indices = False, 
    is_complex = complex,
    single_precision = True
)

co.co_run(seismic = False)

co.gain = 800
co.exaggeration = 0.05
co.sectionplot(plot_complex = complex)

