from seidart.simulations.common_offset import CommonOffset
from seidart.visualization.im2anim import build_animation

prjfile = 'common_offset.prj' 
rcxfile = 'common_offset_receivers.xyz'
srcfile = 'common_offset_sources.xyz'

channel = 'Ex'
complex = False

co = CommonOffset(
    srcfile, 
    channel, 
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

