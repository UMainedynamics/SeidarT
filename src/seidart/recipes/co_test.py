from seidart.simulations.common_offset import CommonOffset
from seidart.visualization.im2anim import build_animation

#
print('Import Modules')
print()
print('The class definition is in the common_offset.py script in the simulations module.')
print('from seidart.simulations.common_offset import CommonOffset')
print()
print()
print('Define the project, receiver, and source file')
print("prjfile = 'glacial_water1.prj'")
print("rcxfile = 'common_offset_receivers.xyz'")
print("srcfile = 'sources.xyz'")
print()
print("Define other inputs we would like for the CommonOffset class")
print("channel = 'Ex'")
print("complex = False")
print()
print('Initialize the CommonOffset object.')
print(
    '''
    co = CommonOffset(
        srcfile, 
        channel, 
        prjfile, 
        rcxfile, 
        receiver_indices = False, 
        is_complex = complex,
        single_precision = True
    )
    '''
)

prjfile = 'glacial_water1.prj' 
rcxfile = 'common_offset_receivers.xyz' #'../../src/seidart/recipes/receivers.xyz'
srcfile = 'sources.xyz'

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

