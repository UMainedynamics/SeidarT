from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
    config = Configuration(None, parent_package, top_path)

    # Define Fortran extensions
    fortran_sources = [
        'src/seidart/fortran/cpmlfdtd.f95',
        # 'src/seidart/fortran/orientsynth.f95'
    ]
    
    config.add_extension(
        name='seidart.fortran.cpmlfdtd',
        sources=fortran_sources[0],
        # extra_f90_compile_args=['m64'],
    )
    # config.add_extension(
    #     name='seidart.fortran.orientsynth',
    #     sources = fortran_sources[1],
    #     # extra_f90_compile_args=['m64'],
    # )
    
    return config

if __name__ == "__main__":
    setup(
        name='seidart',
        version='2.2.3',
        packages=[
            'seidart', 
            'seidart.fortran', 
            'seidart.routines', 
            'seidart.simulations', 
            'seidart.visualization'
        ],
        configuration=configuration,
        entry_points = {
            'console_scripts': [
                'prjbuild=seidart.routines.prjbuild:main',
                'prjrun=seidart.routines.prjrun:main',
                'arraybuild=seidart.routines.arraybuild:main',
                'sourcefunction=seidart.routines.sourcefunction:main',
                'rcxdisplay=seidart.visualization.rcxdisplay:main',
                'im2anim=seidart.visualiztion.im2anim:build_animation',
            ]
        },
        install_requires=[
            'numpy',
            'setuptools',
            'wheel',
            'matplotlib',
            'seaborn'
        ]
    )
