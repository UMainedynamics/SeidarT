from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

class BuildExt(build_ext):
    def build_extensions(self):
        super().build_extensions()

def configuration():
    config = {
        'name': 'seidart',
        'version': '2.2.3',
        'packages': [
            'seidart', 
            'seidart.fortran', 
            'seidart.routines', 
            'seidart.simulations', 
            'seidart.visualization'
        ],
        # 'entry_points': {
        #     'console_scripts': [
        #         'prjbuild=seidart.routines.prjbuild:main',
        #         'prjrun=seidart.routines.prjrun:main',
        #         'arraybuild=seidart.routines.arraybuild:main',
        #         'sourcefunction=seidart.routines.sourcefunction:main',
        #         'rcxdisplay=seidart.visualization.rcxdisplay:main',
        #         'im2anim=seidart.visualization.im2anim:build_animation',
        #     ]
        # },
        'install_requires': [
            'numpy',
            'setuptools',
            'wheel',
            'matplotlib',
            'seaborn'
        ],
        'ext_modules': [
            Extension(
                name='seidart.fortran.cpmlfdtd',
                sources=['src/seidart/fortran/cpmlfdtd.f95'],
                # extra_f90_compile_args=['-m64'],
            ),
            # Uncomment if needed
            # Extension(
            #     name='seidart.fortran.orientsynth',
            #     sources=['src/seidart/fortran/orientsynth.f95'],
            #     # extra_f90_compile_args=['-m64'],
            # ),
        ],
        'cmdclass': {'build_ext': BuildExt},
    }
    return config

if __name__ == "__main__":
    setup(**configuration())