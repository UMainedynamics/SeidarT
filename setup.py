from setuptools import setup, find_packages

setup(
    name='seidart',
    version='2.4.1',
    license='GPL-3.0-or-later',
    author='Steven Bernsen, Christopher Gerbi',
    author_email='steven.bernsen@maine.edu, christopher.gerbi@maine.edu',
    maintainer='Steven Bernsen, Christopher Gerbi',
    maintainer_email='steven.bernsen@maine.edu, christopher.gerbi@maine.edu',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,  # Ensures that non-Python files are included
    entry_points={
        'console_scripts': [
            'prjbuild=seidart.routines.prjbuild:main',
            'prjrun=seidart.routines.prjrun:main',
            'arraybuild=seidart.routines.arraybuild:main',
            'sourcefunction=seidart.routines.sourcefunction:main',
            'rcxdisplay=seidart.visualization.rcxdisplay:main',
            'im2anim=seidart.visualization.im2anim:build_animation',
        ],
    },
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'pandas',
        'mplstereonet',
        'pyevtk',
        'glob2',
        'dill',
        'seaborn'
    ],
    project_urls={
        'Homepage': 'https://github.com/UMainedynamics/SeidarT',
        'Documentation': 'https://umainedynamics.github.io/SeidarT/index.html',
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: OS Independent',
        'Framework :: Matplotlib',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Development Status :: 4 - Beta'
    ],
    keywords=['wave propagation', 'seismic', 'radar', 'snow', 'ice']
)