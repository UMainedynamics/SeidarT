from setuptools import setup, find_packages
import os 
import shutil 
import platform
from setuptools.command.install import install

class CustomInstallCommand(install):
    def run(self):
        install.run(self)  # Run the original install command
        self._post_install()
    
    def _post_install(self):
        env_bin_path = os.path.join(os.getenv('CONDA_PREFIX'), 'bin')
        if os.path.exists(env_bin_path):
            system = platform.system().lower()
            arch = platform.machine().lower()
            
            binary_name = f'seidartfdtd-{system}-{arch}'
            if system == 'darwin' and arch == 'arm64':
                binary_name = 'seidartfdtd-darwin-arm64'
            elif system == 'darwin' and arch == 'x86_64':
                binary_name = 'seidartfdtd-darwin-x86_64'
            elif system == 'linux' and arch == 'x86_64':
                binary_name = 'seidartfdtd-linux-x86_64'
            else:
                raise RuntimeError("Unsupported platform {}-{}".format(system,arch))
            
            src_binary_path = os.path.join('src', 'seidart', 'binaries', binary_name)
            dest_binary_path = os.path.join(env_bin_path, 'seidartfdtd')
            shutil.copy(src_binary_path, dest_binary_path)


setup(
    name='seidart',
    version='3.0.0',
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
    keywords=['wave propagation', 'seismic', 'radar', 'snow', 'ice'],
    cmdclass={
        'install': CustomInstallCommand,
    }
)