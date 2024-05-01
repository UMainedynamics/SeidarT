# Create the anaconda environment using the yml file. Nothing will
# happen if it exists
conda env create -f seidart-environment.yml
conda activate SeidarT

# Compile the fortran modules
cd seidart/fortran 
#f2py -c --fcompiler=gnu95 -m seismicfdtd2d seismicFDTD2D.f95 
#f2py -c --fcompiler=gnu95 -m seismicfdtd25d seismicFDTD25D.f95 
#f2py -c --fcompiler=gnu95 -m emfdtd2d emFDTD2D.f95 
#f2py -c --fcompiler=gnu95 -m emfdtd25d emFDTD25D.f95 
#f2py -c --fcompiler=gnu95 cpmlfdtd.f95 seismicFDTD2D.f95 seismicFDTD25D.f95 emFDTD2D.f95 emFDTD25D.f95 -m fdtdmodeling
f2py -c --fcompiler=gnu95 cpmlfdtd.f95 -m cpmlfdtd
cd ../../

# Install either the editable or the pypi version
pip install -e seidart
# pip install -i https://test.pypi.org/simple/ seidart
