
#!/bin/bash

# Install script for SEIDART toolbox

# Anaconda has issues.
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    f2pypath=`which f2py3`
else 
    f2pypath=`which pythonw3`
fi
f2pypath=`dirname $f2pypath`

# Make sure we have a folder to put everything in
if [ ! -d "bin" ]; then
	mkdir bin
fi

# Compile the fortran code
#2D
cd fdtd
$f2pypath/f2py3 -c --fcompiler=gnu95 -m emfdtd2d emFDTD2D.f95
$f2pypath/f2py3 -c --fcompiler=gnu95 -m seismicfdtd2d seismicFDTD2D.f95
$f2pypath/f2py3 -c --fcompiler=gnu95 -m emfdtd25d emFDTD25D.f95
$f2pypath/f2py3 -c --fcompiler=gnu95 -m seismicfdtd25d seismicFDTD25D.f95
mv *.so ../bin
cd ..

# Synthetic microstructure
cd materials
$f2pypath/f2py3 -c --fcompiler=gnu95 -m orientsynth orientsynth.f95
mv *.so ../bin
cd ..

# --------------------------- Create the executables --------------------------
# Start with the python scripts
cp exe/prjbuild.py bin/prjbuild
cp exe/prjrun.py bin/prjrun
cp exe/sourcefunction.py bin/sourcefunction
cp materials/orientation_tensor.py bin/orientation_tensor
cp materials/emcpml.py bin/emcpml 
cp materials/seiscpml.py bin/seiscpml 

# Move the visualization tools
cp vis/arraybuild.py bin/arraybuild
cp vis/rcxdisplay.py bin/rcxdisplay
cp vis/im2anim.py bin/im2anim
cp vis/vtkbuild.py bin/vtkbuild 
cp vis/wiggleplot.py bin/wiggleplot 
cp vis/imgen.py bin/imgen.py # The generalized image functions module
cp vis/imvector.py bin/imvector
cp vis/vectoranim.py bin/vectoranim
cp vis/implot.py bin/implot
cp vis/slice25d.py bin/slice25d
cp vis/rcxgen.py bin/rcxgen

# Change them to executables
chmod +x bin/prjbuild \
        bin/prjrun \
        bin/sourcefunction \
        bin/arraybuild \
        bin/rcxdisplay \
        bin/wiggleplot \
        bin/im2anim \
        bin/orientation_tensor \
        bin/vtkbuild \
        bin/imvector \
        bin/vectoranim \
        bin/implot \
        bin/slice25d \
        bin/emcpml \
        bin/seiscpml \
        bin/rcxgen



# Now do the bash scripts
cp survey_wrappers/common_offset bin/common_offset
cp survey_wrappers/common_midpoint bin/common_midpoint

cp io/array2sac bin/array2sac

chmod +x bin/common_offset bin/common_midpoint bin/array2sac


# ---------------- Move all other required files to bin folder ----------------
cp materials/material_functions.py bin/material_functions.py
cp materials/definitions.py bin/definitions.py