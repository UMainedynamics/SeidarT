
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
cd source
$f2pypath/f2py3 -c --fcompiler=gnu95 -m emfdtd2d emFDTD2D.f95
$f2pypath/f2py3 -c --fcompiler=gnu95 -m seismicfdtd2d seismicFDTD2D.f95
$f2pypath/f2py3 -c --fcompiler=gnu95 -m emfdtd25d emFDTD25D.f95
$f2pypath/f2py3 -c --fcompiler=gnu95 -m seismicfdtd25d seismicFDTD25D.f95
$f2pypath/f2py3 -c --fcompiler=gnu95 -m orientsynth orientsynth.f95
mv *.so ../bin
cd ..

# Synthetic microstructure
# cd materials
# $f2pypath/f2py3 -c --fcompiler=gnu95 -m orientsynth orientsynth.f95
# mv *.so ../bin
# cd ..

# --------------------------- Create the executables --------------------------
# Start with the python scripts
cp source/routines/prjbuild.py bin/prjbuild
cp source/routines/sourcefunction.py bin/sourcefunction
cp source/routines/orientation_tensor.py bin/orientation_tensor

# Move the visualization tools
cp source/visualization/rcxdisplay.py bin/rcxdisplay
cp source/visualization/im2anim.py bin/im2anim
cp source/visualization/vtkbuild.py bin/vtkbuild 
cp source/visualization/wiggleplot.py bin/wiggleplot 
cp source/visualization/imgen.py bin/imgen.py # The generalized image functions module
cp source/visualization/imvector.py bin/imvector
cp source/visualization/vectoranim.py bin/vectoranim
cp source/visualization/implot.py bin/implot
cp source/visualization/slice25d.py bin/slice25d
cp source/visualization/rcxgen.py bin/rcxgen

# Change them to executables
chmod +x bin/prjbuild \
        bin/sourcefunction \
        bin/rcxdisplay \
        bin/wiggleplot \
        bin/im2anim \
        bin/orientation_tensor \
        bin/vtkbuild \
        bin/imvector \
        bin/vectoranim \
        bin/implot \
        bin/slice25d \
        bin/rcxgen


# Now do the bash scripts
cp survey_wrappers/common_offset bin/common_offset
cp survey_wrappers/common_midpoint bin/common_midpoint

chmod +x bin/common_offset bin/common_midpoint

# ---------------- Move all other required files to bin folder ----------------
cp source/routines/prjrun.py bin/
cp source/routines/arraybuild.py bin/
cp source/routines/material_functions.py bin/
cp source/routines/definitions.py bin/