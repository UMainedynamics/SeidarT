#!/bin/bash

# Set the Fortran compiler (adjust as needed)
FC=gfortran
fortran_path="src/seidart/fortran" 

# Compile each Fortran module file
$FC -c $fortran_path/cpmlfdtd.f95 -J$fortran_path -o $fortran_path/cpmlfdtd.o
$FC -c $fortran_path/cpmlfdtd.f95 -J"builddir" -o builddir/cpmlfdtd.o

# You can add error checking here if desired
if [ $? -ne 0 ]; then
    echo "Module compilation failed"
    exit 1
fi

echo "Modules compiled successfully"
