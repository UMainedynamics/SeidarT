#!/bin/bash
set -e -x

# Define the directory to place the built wheels.
WHEELHOUSE=/io/wheelhouse

# Upgrade pip and install build dependencies in all Python environments.
for PYBIN in /opt/python/cp3*/bin; do
    "${PYBIN}/pip" install --upgrade pip setuptools wheel
    # Explicitly install numpy, though it should be picked up from pyproject.toml
    "${PYBIN}/pip" install numpy
done

# Compile wheels. Adjust the path to your package if necessary.
for PYBIN in /opt/python/cp3*/bin; do
    "${PYBIN}/pip" wheel /io/ --no-deps -w ${WHEELHOUSE}
done

# Repair wheels using auditwheel. This makes the wheels manylinux compatible.
for whl in ${WHEELHOUSE}/*.whl; do
    auditwheel repair "$whl" --plat manylinux2014_x86_64 -w ${WHEELHOUSE}
done

# (Optional) Adjust permissions of the wheels.
chmod -R a+rw ${WHEELHOUSE}
