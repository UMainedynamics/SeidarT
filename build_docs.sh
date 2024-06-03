#!/bin/bash 

cp README.rst docs/source 

cd ./docs
sphinx-apidoc -o ./source/ ../src -e 
make html 
cd ../
