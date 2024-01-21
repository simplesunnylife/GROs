#!/bin/sh

#mpirun -N 1 -n 1 -c 1 
./build/gar_analysis --yaml_file  ${1}
