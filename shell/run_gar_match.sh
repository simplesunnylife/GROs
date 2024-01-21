#!/bin/sh

#mpirun -n 1 -cpus-per-proc 4 
./build/gar_match --yaml_file  ${1}
