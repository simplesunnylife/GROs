#!/bin/sh

#mpirun -n 1 -cpus-per-proc 4 
./build/prob_gar_chase --yaml_file  ${1}
