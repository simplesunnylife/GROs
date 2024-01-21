#!/bin/sh

mpirun -n 4 -cpus-per-proc 1 ./build/pattern_match --yaml_file  ${1}
