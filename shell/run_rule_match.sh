#!/bin/sh

# mpirun -n 2 -cpus-per-proc 2 
./build/rule_match --yaml_file  ${1}
