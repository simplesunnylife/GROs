#!/bin/sh

# mpirun -n 1 -cpus-per-proc 4 
./build/run_app  --vfile ./dataset/test.v --efile ./dataset/test.e --application select_edge --pattern_v_file ./dataset/pattern_v.csv --pattern_e_file ./dataset/pattern_e.csv --out_prefix ./output_dpiso --directed --yaml_file  ${1}
