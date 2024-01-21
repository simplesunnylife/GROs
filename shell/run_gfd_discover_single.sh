#!/bin/sh

#mpirun -N 1 -n 1 -c 1 
./build/run_app --vfile ./dataset/test.v --efile ./dataset/test.e --application gfd_discover --pattern_v_file ./dataset/pattern_v.csv --pattern_e_file ./dataset/pattern_e.csv --out_prefix ./output_dpiso --directed --yaml_file  ${1}
