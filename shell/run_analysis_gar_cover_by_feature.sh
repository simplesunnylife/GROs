#!/bin/sh

mpirun -n 1 -cpus-per-proc 1 ./build/run_app --vfile ./dataset/test.v --efile ./dataset/test.e --application gar_cover_by_feature --pattern_v_file ./dataset/pattern_v.csv --pattern_e_file ./dataset/pattern_e.csv --out_prefix ./output_dpiso --directed --yaml_file  ${1}

./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_4_supp_50_conf_0.5.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_4_supp_50_conf_0.5.log &

./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.1.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.1.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.2.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.2.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.2.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.2.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.3.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.3.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.3.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.3.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.5.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.5.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.1.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.1.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.1.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.1.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.2.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.2.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.2.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.2.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.3.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.3.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.3.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.3.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.5.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.5.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.1.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.1.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.1.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.1.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.2.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.2.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.2.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.2.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.3.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.3.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.3.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.3.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.5.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.5.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.1.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.1.log &

./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.4.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.4.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.4.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.4.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.4.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.4.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.4.log &

./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.6.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.6.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.6.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.6.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.6.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.6.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.6.log &

./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.7.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.7.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.7.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.7.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.7.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.7.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.7.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.7.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.7.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.7.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.7.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.7.log &

./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.8.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_100_conf_0.8.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.8.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_100_conf_0.8.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.8.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_150_conf_0.8.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.8.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_150_conf_0.8.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.8.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_3_supp_50_conf_0.8.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.8.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_2_supp_50_conf_0.8.log &

./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_5_supp_50_conf_0.4.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_5_supp_50_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_5_supp_50_conf_0.5.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_5_supp_50_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_6_supp_50_conf_0.4.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_6_supp_50_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_movielen_yaml/gcr_discover_movielen_10_6_supp_50_conf_0.5.yaml &> ./gcr_movielen_yaml/gcr_discover_movielen_10_6_supp_50_conf_0.5.log &


./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_100_conf_0.4.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_100_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_100_conf_0.4.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_100_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_150_conf_0.4.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_150_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_150_conf_0.4.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_150_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_50_conf_0.4.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_50_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_50_conf_0.4.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_50_conf_0.4.log &

./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_100_conf_0.5.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_100_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_100_conf_0.5.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_100_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_150_conf_0.5.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_150_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_150_conf_0.5.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_150_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_50_conf_0.5.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_50_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_50_conf_0.5.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_50_conf_0.5.log &

./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_100_conf_0.6.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_100_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_100_conf_0.6.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_100_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_150_conf_0.6.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_150_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_150_conf_0.6.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_150_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_50_conf_0.6.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_50_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_50_conf_0.6.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_50_conf_0.6.log &

./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_50_conf_0.4.yaml  &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_50_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_100_conf_0.4.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_100_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_150_conf_0.4.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_150_conf_0.4.log &

./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_50_conf_0.4.yaml  &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_50_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_100_conf_0.4.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_100_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_150_conf_0.4.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_150_conf_0.4.log &

./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_5_supp_50_conf_0.4.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_5_supp_50_conf_0.4.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_6_supp_50_conf_0.4.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_6_supp_50_conf_0.4.log &

./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_50_conf_0.5.yaml  &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_50_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_100_conf_0.5.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_100_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_150_conf_0.5.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_150_conf_0.5.log &

./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_50_conf_0.5.yaml  &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_50_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_100_conf_0.5.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_100_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_150_conf_0.5.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_150_conf_0.5.log &

./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_5_supp_50_conf_0.5.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_5_supp_50_conf_0.5.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_6_supp_50_conf_0.5.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_6_supp_50_conf_0.5.log &

./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_50_conf_0.6.yaml  &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_50_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_100_conf_0.6.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_100_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_150_conf_0.6.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_2_supp_150_conf_0.6.log &

./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_50_conf_0.6.yaml  &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_50_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_100_conf_0.6.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_100_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_150_conf_0.6.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_3_supp_150_conf_0.6.log &

./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_5_supp_50_conf_0.6.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_5_supp_50_conf_0.6.log &
./shell/run_gar_discover_single_thread.sh ./gcr_ciao_yaml/gcr_discover_ciao_10_6_supp_50_conf_0.6.yaml &> ./gcr_ciao_yaml/gcr_discover_ciao_10_6_supp_50_conf_0.6.log &
