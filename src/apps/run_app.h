/** Copyright 2020 Alibaba Group Holding Limited.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef EXAMPLES_ANALYTICAL_APPS_RUN_APP_H_
#define EXAMPLES_ANALYTICAL_APPS_RUN_APP_H_

#include <gflags/gflags.h>
#include <gflags/gflags_declare.h>
#include <glog/logging.h>
#include <grape/fragment/immutable_edgecut_fragment.h>
#include <grape/fragment/loader.h>
#include <grape/grape.h>
#include <grape/util.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef GRANULA
#include "thirdparty/atlarge-research-granula/granula.hpp"
#endif

// #include "distributed_match/distributed_match.h"
// #include "dpiso/dpiso.h"
#include "flags.h"
#include "rule_discover/gar_discover.h"
#include "select_edge/select_edge.h"
#include "gar_exist/gar_exist.h"
#include "pattern_exist/pattern_exist.h"
#include "pattern_match/pattern_match.h"
#include "cluster_to_support/cluster_to_support.h"
#include "constant_rule_analysis/constant_rule_analysis.h"
#include "gar_cover_by_feature/gar_cover_by_feature.h"
#include "discovered_gfd_analysis/discovered_gfd_analysis.h"
#include "gar_overlap/gar_overlap.h"
#include "graph_sample/graph_sample.h"
#include "graph_convert/graph_convert.h"
#include "three_tuple_to_csv/three_tuple_to_csv.h"
#include "merge_graph/merge_graph.h"
#include "test/test.h"
#include "timer.h"
namespace grape {

void Init() {
  if (FLAGS_out_prefix.empty()) {
    LOG(FATAL) << "Please assign an output prefix.";
  }
  if (FLAGS_deserialize && FLAGS_serialization_prefix.empty()) {
    LOG(FATAL) << "Please assign a serialization prefix.";
  } else if (FLAGS_vfile.empty() || FLAGS_efile.empty()) {
    LOG(FATAL) << "Please assign input vertex/edge files.";
  }

  if (access(FLAGS_out_prefix.c_str(), 0) != 0) {
    mkdir(FLAGS_out_prefix.c_str(), 0777);
  }

  InitMPIComm();
  CommSpec comm_spec;
  comm_spec.Init(MPI_COMM_WORLD);
  if (comm_spec.worker_id() == kCoordinatorRank) {
    VLOG(1) << "Workers of libgrape-lite initialized.";
  }
}

void Finalize() {
  FinalizeMPIComm();
  VLOG(1) << "Workers finalized.";
}

template <typename FRAG_T, typename APP_T, typename... Args>
void CreateAndQuery(const CommSpec& comm_spec, const bool load_whole_graph,
                    const std::string& efile, const std::string& vfile,
                    const std::string& out_prefix, int fnum,
                    const ParallelEngineSpec& spec, Args... args) {
  assert((load_whole_graph &&
          FRAG_T::load_strategy == LoadStrategy::kLoadWholeGraph) ||
         (!load_whole_graph &&
          FRAG_T::load_strategy != LoadStrategy::kLoadWholeGraph &&
          FRAG_T::load_strategy != LoadStrategy::kNullLoadStrategy));
  timer_next("load graph");
  LoadGraphSpec graph_spec = DefaultLoadGraphSpec();
  graph_spec.set_directed(FLAGS_directed);
  graph_spec.set_rebalance(FLAGS_rebalance, FLAGS_rebalance_vertex_factor);
  if (FLAGS_deserialize) {
    graph_spec.set_deserialize(true, FLAGS_serialization_prefix);
  } else if (FLAGS_serialize) {
    graph_spec.set_serialize(true, FLAGS_serialization_prefix);
  }
  graph_spec.set_load_whole_graph(load_whole_graph);
  std::shared_ptr<FRAG_T> fragment;
  if (FLAGS_segmented_partition) {
    fragment = LoadGraph<FRAG_T, SegmentedPartitioner<typename FRAG_T::oid_t>>(
        efile, vfile, comm_spec, graph_spec);
  } else {
    fragment = LoadGraph<FRAG_T, HashPartitioner<typename FRAG_T::oid_t>>(
        efile, vfile, comm_spec, graph_spec);
  }
  auto app = std::make_shared<APP_T>();
  timer_next("load application");
  auto worker = APP_T::CreateWorker(app, fragment);
  worker->Init(comm_spec, spec);
  timer_next("run algorithm");
  worker->Query(std::forward<Args>(args)...);
  timer_next("print output");

  std::ofstream ostream;
  std::string output_path =
      grape::GetResultFilename(out_prefix, fragment->fid());
  ostream.open(output_path);
  worker->Output(ostream);
  ostream.close();
  worker->Finalize();
  timer_end();
  VLOG(1) << "Worker-" << comm_spec.worker_id() << " finished: " << output_path;
}

template <typename OID_T, typename VID_T, typename VDATA_T, typename EDATA_T>
void Run() {
  CommSpec comm_spec;
  comm_spec.Init(MPI_COMM_WORLD);

  bool is_coordinator = comm_spec.worker_id() == kCoordinatorRank;
  timer_start(is_coordinator);
#ifdef GRANULA
  std::string job_id = FLAGS_jobid;
  granula::startMonitorProcess(getpid());
  granula::operation grapeJob("grape", "Id.Unique", "Job", "Id.Unique");
  granula::operation loadGraph("grape", "Id.Unique", "LoadGraph", "Id.Unique");
  if (comm_spec.worker_id() == kCoordinatorRank) {
    std::cout << grapeJob.getOperationInfo("StartTime", grapeJob.getEpoch())
              << std::endl;
    std::cout << loadGraph.getOperationInfo("StartTime", loadGraph.getEpoch())
              << std::endl;
  }

  granula::linkNode(job_id);
  granula::linkProcess(getpid(), job_id);
#endif

#ifdef GRANULA
  if (comm_spec.worker_id() == kCoordinatorRank) {
    std::cout << loadGraph.getOperationInfo("EndTime", loadGraph.getEpoch())
              << std::endl;
  }
#endif
  // FIXME: no barrier apps. more manager? or use a dynamic-cast.
  std::string efile = FLAGS_efile;
  std::string vfile = FLAGS_vfile;
  std::string out_prefix = FLAGS_out_prefix;
  auto spec = DefaultParallelEngineSpec();
  if (FLAGS_app_concurrency != -1) {
    spec.thread_num = FLAGS_app_concurrency;
  } else {
    spec = MultiProcessSpec(comm_spec, false);
  }
  int fnum = comm_spec.fnum();
  std::string name = FLAGS_application;
  if (name.find("test") != std::string::npos) {
    constexpr bool load_whole_graph = false;
    using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
                                               LoadStrategy::kBothOutIn>;
    using AppType = Test<GraphType>;
    std::string yaml_file = FLAGS_yaml_file;
    CreateAndQuery<GraphType, AppType, std::string>(
        comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
        yaml_file);
  // } else if (name.find("dpiso") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = DPISO<GraphType>;
  //   std::string pattern_v_file = FLAGS_pattern_v_file;
  //   std::string pattern_e_file = FLAGS_pattern_e_file;
  //   CreateAndQuery<GraphType, AppType, std::string, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       pattern_v_file, pattern_e_file);
  // } else if (name.find("distributedmatch") != std::string::npos) {
  //   constexpr bool load_whole_graph = false;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kBothOutIn>;
  //   using AppType = DistributedMatch<GraphType>;
  //   std::string pattern_v_file = FLAGS_pattern_v_file;
  //   std::string pattern_e_file = FLAGS_pattern_e_file;
  //   CreateAndQuery<GraphType, AppType, std::string, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       pattern_v_file, pattern_e_file);
  // } else if (name.find("gar_cover_by_feature") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = GarCoverByFeature<GraphType>;
  //   std::string yaml_file = FLAGS_yaml_file;
  //   std::cout << "comm_spec: " << comm_spec.fid()
  //             << " fnum:" << comm_spec.fnum() << std::endl;
  //   CreateAndQuery<GraphType, AppType, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       yaml_file, comm_spec.fnum());
  // } else if (name.find("gar_exist") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = GarExist<GraphType>;
  //   std::string yaml_file = FLAGS_yaml_file;
  //   std::cout << "comm_spec: " << comm_spec.fid()
  //             << " fnum:" << comm_spec.fnum() << std::endl;
  //   CreateAndQuery<GraphType, AppType, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       yaml_file, comm_spec.fnum());
  // } else if (name.find("pattern_exist") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = PatternExist<GraphType>;
  //   std::string yaml_file = FLAGS_yaml_file;
  //   std::cout << "comm_spec: " << comm_spec.fid()
  //             << " fnum:" << comm_spec.fnum() << std::endl;
  //   CreateAndQuery<GraphType, AppType, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       yaml_file, comm_spec.fnum());
  // } else if (name.find("pattern_match") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = PatternMatch<GraphType>;
  //   std::string yaml_file = FLAGS_yaml_file;
  //   std::cout << "comm_spec: " << comm_spec.fid()
  //             << " fnum:" << comm_spec.fnum() << std::endl;
  //   CreateAndQuery<GraphType, AppType, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       yaml_file, comm_spec.fnum());
  // } else if (name.find("discovered_gfd_analysis") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = DiscoveredGfdAnalysis<GraphType>;
  //   std::string yaml_file = FLAGS_yaml_file;
  //   std::cout << "comm_spec: " << comm_spec.fid()
  //             << " fnum:" << comm_spec.fnum() << std::endl;
  //   CreateAndQuery<GraphType, AppType, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       yaml_file, comm_spec.fnum());
  // } else if (name.find("cluster_to_support") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = ClusterToSupport<GraphType>;
  //   std::string yaml_file = FLAGS_yaml_file;
  //   std::cout << "comm_spec: " << comm_spec.fid()
  //             << " fnum:" << comm_spec.fnum() << std::endl;
  //   CreateAndQuery<GraphType, AppType, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       yaml_file, comm_spec.fnum());
  // } else if (name.find("graph_sample") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = GraphSample<GraphType>;
  //   std::string yaml_file = FLAGS_yaml_file;
  //   std::cout << "comm_spec: " << comm_spec.fid()
  //             << " fnum:" << comm_spec.fnum() << std::endl;
  //   CreateAndQuery<GraphType, AppType, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       yaml_file, comm_spec.fnum());
  // } else if (name.find("graph_convert") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = GraphConvert<GraphType>;
  //   std::string yaml_file = FLAGS_yaml_file;
  //   std::cout << "comm_spec: " << comm_spec.fid()
  //             << " fnum:" << comm_spec.fnum() << std::endl;
  //   CreateAndQuery<GraphType, AppType, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       yaml_file, comm_spec.fnum());
  // } else if (name.find("three_tuple_to_csv") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = ThreeTupleToCsv<GraphType>;
  //   std::string yaml_file = FLAGS_yaml_file;
  //   std::cout << "comm_spec: " << comm_spec.fid()
  //             << " fnum:" << comm_spec.fnum() << std::endl;
  //   CreateAndQuery<GraphType, AppType, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       yaml_file, comm_spec.fnum());
  // } else if (name.find("constant_rule_analysis") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = ConstantRuleAnalysis<GraphType>;
  //   std::string yaml_file = FLAGS_yaml_file;
  //   std::cout << "comm_spec: " << comm_spec.fid()
  //             << " fnum:" << comm_spec.fnum() << std::endl;
  //   CreateAndQuery<GraphType, AppType, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       yaml_file, comm_spec.fnum());
  // } else if (name.find("select_edge") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = SelectEdge<GraphType>;
  //   std::string yaml_file = FLAGS_yaml_file;
  //   std::cout << "comm_spec: " << comm_spec.fid()
  //             << " fnum:" << comm_spec.fnum() << std::endl;
  //   CreateAndQuery<GraphType, AppType, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       yaml_file, comm_spec.fnum());
  // } else if (name.find("merge_graph") != std::string::npos) {
  //   constexpr bool load_whole_graph = true;
  //   using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
  //                                              LoadStrategy::kLoadWholeGraph>;
  //   using AppType = MergeGraph<GraphType>;
  //   std::string yaml_file = FLAGS_yaml_file;
  //   std::cout << "comm_spec: " << comm_spec.fid()
  //             << " fnum:" << comm_spec.fnum() << std::endl;
  //   CreateAndQuery<GraphType, AppType, std::string>(
  //       comm_spec, load_whole_graph, efile, vfile, out_prefix, fnum, spec,
  //       yaml_file, comm_spec.fnum());
  } else {
    LOG(FATAL) << "No avaiable application named [" << name << "].";
  }
#ifdef GRANULA
  granula::operation offloadGraph("grape", "Id.Unique", "OffloadGraph",
                                  "Id.Unique");
#endif

#ifdef GRANULA
  if (comm_spec.worker_id() == kCoordinatorRank) {
    std::cout << offloadGraph.getOperationInfo("StartTime",
                                               offloadGraph.getEpoch())
              << std::endl;

    std::cout << grapeJob.getOperationInfo("EndTime", grapeJob.getEpoch())
              << std::endl;
  }

  granula::stopMonitorProcess(getpid());
#endif
}

}  // namespace grape

#endif  // EXAMPLES_ANALYTICAL_APPS_RUN_APP_H_
