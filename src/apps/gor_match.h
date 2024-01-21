#ifndef EXAMPLES_ANALYTICAL_APPS_GOR_MATCH_H_
#define EXAMPLES_ANALYTICAL_APPS_GOR_MATCH_H_

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

#include "flags.h"
#include "gor_match/gor_match.h"
#include "timer.h"

namespace grape {

void Init() {
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
                    int fnum,
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
  std::shared_ptr<FRAG_T> fragment(new FRAG_T());
  std::vector<typename FRAG_T::internal_vertex_t> vertices;
  std::vector<typename FRAG_T::edge_t> edges;
  fragment->Init(comm_spec.fid(), vertices, edges);
  auto app = std::make_shared<APP_T>();
  timer_next("load application");
  auto worker = APP_T::CreateWorker(app, fragment);
  worker->Init(comm_spec, spec);
  timer_next("run algorithm");
  worker->Query(std::forward<Args>(args)...);
  timer_next("print output");
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
  auto spec = DefaultParallelEngineSpec();
  if (FLAGS_app_concurrency != -1) {
    spec.thread_num = FLAGS_app_concurrency;
  } else {
    spec = MultiProcessSpec(comm_spec, false);
  }
  int fnum = comm_spec.fnum();
  constexpr bool load_whole_graph = true;
  using GraphType = ImmutableEdgecutFragment<OID_T, VID_T, VDATA_T, EDATA_T,
                                            LoadStrategy::kLoadWholeGraph>;
  using AppType = GorMatch<GraphType>;
  std::string yaml_file = FLAGS_yaml_file;
  std::cout << "comm_spec: " << comm_spec.fid()
            << " fnum:" << comm_spec.fnum() << std::endl;
  CreateAndQuery<GraphType, AppType, std::string>(
      comm_spec, load_whole_graph, fnum, spec,
      yaml_file, comm_spec.fnum());
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

#endif  // EXAMPLES_ANALYTICAL_APPS_GOR_MATCH_H_
