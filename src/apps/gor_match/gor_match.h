#ifndef EXAMPLES_ANALYTICAL_APPS_GOR_MATCH_GOR_MATCH_H_
#define EXAMPLES_ANALYTICAL_APPS_GOR_MATCH_GOR_MATCH_H_

#include <grape/grape.h>
#include <omp.h>

#include <functional>

#include <yaml-cpp/yaml.h>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "include/gundam/tool/radius.h"

#include "gor_match/gor_match_context.h"

#include "gor/gor.h"
#include "gor/gor_match.h"
#include "gor/csv_gor.h"

#include "util/log.h"

#include "../timer.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

namespace grape {

template <typename FRAG_T>
class GorMatch : public ParallelAppBase<FRAG_T, GorMatchContext<FRAG_T>>,
                 public ParallelEngine {
 public:
  static constexpr LoadStrategy load_strategy = LoadStrategy::kLoadWholeGraph;

 private:
  inline std::pair<std::tuple<std::string,
                              std::string,
                              std::string,
                              std::string,
                              std::string>,
                   bool> GetGorPathInfo(const YAML::Node& gor_path_config) const {
    if (!gor_path_config["VFile"]) {
      util::Error( "cannot get gor v file" );
      return std::pair(std::tuple("","","","",""), false);
    }
    const std::string kGorPathVFile 
           = gor_path_config["VFile"].as<std::string>();

    if (!gor_path_config["EFile"]) {
      util::Error( "cannot get gor e file" );
      return std::pair(std::tuple("","","","",""), false);
    }
    const std::string kGorPathEFile  
           = gor_path_config["EFile"].as<std::string>();

    if (!gor_path_config["XFile"]) {
      util::Error( "cannot get gor x file" );
      return std::pair(std::tuple("","","","",""), false);
    }
    const std::string kGorPathXFile 
           = gor_path_config["XFile"].as<std::string>();

    if (!gor_path_config["YFile"]) {
      util::Error( "cannot get gor y file" );
      return std::pair(std::tuple("","","","",""), false);
    }
    const std::string kGorPathYFile  
           = gor_path_config["YFile"].as<std::string>();

    std::string gor_pivot_id_str = "";
    if (gor_path_config["PivotId"]) {
      gor_pivot_id_str = gor_path_config["PivotId"].as<std::string>();
    }

    return std::pair(std::tuple(kGorPathVFile,
                                kGorPathEFile, 
                                kGorPathXFile,
                                kGorPathYFile,
                                 gor_pivot_id_str), true);
  }

  static constexpr std::string_view kWorkerGorMatchPosfix = "@worker";

  static constexpr std::string_view kRulePosfix    = GorMatchContext<FRAG_T>::kRulePosfix;
  static constexpr std::string_view kRuleSetPosfix = GorMatchContext<FRAG_T>::kRuleSetPosfix;

  static constexpr std::string_view kSuppPosfix    = GorMatchContext<FRAG_T>::kSuppPosfix;
  static constexpr std::string_view kSuppSetPosfix = GorMatchContext<FRAG_T>::kSuppSetPosfix;
                  
  using    VertexIDType = typename GorMatchContext<FRAG_T>::   VertexIDType;
  using VertexLabelType = typename GorMatchContext<FRAG_T>::VertexLabelType;
  using      EdgeIDType = typename GorMatchContext<FRAG_T>::     EdgeIDType;
  using   EdgeLabelType = typename GorMatchContext<FRAG_T>::  EdgeLabelType;
  using    DataGraph    = typename GorMatchContext<FRAG_T>::DataGraph;
  using      Pattern    = typename GorMatchContext<FRAG_T>::  Pattern;
  using      GorType    = typename GorMatchContext<FRAG_T>::  GorType;

  using   PatternVertexHandle = typename GUNDAM::VertexHandle<  Pattern>::type;
  using DataGraphVertexHandle = typename GUNDAM::VertexHandle<DataGraph>::type;

  using CandidateSet = typename GorMatchContext<FRAG_T>::CandidateSet;

  static constexpr auto compare_by_vertex_id 
  =[](const DataGraphVertexHandle& vertex_handle_0,
      const DataGraphVertexHandle& vertex_handle_1) {
    return vertex_handle_0->id() 
         < vertex_handle_1->id();
  };

  // the fragment calculate confidence and deliever the file name of the gors
  // to other workers
  static constexpr int kTimmerFragID = 0;

  // from ordinary workers to kTimmerFragID
  //    info the kTimmerFragID how many process each worker have
  static constexpr std::string_view kInfoProcessPrefix = "#info_process";
  // from ordinary workers to kTimmerFragID
  //    inform that all support has been calcluated, it would be ok to calcluate the
  //    confidence
  static constexpr std::string_view kInfoFinishPrefix = "#info_finish";

  // ###########################################
  // ##   protocol without workload balance   ##
  // ###########################################
  // from kTimmerFragID to ordinary workers
  //    info the ordinary workers begin to run the matching
  static constexpr std::string_view kBeginMatchPrefix = "#begin_match";

  // ########################################
  // ##   protocol with workload balance   ##
  // ########################################
  // from kTimmerFragID to ordinary workers
  //    alloc pivot set to each worker
  static constexpr std::string_view kAllocPivotSetPrefix = "#alloc_pivot_set";
  // from ordinary workers to kTimmerFragID
  //    info how many pivots are processed
  static constexpr std::string_view kInfoStatusPrefix = "#info_status_set";

  // #################################
  // ##   protocol with multi gor   ##
  // #################################
  // from kTimmerFragID to ordinary workers
  //    info them to reset parameter
  static constexpr std::string_view kResetParameterPrefix = "#reset_parameter";

  static constexpr size_t kInitialPivotSetNum = 5;

 public:
  // specialize the templated worker.
  INSTALL_PARALLEL_WORKER(GorMatch<FRAG_T>, 
                          GorMatchContext<FRAG_T>,
                          FRAG_T)
  using vertex_t = typename fragment_t::vertex_t;

  /**
   * @brief Partial evaluation for GorMatch.
   *
   * @param frag
   * @param ctx
   * @param messages
   */
  void PEval(const fragment_t& frag, 
                    context_t& ctx,
            message_manager_t& messages) {
    
    messages.InitChannels(thread_num());
    auto begin = timer();
    timer_next("load graph");

    util::Info("yaml file: " + ctx.yaml_file_);

    YAML::Node config = YAML::LoadFile(ctx.yaml_file_);

    if (!config["DataGraphPath"]) {
      util::Error("cannot get data graph path");
      return;
    }
    YAML::Node data_graph_path_config = config["DataGraphPath"];

    if (!data_graph_path_config["VFile"]) {
      util::Error("cannot get data graph v file");
      return;
    }
    const std::string kDataGraphPathVFile 
          = data_graph_path_config["VFile"].as<std::string>();

    if (!data_graph_path_config["EFile"]) {
      util::Error("cannot get data graph e file");
      return;
    }
    const std::string kDataGraphPathEFile 
          = data_graph_path_config["EFile"].as<std::string>();

    if (config["TimeLogFile"]) {
      assert(!ctx.time_log_file_.is_open());
      const std::string kTimeLogFileName = config["TimeLogFile"].as<std::string>();
      if (ctx.fid_ == kTimmerFragID) {
        ctx.time_log_file_.open(kTimeLogFileName, std::ios::app);
        ctx.time_log_file_ << "time log of: " << ctx.yaml_file_ << std::endl;
        if (!ctx.time_log_file_.is_open()) {
          util::Error("open time log file failed: " 
                     + kTimeLogFileName);
          return;
        }
      }
    }

    if (!config["GorPath"]) {
      util::Error("cannot get gor path");
      return;
    }
    YAML::Node gor_path_config = config["GorPath"];

    if (gor_path_config.IsSequence()) {
      ctx.gor_path_set_.reserve(gor_path_config.size());
      for (int i = 0; i < gor_path_config.size(); i++) {
        auto [ gor_path_string, 
               gor_path_ret ] = GetGorPathInfo(gor_path_config[i]);
        if (!gor_path_ret) {
          util::Error("load gor path error!");
          return;
        }
        ctx.gor_path_set_.emplace_back(gor_path_string);
      }
    }
    else {
      auto [ gor_path_string, 
             gor_path_ret ] = GetGorPathInfo(gor_path_config);
      if (!gor_path_ret) {
        util::Error("load gor path error!");
        return;
      }
      ctx.gor_path_set_.emplace_back(gor_path_string);
    }

    if (ctx.gor_path_set_.size() == 0) {
      util::Error("empty gor set");
      return;
    }
    
    ctx.gor_set_.reserve(ctx.gor_path_set_.size());
    assert(ctx.gor_set_.empty());

    for (const auto& [kGorPathVFile, 
                      kGorPathEFile, 
                      kGorPathXFile, 
                      kGorPathYFile, 
                           kPivotId] : ctx.gor_path_set_) {

      GorType gor;

      auto read_gor_res = gor::ReadGOR(gor, kGorPathVFile, 
                                            kGorPathEFile, 
                                            kGorPathXFile, 
                                            kGorPathYFile);

      if (read_gor_res < 0) {
        util::Error("read gor error! v file: " + kGorPathVFile
                                 + " e file: " + kGorPathEFile
                                 + " x file: " + kGorPathXFile
                                 + " y file: " + kGorPathYFile);
        return;
      }

      ctx.gor_set_.emplace_back(std::move(gor), PatternVertexHandle());

      assert(!ctx.gor_set_.back().second);

      auto& gor_ref = ctx.gor_set_.back().first;

      if (kPivotId != "") {
        const VertexIDType kGorPivotId = GUNDAM::StringToDataType<VertexIDType>(kPivotId);

        Pattern rhs_pattern(gar::_gar_supp::LiteralPattern(gor_ref, true));

        if (!rhs_pattern.FindVertex(kGorPivotId)) {
          // cannot find the pivot vertex in rhs pattern
          std::string rhs_pattern_vertex_set;
          for (auto vertex_it = rhs_pattern.VertexBegin();
                   !vertex_it.IsDone();
                    vertex_it++) {
            rhs_pattern_vertex_set += " " + GUNDAM::ToString(vertex_it->id());
          }
          util::Error("cannot find the pivot vertex with id: " + GUNDAM::ToString(kGorPivotId)
                      + " in rhs pattern contains of vertex: " + rhs_pattern_vertex_set);
          return;
        }
        // gor_ref.pattern().FindVertex(kGorPivotId)
        //     != rhs_pattern.FindVertex(kGorPivotId)
        // need to find it from graph pattern again
        ctx.gor_set_.back().second = gor_ref.pattern().FindVertex(kGorPivotId);
      }
      else {
        // randomly select one from rhs_pattern
        Pattern rhs_pattern(gar::_gar_supp::LiteralPattern(gor_ref, true));

        std::vector<VertexIDType> rhs_gor_vertex_id_set;
        for (auto vertex_it = rhs_pattern.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          rhs_gor_vertex_id_set.emplace_back(vertex_it->id());
        }
        std::sort(rhs_gor_vertex_id_set.begin(),
                  rhs_gor_vertex_id_set.end());
        assert(!rhs_gor_vertex_id_set.empty());

        ctx.gor_set_.back().second = gor_ref.pattern().FindVertex(*rhs_gor_vertex_id_set.begin());
      }
      assert(ctx.gor_set_.back().second);
    }

    #ifndef NDEBUG
    for (const auto& [gor, pivot] : ctx.gor_set_) {
      assert(pivot);
      assert(pivot == gor.pattern().FindVertex(pivot->id()));
    }
    #endif // NDEBUG

    if (GUNDAM::ReadCSVGraph(ctx.data_graph_, 
                            kDataGraphPathVFile,
                            kDataGraphPathEFile) < 0) {
      util::Error("load data graph failed!");
      return;
    }
     
    ctx.is_incremental_ = false;
    std::vector<VertexIDType> delta_vertex_id_set;
    if (data_graph_path_config["DeltaVFile"]) {
      ctx.is_incremental_ = true;
      const std::string kDataGraphDeltaVFile
        = data_graph_path_config["DeltaVFile"].as<std::string>();

      try {
        std::ifstream delta_vertex_file(kDataGraphDeltaVFile);
        
        if (!delta_vertex_file.good()) {
          util::Error("delta vertex file: "
                    + kDataGraphDeltaVFile + " is not good!");
          return;
        }

        bool has_found_empty_line = false;
        while (delta_vertex_file) {
          std::string s;
          if (!std::getline( delta_vertex_file, s )) 
            break;
          if (std::all_of(s.begin(),s.end(), 
                          [](char c){ 
                            return std::isspace(static_cast<unsigned char>(c));
                          })) {
            // allow empty lines at the last
            has_found_empty_line = true;
            continue;
          }
          // util::Debug(s);
          if (has_found_empty_line) {
            util::Error("has an empty line(s) but not at the last!");
            return;
          }
          const VertexIDType kVertexId
              = GUNDAM::StringToDataType<VertexIDType>(s);
          delta_vertex_id_set.emplace_back(kVertexId);
        }
        delta_vertex_file.close();
      }
      catch (std::exception& e) {
        util::Error("read delta vertex file: " 
                  + kDataGraphDeltaVFile 
                  + " illegal! : "
                  + e.what());
        return;
      }

      if (delta_vertex_id_set.empty()) {
        util::Error("empty delta vertex set");
        return;
      }
    }
    
    bool unbalanced_work_load = false;
    std::vector<size_t> worker_set;
    if (config[ "WorkloadUnbalanceRatio" ]) {
      const int kWorkloadUnbalanceRatio
       = config["WorkloadUnbalanceRatio"].as<int>();
      if (kWorkloadUnbalanceRatio <= 0) {
        util::Error("illegal WorkloadUnbalanceRatio: "
           + std::to_string(kWorkloadUnbalanceRatio));
        return;
      }
      if (ctx.frag_num_ == 1) {
        util::Error("needs to be more than one worker!");
        return;
      }
      unbalanced_work_load = true;
      worker_set.reserve(kWorkloadUnbalanceRatio + ctx.frag_num_ - 1);
      for (int worker_id = 0; worker_id < ctx.frag_num_; worker_id++) {
        if (worker_id == kTimmerFragID) {
          for (int i = 0; i < kWorkloadUnbalanceRatio; i++) {
            worker_set.emplace_back(worker_id);
          }
          continue;
        }
        worker_set.emplace_back(worker_id);
      }
      assert(worker_set.size() == kWorkloadUnbalanceRatio + ctx.frag_num_ - 1);
    }

    ctx.work_load_balance_ = false;
    if (config["WorkloadBalance"]) {
      if (unbalanced_work_load) {
        util::Error("specified both WorkloadBalance and WorkloadUnbalanceRatio");
        return;
      }
      if (ctx.is_incremental_) {
        util::Error("does not support workload balance for incremental matching");
        return;
      }
      YAML::Node workload_balance_config = config["WorkloadBalance"]; 

      if (!workload_balance_config["BalancePeriod"]) {
        util::Error("does not specify balance period");
        return;
      }
      ctx.balance_period_ms_ = workload_balance_config["BalancePeriod"].as<int32_t>(); 
      
      if (ctx.balance_period_ms_ <= 0 ) {
        util::Error("Illegal BalancePeriod: "
                  + std::to_string(ctx.balance_period_ms_)
                  + " ms, needs > 0");
        return;
      }
      // in order not to send all vertex id to each worker and let they
      // find corresponding handle in data graph again, each holds the 
      // set of candidate for pivot

      if (!workload_balance_config["PivotedCandidateBalanceSize"]) {
        util::Error("does not specify pivoted candidate size to balance");
        return;
      }
      ctx.pivoted_candidate_balance_size_ 
      = workload_balance_config["PivotedCandidateBalanceSize"].as<int>();
      if (ctx.pivoted_candidate_balance_size_ <= 0) {
        util::Error("Illegal PivotedCandidateBalanceSize: " 
                  + std::to_string(ctx.pivoted_candidate_balance_size_)
                  + ", needs > 0");
        return;
      }
      ctx.work_load_balance_ = true;
    }

    ctx.time_limit_ = -1.0;
    if (config["TimeLimit"]) {
      ctx.time_limit_ = config["TimeLimit"].as<double>();
      util::Info("TimeLimit: " + std::to_string(ctx.time_limit_) + "s");
    }

    if (config["Partition"]) {
      // specify partition
      YAML::Node partition_config = config["Partition"];
      if (partition_config.size() != ctx.frag_num_) {
        util::Error("partition number miss matches the worker number!");
        return;
      }
      const std::string kVertexIdSetFile = partition_config[(int)ctx.fid_].as<std::string>();
      util::Info("vertex id set file: " + kVertexIdSetFile);
      // wenzhi:
      //    ToDo: 
      //      allow user to specify graph partition, each worker
      //      loads the vertex id set from kVertexIdSetFile, find the
      //      corresponding vertex handle for each id in the data graph, 
      //      add those corresponding vertex handles into 
      //      ctx.data_graph_pivot_range_
      if (ctx.is_incremental_) {
        // wenzhi: 
        //  needs to partition the delta vertex set 
        //  according to the partition at the same time
      }
      util::Error("current version does not support user specify partition now! to be completed!");
      return;
    }
    else {
      // does not specify the partition, randomly partition 
      std::vector<DataGraphVertexHandle> data_graph_vertex_handle_set;
      if (!ctx.is_incremental_) {
        // partition the data graph
        ctx.data_graph_pivot_range_.reserve((ctx.data_graph_.CountVertex() 
                                           / ctx.frag_num_) + 1);
        // based on vertex id, to make sure the partition at
        // each worker are the same
        for (auto vertex_it = ctx.data_graph_.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          data_graph_vertex_handle_set.emplace_back(vertex_it);
        }
      }
      else {
        // partition both the delta vertex set and the 
        // vertex set in data graph
        ctx.data_graph_pivot_range_.reserve(
                ((ctx.data_graph_.CountVertex() - delta_vertex_id_set.size())
                / ctx.frag_num_) + 1);

        assert(!delta_vertex_id_set.empty());
        std::sort(delta_vertex_id_set.begin(),
                  delta_vertex_id_set.end());
        for (const auto& vertex_id : delta_vertex_id_set) {
          ctx.delta_vertex_handle_set_.emplace_back(ctx.data_graph_.FindVertex(vertex_id));
          // should have found this vertex
          assert(ctx.delta_vertex_handle_set_.back());
        }
        assert(ctx.delta_vertex_handle_set_.size()
                    == delta_vertex_id_set .size());
        // should have been sorted by vertex id
        // std::sort(ctx.delta_vertex_handle_set_.begin(),
        //           ctx.delta_vertex_handle_set_. end (),
        //           compare_by_vertex_id);
        for (auto vertex_it = ctx.data_graph_.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          if (std::binary_search(delta_vertex_id_set.begin(),
                                 delta_vertex_id_set.end(),
                                 vertex_it->id())) {
            // this vertex is in the delta vertex set
            // since ctx.delta_vertex_handle_set_ is also sorted by 
            // vertex id, it should also find this vertex
            #ifndef NDEBUG
            DataGraphVertexHandle vertex_handle = vertex_it;
            assert(std::binary_search(ctx.delta_vertex_handle_set_.begin(),
                                      ctx.delta_vertex_handle_set_.end(),
                                      vertex_handle, compare_by_vertex_id));
            #endif // NDEBUG
            continue;
          }
          // this vertex is not contained in the delta vertex set
          data_graph_vertex_handle_set.emplace_back(vertex_it);
        }
        assert(data_graph_vertex_handle_set .size()
              + ctx.delta_vertex_handle_set_.size()
             == ctx.data_graph_.CountVertex());
      }
      std::sort(data_graph_vertex_handle_set.begin(),
                data_graph_vertex_handle_set. end (),
                compare_by_vertex_id);

      if (!ctx.work_load_balance_) {
        // without workload balance, divide the pivot candidate here
        for (size_t vertex_handle_idx = 0;
                    vertex_handle_idx < data_graph_vertex_handle_set.size();
                    vertex_handle_idx++) {
          const size_t kWorkerToAlloc = unbalanced_work_load?
                        worker_set[vertex_handle_idx % worker_set.size()]
                                 : vertex_handle_idx % ctx.frag_num_;
          if (kWorkerToAlloc != ctx.fid_) {
            // would not be considered in this worker
            continue;
          }
          // this vertex is in the fragment of this node
          // should be in the data graph
          assert(data_graph_vertex_handle_set[vertex_handle_idx]);
          assert(data_graph_vertex_handle_set[vertex_handle_idx]
              == ctx.data_graph_.FindVertex(data_graph_vertex_handle_set[vertex_handle_idx]->id()));
          ctx.data_graph_pivot_range_.emplace_back(data_graph_vertex_handle_set[vertex_handle_idx]);
        }
        assert(std::is_sorted(ctx.data_graph_pivot_range_.begin(),
                              ctx.data_graph_pivot_range_.end(),
                              compare_by_vertex_id));

        std::sort(ctx.data_graph_pivot_range_.begin(),
                  ctx.data_graph_pivot_range_.end());

        util::Info("without workload balance, worker: "
                 +  std::to_string(ctx.fid_) 
                 + " pivots vertexes set size: "
                 +  std::to_string(ctx.data_graph_pivot_range_.size()));

        if (ctx.is_incremental_) {
          util::Debug("ctx.delta_vertex_handle_set_.size(): "
                     + std::to_string(ctx.delta_vertex_handle_set_.size()));
          util::Debug("ctx.partition_delta_vertex_handle_set_.size(): "
                     + std::to_string(ctx.partition_delta_vertex_handle_set_.size()));
          // also need to partition the delta vertex set
          for (size_t vertex_handle_idx = 0;
                      vertex_handle_idx < ctx.delta_vertex_handle_set_.size();
                      vertex_handle_idx++) {
            const size_t kWorkerToAlloc = unbalanced_work_load?
                         worker_set[vertex_handle_idx % worker_set.size()]
                                  : vertex_handle_idx % ctx.frag_num_;
            if (kWorkerToAlloc != ctx.fid_) {
              // would not be considered in this worker
              continue;
            }
            // this vertex is in the fragment of this node
            // should be in the data graph
            assert(ctx.delta_vertex_handle_set_[vertex_handle_idx]);
            assert(ctx.delta_vertex_handle_set_[vertex_handle_idx]
                == ctx.data_graph_.FindVertex(ctx.delta_vertex_handle_set_[vertex_handle_idx]->id()));
            ctx.partition_delta_vertex_handle_set_
               .emplace_back(ctx.delta_vertex_handle_set_[vertex_handle_idx]);
          }
          util::Debug("ctx.partition_delta_vertex_handle_set_.size(): "
                     + std::to_string(ctx.partition_delta_vertex_handle_set_.size()));
          assert(std::is_sorted(ctx.partition_delta_vertex_handle_set_.begin(),
                                ctx.partition_delta_vertex_handle_set_.end(),
                                compare_by_vertex_id));
          std::sort(ctx.partition_delta_vertex_handle_set_.begin(),
                    ctx.partition_delta_vertex_handle_set_.end());
          std::sort(ctx.delta_vertex_handle_set_.begin(),
                    ctx.delta_vertex_handle_set_.end());
          util::Info("without workload balance, worker: "
                  +  std::to_string(ctx.fid_) 
                  + " delta pivots vertexes set size: "
                  +  std::to_string(ctx.partition_delta_vertex_handle_set_.size()));
        }
      }
      else {
        ctx.sorted_vertex_handle_set_.swap(data_graph_vertex_handle_set);
        util::Info("with workload balance, worker: "
                 +  std::to_string(ctx.fid_) 
                 + " holds all pivots vertexes with size: "
                 +  std::to_string(ctx.sorted_vertex_handle_set_.size()));
      }

      #ifndef NDEBUG
      if (!ctx.is_incremental_) {
        if (!ctx.work_load_balance_) {
          assert(ctx.data_graph_pivot_range_.size() 
              >= ctx.data_graph_.CountVertex() / ctx.frag_num_);
          assert(ctx.data_graph_pivot_range_.size() 
              <= ctx.data_graph_.CountVertex() / ctx.frag_num_ + 1);
        }
      }
      else {
        assert(ctx.data_graph_pivot_range_ .size()   
           >= (ctx.data_graph_.CountVertex() - ctx.delta_vertex_handle_set_.size()) / ctx.frag_num_);
        assert(ctx.data_graph_pivot_range_ .size() 
           <= (ctx.data_graph_.CountVertex() - ctx.delta_vertex_handle_set_.size()) / ctx.frag_num_ + 1);
      }
      #endif // NDEBUG
    }

    std::string msg(kInfoProcessPrefix);
    msg += " " + std::to_string(ctx.fid_)
         + " " + std::to_string(omp_get_num_procs());
    auto& channel_0 = messages.Channels()[0];
    channel_0.SendToFragment(kTimmerFragID, msg);
    
    if (ctx.fid_ == kTimmerFragID
     && ctx.time_log_file_.is_open()) {
      ctx.time_log_file_ << timer_now() << std::endl;
    }

    #ifdef PROFILING
    ctx.exec_time -= GetCurrentTime();
    #endif

    #ifdef PROFILING
    ctx.exec_time += GetCurrentTime();
    ctx.postprocess_time -= GetCurrentTime();
    #endif
    // messages.ForceContinue();

    #ifdef PROFILING
    ctx.postprocess_time += GetCurrentTime();
    #endif
  }

  /**
   * @brief Incremental evaluation for GorMatch.
   *
   * @param frag
   * @param ctx
   * @param messageskAllocPivotSetPrefix
   */
  void IncEval(const fragment_t& frag, 
                      context_t& ctx,
              message_manager_t& messages) {

    bool receive_message      = false,
         receive_info_process = false,
         receive_info_finish  = false,
        // ###########################################
        // ##   protocol without workload balance   ##
        // ###########################################
         receive_begin_match  = false,
        // ########################################
        // ##   protocol with workload balance   ##
        // ########################################
         receive_alloc_pivot_set = false,
         receive_info_status     = false,
        // #################################
        // ##   protocol with multi gor   ##
        // #################################
         receive_reset_parameter = false;

    size_t received_pivot_vertex_set_idx_set_begin = 0,
           received_pivot_vertex_set_idx_set_end   = 0;

    messages.ParallelProcess<std::string>(
        // thread_num(),
        1, [&ctx, 
            &received_pivot_vertex_set_idx_set_begin, 
            &received_pivot_vertex_set_idx_set_end,
            &receive_message,
            &receive_info_process, 
            &receive_info_finish,
          // ###########################################
          // ##   protocol without workload balance   ##
          // ###########################################
            &receive_begin_match,
          // ########################################
          // ##   protocol with workload balance   ##
          // ########################################
            &receive_alloc_pivot_set,
            &receive_info_status,
          // #################################
          // ##   protocol with multi gor   ##
          // #################################
            &receive_reset_parameter](int tid, std::string msg) {
          util::Info("fid: "  + std::to_string(ctx.fid_) 
                   + " receive message: " + msg);
          auto res_info_process 
            = std::mismatch(kInfoProcessPrefix.begin(),
                            kInfoProcessPrefix. end (), msg.begin());
          auto res_info_finish 
            = std::mismatch(kInfoFinishPrefix.begin(),
                            kInfoFinishPrefix. end (), msg.begin());
          // ###########################################
          // ##   protocol without workload balance   ##
          // ###########################################
          auto res_begin_match 
            = std::mismatch(kBeginMatchPrefix.begin(),
                            kBeginMatchPrefix. end (), msg.begin());
          // ########################################
          // ##   protocol with workload balance   ##
          // ########################################
          auto res_alloc_pivot_set
            = std::mismatch(kAllocPivotSetPrefix.begin(),
                            kAllocPivotSetPrefix. end (), msg.begin());
          auto res_info_status
            = std::mismatch(kInfoStatusPrefix.begin(),
                            kInfoStatusPrefix. end (), msg.begin());
          // #################################
          // ##   protocol with multi gor   ##
          // #################################
          auto res_reset_parameter
            = std::mismatch(kResetParameterPrefix.begin(),
                            kResetParameterPrefix. end (), msg.begin());

          receive_message = true;
          if (res_info_finish.first == kInfoFinishPrefix.end()) {
            msg = msg.substr(kInfoFinishPrefix.size());
            // kInfoFinishPrefix is the prefix of msg.
            receive_info_finish = true;

            std::string match_count_str;

            std::stringstream ss;
            ss << msg;
            ss >> match_count_str;
            
            ctx.final_result_ += std::stoi(match_count_str);

            util::Info("## receive_info_finish ##");
          } else if (res_info_process.first == kInfoProcessPrefix.end()) {
            msg = msg.substr(kInfoProcessPrefix.size());
            // kInfoProcessPrefix is the prefix of msg.
            receive_info_process = true;
            util::Info("##  receive_info_process  ##");
          } // ###########################################
            // ##   protocol without workload balance   ##
            // ###########################################
            else if (res_begin_match.first == kBeginMatchPrefix.end()) {
            msg = msg.substr(kBeginMatchPrefix.size());
            // kBeginMatchPrefix is the prefix of msg.
            receive_begin_match = true;
            util::Info("##  receive_begin_match  ##");
          } // ########################################
            // ##   protocol with workload balance   ##
            // ########################################
            else if (res_alloc_pivot_set.first == kAllocPivotSetPrefix.end()) {
            msg = msg.substr(kAllocPivotSetPrefix.size());
            // kAllocPivotSetPrefix is the prefix of msg.
            receive_alloc_pivot_set = true;
            std::stringstream ss;
            std::string begin_str, end_str;
            ss << msg;
            ss >> begin_str;
            ss >>   end_str;
            received_pivot_vertex_set_idx_set_begin = std::stoi(begin_str);
            received_pivot_vertex_set_idx_set_end   = std::stoi(  end_str);
            util::Info("##  receive_alloc_pivot_set  ##");
          } else if (res_info_status.first == kInfoStatusPrefix.end()) {
            msg = msg.substr(kInfoStatusPrefix.size());
            // kInfoStatusPrefix is the prefix of msg.
            receive_info_status = true;

            std::string from_id_str,
                         result_str,
                       remained_str;

            std::stringstream ss;
            ss << msg;
            ss >> from_id_str;
            ss >> result_str;
            ss >> remained_str;

            ctx.final_result_ += std::stoi(result_str);
            ctx.remained_pivot_size_[std::stoi(from_id_str)] 
                                   = std::stoi(remained_str);

            util::Info("##  receive_info_status  ##");
          } else if (res_reset_parameter.first == kResetParameterPrefix.end()) {
            msg = msg.substr(kResetParameterPrefix.size());
            // kResetParameterPrefix is the prefix of msg.
            receive_reset_parameter = true;
            util::Info("##  receive_reset_parameter  ##");
          } else {
            // unknown message type
            assert(false);
          }
        });

    if (!receive_message) {
      return;
    }

    if (receive_reset_parameter) {
      ctx.ResetParameter();
      ctx.ToNextGor();
      if (!ctx.HasGorToProcess()) {
        if (ctx.fid_ == kTimmerFragID
         && ctx.GorSetSize() > 1) {
          // multi gar
          if (ctx.time_log_file_.is_open()) {
            ctx.time_log_file_ << " - " << "total match time" 
                               << ": " << ctx.TotalMatchTime() << " sec" << std::endl;
          }
        }
        return;
      }
      // to begin the next round match
      std::string msg(kInfoProcessPrefix);
      msg += " " + std::to_string(ctx.fid_)
           + " " + std::to_string(omp_get_num_procs());
      auto& channel_0 = messages.Channels()[0];
      channel_0.SendToFragment(kTimmerFragID, msg);
      return;
    }

    if (receive_info_process) {
      // each fragment are loaded, timmer node begin to record
      assert (ctx.fid_ == kTimmerFragID);
      timer_next("distributed gor match");
      if (!ctx.work_load_balance_) {
        // without work load balance, just send message
        // to each fragment to info them begin matching
        auto& channel_0 = messages.Channels()[0];
        for (int dst_fid = 0; dst_fid < ctx.frag_num_; dst_fid++) {
          std::string msg(kBeginMatchPrefix);
          channel_0.SendToFragment(dst_fid, std::move(msg));
        }
        return;
      }
      // with work load balance, need to alloc pivot set to each
      // worker to each fragment to info them begin matching
    }

    if (receive_info_finish) {
      // match end
      assert(ctx.fid_ == kTimmerFragID);
      // only timmer node record the time log
      ctx.AddMatchTime(timer_now_time());
      if (ctx.time_log_file_.is_open()) {
        ctx.time_log_file_ << ctx.CurrentGorName() 
                           << " match count: " 
                           << ctx.final_result_ << std::endl;
        ctx.time_log_file_ << timer_now() << std::endl;
      }
      // distributed match end
      util::Output("match count: " + std::to_string(ctx.final_result_) );
      // has gor to process, info all worker to reset the
      // parameter, begin next round
      auto& channel_0 = messages.Channels()[0];
      for (int dst_fid = 0; dst_fid < ctx.frag_num_; dst_fid++) {
        std::string msg(kResetParameterPrefix);
        channel_0.SendToFragment(dst_fid, std::move(msg));
      }
      return;
    }

    if (receive_info_process || receive_info_status) {
      assert(ctx.fid_ == kTimmerFragID);  
      assert(ctx.work_load_balance_);
      // with work load balance, send work load 
      // to each fragment for their processing
      std::vector<std::pair<size_t,  // pivot_vertex_set_idx_set_begin
                            size_t>  // pivot_vertex_set_idx_set_end
                 > pivot_set_for_each_worker(ctx.frag_num_, 
                                             std::pair(0, 0));

      bool all_empty = true;
      for (auto fid = 0; fid < ctx.frag_num_; fid++) {
        util::Debug("fid: " + std::to_string(fid) + " ctx.remained_pivot_size_[fid]: "
                                     + std::to_string(ctx.remained_pivot_size_[fid]));
        
        if (ctx.current_allocated_pivot_set_idx_ 
         >= ctx.sorted_vertex_handle_set_.size()) {
          assert(ctx.current_allocated_pivot_set_idx_ 
              == ctx.sorted_vertex_handle_set_.size());
          // all pivots have been allocated
          if (ctx.remained_pivot_size_[fid] != 0) {
            // still have pivot to process
            all_empty = false;
          }
          continue;
        }

        assert(pivot_set_for_each_worker[fid].first  == 0
            && pivot_set_for_each_worker[fid].second == 0);

        assert(ctx.current_allocated_pivot_set_idx_
             < ctx.sorted_vertex_handle_set_.size());

        pivot_set_for_each_worker[fid].first = ctx.current_allocated_pivot_set_idx_;

        const auto kWorkloadRequires = ctx.pivoted_candidate_balance_size_ 
                                     - ctx.remained_pivot_size_[fid];
        assert(kWorkloadRequires >= 0
            && kWorkloadRequires <= ctx.pivoted_candidate_balance_size_);

        const auto kWorkloadRemained = ctx.sorted_vertex_handle_set_.size()
                                     - ctx.current_allocated_pivot_set_idx_;
        assert(kWorkloadRemained > 0
            && kWorkloadRemained <= ctx.sorted_vertex_handle_set_.size());

        const auto kWorkloadToSend 
                 = kWorkloadRequires < kWorkloadRemained?
                   kWorkloadRequires : kWorkloadRemained;

        pivot_set_for_each_worker[fid].second = pivot_set_for_each_worker[fid].first
                                              + kWorkloadToSend;

        assert(pivot_set_for_each_worker[fid].second > 0
            &&(pivot_set_for_each_worker[fid].second 
             > pivot_set_for_each_worker[fid].first)
            &&(pivot_set_for_each_worker[fid].second <= ctx.sorted_vertex_handle_set_.size()));

        ctx.remained_pivot_size_[fid] += kWorkloadToSend;
        ctx.current_allocated_pivot_set_idx_ += kWorkloadToSend;

        all_empty = false;
      }

      if (all_empty) {
        // all are processed, end the matching
        assert(ctx.fid_ == kTimmerFragID);
        // only timmer node record the time log
        ctx.AddMatchTime(timer_now_time());
        if (ctx.time_log_file_.is_open()) {
          ctx.time_log_file_ << ctx.CurrentGorName() 
                             << " match count: " 
                             << ctx.final_result_ << std::endl;
          ctx.time_log_file_ << timer_now() << std::endl;
        }
        // distributed match end
        util::Output("match count: " + std::to_string(ctx.final_result_) );
        // has gor to process, begin next round
        auto& channel_0 = messages.Channels()[0];
        for (int dst_fid = 0; dst_fid < ctx.frag_num_; dst_fid++) {
          std::string msg(kResetParameterPrefix);
          channel_0.SendToFragment(dst_fid, std::move(msg));
        }
        return;
      }
      auto& channel_0 = messages.Channels()[0];
      for (int dst_fid = 0; dst_fid < ctx.frag_num_; dst_fid++) {
        std::string msg(kAllocPivotSetPrefix);
        msg = std::move(msg) + " " + std::to_string(pivot_set_for_each_worker[dst_fid].first)
                             + " " + std::to_string(pivot_set_for_each_worker[dst_fid].second);
        channel_0.SendToFragment(dst_fid, std::move(msg));
      }
      return;
    }

    assert(( ctx.work_load_balance_ && receive_alloc_pivot_set)
        || (!ctx.work_load_balance_ && receive_begin_match));
    if (!ctx.candidate_set_initialized_) {
      ctx.candidate_set_.clear();
      if (!GUNDAM::_dp_iso_using_match::InitCandidateSet<GUNDAM::MatchSemantics::kIsomorphism>(
              ctx.CurrentGor().pattern(), 
              ctx.data_graph_, 
              ctx.candidate_set_)) {
        util::Debug("InitCandidateSet fail!");
        std::string msg(kInfoFinishPrefix);
        auto& channel_0 = messages.Channels()[0];
        channel_0.SendToFragment(kTimmerFragID, msg + " 0 0"); // supp_x = 0 && supp_y = 0
        return;
      }
      //  if need to workload balance:
      //    refine the entire candidate set at once
      //  if do not need to workload balance:
      //    receive the partition, then refine the candidate set
      if (ctx.work_load_balance_) {
        // with workload balance, refine the entire candidate set at once
        if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(
                  ctx.CurrentGor().pattern(), 
                  ctx.data_graph_, 
                  ctx.candidate_set_)) {
          util::Debug("RefineCandidateSet fail!");
          std::string msg(kInfoFinishPrefix);
          auto& channel_0 = messages.Channels()[0];
          channel_0.SendToFragment(kTimmerFragID, msg + " 0 0"); // supp_x = 0 && supp_y = 0
          return;
        }
      }
      ctx.candidate_set_initialized_ = true;
    }

    if (receive_alloc_pivot_set) {
      // #############################
      // ##  receive the pivot set  ##
      // #############################
      assert( ctx.work_load_balance_);
      assert(!ctx.is_incremental_); // incremental does not support workload balance
      // add allocated pivot set to this worker

      #ifndef NDEBUG
      const auto kSizeBeforeInsert = ctx.data_graph_pivot_range_.size();
      #endif

      ctx.data_graph_pivot_range_.insert(
      ctx.data_graph_pivot_range_  .end(),
      ctx.sorted_vertex_handle_set_.begin() + received_pivot_vertex_set_idx_set_begin,
      ctx.sorted_vertex_handle_set_.begin() + received_pivot_vertex_set_idx_set_end);

      #ifndef NDEBUG
      if (received_pivot_vertex_set_idx_set_begin 
       == received_pivot_vertex_set_idx_set_end) {
        assert(kSizeBeforeInsert == ctx.data_graph_pivot_range_.size());
      }
      #endif

      // receive all pivot vertexes, needs to sort it 
      std::sort(ctx.data_graph_pivot_range_.begin(),
                ctx.data_graph_pivot_range_.end());
    }

    assert(receive_alloc_pivot_set
        || receive_begin_match);

    util::Debug("wenzhi 1 here!");

    auto pivot_candidate_set_it 
         = ctx.candidate_set_.find(ctx.CurrentPivot());
    assert(pivot_candidate_set_it 
          != ctx.candidate_set_.end());
    
    std::vector<DataGraphVertexHandle> selected_pivot_candidate_set;

    util::Debug("pivot_candidate_set_it->second.size(): "
                + std::to_string(pivot_candidate_set_it->second.size()));
    util::Debug("ctx.data_graph_pivot_range_.size(): "
                + std::to_string(ctx.data_graph_pivot_range_.size()));

    std::set_intersection(pivot_candidate_set_it->second.begin(),
                          pivot_candidate_set_it->second. end (),
                          ctx.data_graph_pivot_range_.begin(),
                          ctx.data_graph_pivot_range_. end (),
                          std::back_inserter(
                                selected_pivot_candidate_set));

    util::Debug("selected_pivot_candidate_set.size(): "
                + std::to_string(selected_pivot_candidate_set.size()));

    #ifndef NDEBUG
    const size_t kPivotCandidateSetSizeBeforeSwap
                = pivot_candidate_set_it->second.size();
    #endif // NDEBUG

    if (!ctx.work_load_balance_) {

      if (ctx.is_incremental_) {
        // ################################################
        // ##  first find how many existed pivot vertex  ##
        // ##  can be updated by the delta vertex set    ##
        // ################################################
        const auto kPatternRadius = GUNDAM::Radius(ctx.CurrentGor().pattern(),
                                                   ctx.CurrentPivot());

        std::set<VertexLabelType> vertex_label_set;
        for (auto vertex_it = ctx.CurrentGor().pattern().VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          vertex_label_set.emplace(vertex_it->label());
        } 

        std::  set <DataGraphVertexHandle> src_vertex_handle_set;
        std::vector<DataGraphVertexHandle> src_vertex_handle_set_in_vector;
        assert(std::is_sorted(ctx.delta_vertex_handle_set_.begin(),
                              ctx.delta_vertex_handle_set_.end()));
        for (const auto& delta_vertex_handle 
                   : ctx.delta_vertex_handle_set_) {
          if (vertex_label_set.find(delta_vertex_handle->label())
           == vertex_label_set.end()) {
            continue;
          }
          src_vertex_handle_set.emplace(delta_vertex_handle);
          src_vertex_handle_set_in_vector.emplace_back(delta_vertex_handle);
        }
        std::sort(src_vertex_handle_set_in_vector.begin(),
                  src_vertex_handle_set_in_vector.end());

        std::vector<DataGraphVertexHandle> vertex_handle_set_in_range;
        auto find_vertex_handle_callback 
              = [&vertex_handle_set_in_range](
              const DataGraphVertexHandle& vertex_handle) {
          vertex_handle_set_in_range.emplace_back(vertex_handle);
          return true;
        };

        auto bfs_ret = GUNDAM::Bfs(ctx.data_graph_,
                                   src_vertex_handle_set,
                                  find_vertex_handle_callback,
                                  kPatternRadius);

        assert(vertex_handle_set_in_range.size()
            == bfs_ret);

        std::sort(vertex_handle_set_in_range.begin(),
                  vertex_handle_set_in_range.end());

        std::vector<DataGraphVertexHandle> vertex_handle_set_involved;
        for (auto vertex_it = ctx.CurrentGor().pattern().VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          CandidateSet temp_candidate_set = ctx.candidate_set_;
          assert(temp_candidate_set.find(vertex_it)
              != temp_candidate_set.end());
          temp_candidate_set[vertex_it] = src_vertex_handle_set_in_vector;
          
          if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(
                  ctx.CurrentGor().pattern(), 
                  ctx.data_graph_, 
                temp_candidate_set)) {
            continue;
          }
          assert(temp_candidate_set.find(ctx.CurrentPivot())
              != temp_candidate_set.end());
            
          std::vector<DataGraphVertexHandle> temp_vertex_handle_set_involved;

          std::set_union(vertex_handle_set_involved.begin(), 
                         vertex_handle_set_involved.end(),
                        temp_candidate_set[ctx.CurrentPivot()].begin(), 
                        temp_candidate_set[ctx.CurrentPivot()].end(),
                        std::back_inserter(temp_vertex_handle_set_involved));

          assert(temp_vertex_handle_set_involved.size()
                   >= vertex_handle_set_involved.size());

          vertex_handle_set_involved.swap(temp_vertex_handle_set_involved);
        }

        std::vector<DataGraphVertexHandle> temp_selected_pivot_candidate_set;

        // add all selected_pivot_candidate_set in it at once
        // for candidate set contained in this partition
        pivot_candidate_set_it->second.swap(selected_pivot_candidate_set);

        // vertex contained in the radius
        std::set_intersection(pivot_candidate_set_it->second.begin(),
                              pivot_candidate_set_it->second. end (),
                              vertex_handle_set_involved.begin(),
                              vertex_handle_set_involved. end (),
                              std::back_inserter(
                               temp_selected_pivot_candidate_set));

        temp_selected_pivot_candidate_set.swap(vertex_handle_set_involved);
        temp_selected_pivot_candidate_set.clear();

        std::set_intersection(vertex_handle_set_in_range.begin(),
                              vertex_handle_set_in_range. end (),
                              vertex_handle_set_involved.begin(),
                              vertex_handle_set_involved. end (),
                              std::back_inserter(
                               temp_selected_pivot_candidate_set));

        pivot_candidate_set_it->second.swap(temp_selected_pivot_candidate_set);

        util::Debug("before partition vertex GORMatch pivot candidate set: " 
                                     + std::to_string(pivot_candidate_set_it->second.size()));

        CandidateSet temp_candidate_set = ctx.candidate_set_;

        size_t result = 0;

        if (GUNDAM::_dp_iso_using_match::RefineCandidateSet(
                ctx.CurrentGor().pattern(), 
                ctx.data_graph_, 
               temp_candidate_set)) {
          std::vector<CandidateSet> match_result;
          util::Debug("before partition vertex GORMatch");
          result += gor::GORMatch(
                        ctx.CurrentGor(), 
                        ctx.CurrentPivot(), 
                        ctx.data_graph_,
                       temp_candidate_set,
                        match_result,
                        ctx.time_limit_);
        }
        #ifndef NDEBUG 
        else {
          util::Debug("RefineCandidateSet fail of match in partition vertex set in worker: "
                    + std::to_string(ctx.fid_));
        }
        #endif // NDEBUG 

        util::Debug("result after partition vertex GORMatch: " + std::to_string(result));

        // #################################################
        // ##  then to find the pivot vertex contained    ##
        // ##  in the partition_delta_vertex_handle_set_  ##
        // #################################################
        
        // swap back
        pivot_candidate_set_it->second.swap(selected_pivot_candidate_set);
        assert(pivot_candidate_set_it->second.size()
           == kPivotCandidateSetSizeBeforeSwap);

        selected_pivot_candidate_set.clear();

        std::set_intersection(pivot_candidate_set_it->second.begin(),
                              pivot_candidate_set_it->second. end (),
                              ctx.partition_delta_vertex_handle_set_.begin(),
                              ctx.partition_delta_vertex_handle_set_. end (),
                              std::back_inserter(
                                    selected_pivot_candidate_set));

        util::Debug("selected_pivot_candidate_set.size(): "
                    + std::to_string(selected_pivot_candidate_set.size()));

        pivot_candidate_set_it->second.swap(selected_pivot_candidate_set);

        util::Debug("before delta vertex GORMatch pivot candidate set: " 
                                     + std::to_string(pivot_candidate_set_it->second.size()));

        if (GUNDAM::_dp_iso_using_match::RefineCandidateSet(
                ctx.CurrentGor().pattern(), 
                ctx.data_graph_, 
                ctx.candidate_set_)) {
          std::vector<CandidateSet> match_result;
          util::Debug("before delta vertex GORMatch");
          result += gor::GORMatch(
                        ctx.CurrentGor(), 
                        ctx.CurrentPivot(), 
                        ctx.data_graph_,
                        ctx.candidate_set_,
                        match_result,
                        ctx.time_limit_);
        }
        #ifndef NDEBUG 
        else {
          util::Debug("RefineCandidateSet fail of match in delta vertex set in worker: "
                    + std::to_string(ctx.fid_));
        }
        #endif // NDEBUG 

        util::Debug("result after delta vertex GORMatch: " + std::to_string(result));

        std::string msg(kInfoFinishPrefix);
        auto& channel_0 = messages.Channels()[0];
        channel_0.SendToFragment(kTimmerFragID, msg + " " + std::to_string(result));
        return;
      }

      // is not incremental

      // add all selected_pivot_candidate_set in it at once
      // for candidate set contained in this partition
      pivot_candidate_set_it->second.swap(selected_pivot_candidate_set);

      size_t result = 0;
      if (GUNDAM::_dp_iso_using_match::RefineCandidateSet(
              ctx.CurrentGor().pattern(), 
              ctx.data_graph_, 
              ctx.candidate_set_)) {
        std::vector<CandidateSet> match_result;
        result = gor::GORMatch(
                      ctx.CurrentGor(), 
                      ctx.CurrentPivot(), 
                      ctx.data_graph_,
                      ctx.candidate_set_,
                      match_result,
                      ctx.time_limit_);
      }
      else {
        util::Debug("RefineCandidateSet fail!");
      }

      // unnecessary to swap back, the kTimmerFragID
      // would info each worker to reset the parameter
      std::string msg(kInfoFinishPrefix);
      auto& channel_0 = messages.Channels()[0];
      channel_0.SendToFragment(kTimmerFragID, msg + " " + std::to_string(result));
      return;
    }

    assert(!ctx.is_incremental_);
    assert(ctx.work_load_balance_);
    // with workload balance
    // add all selected_pivot_candidate_set one-by-one

    auto kBeginTime = std::chrono::high_resolution_clock::now();

    ctx.data_graph_pivot_range_.swap(selected_pivot_candidate_set);

    std::vector<DataGraphVertexHandle>  temp_pivot_candidate_set;
    pivot_candidate_set_it->second.swap(temp_pivot_candidate_set);
    pivot_candidate_set_it->second.clear();
    pivot_candidate_set_it->second.resize(1);
    assert(pivot_candidate_set_it->second.size() == 1);

    size_t result = 0;

    while (!ctx.data_graph_pivot_range_.empty()) {

      assert(pivot_candidate_set_it->second.size() == 1);

      pivot_candidate_set_it->second[0] = ctx.data_graph_pivot_range_.back();
      
      std::vector<CandidateSet> match_result;
      
      result += gor::GORMatch(
                    ctx.CurrentGor(), 
                    ctx.CurrentPivot(), 
                    ctx.data_graph_,
                    ctx.candidate_set_,
                    match_result,
                    ctx.time_limit_);

      // this pivot is processed
      ctx.data_graph_pivot_range_.pop_back();

      if ((std::chrono::duration_cast<std::chrono::milliseconds>(
           std::chrono::high_resolution_clock::now() - kBeginTime).count()
              > ctx.balance_period_ms_)) {
        break;
      }
    }

    // store back
    pivot_candidate_set_it->second.swap(temp_pivot_candidate_set);
    
    std::string msg(kInfoStatusPrefix);
    auto& channel_0 = messages.Channels()[0];
    channel_0.SendToFragment(kTimmerFragID, msg + " " + std::to_string(ctx.fid_) 
                                                + " " + std::to_string(result) 
                                                + " " + std::to_string(ctx.data_graph_pivot_range_.size()));
  
#ifdef PROFILING
    ctx.preprocess_time -= GetCurrentTime();
#endif

#ifdef PROFILING
    ctx.preprocess_time += GetCurrentTime();
    ctx.exec_time -= GetCurrentTime();
#endif

#ifdef PROFILING
    ctx.exec_time += GetCurrentTime();
    ctx.postprocess_time -= GetCurrentTime();
#endif

#ifdef PROFILING
    ctx.postprocess_time += GetCurrentTime();
#endif

    return;
  }
};

}  // namespace grape

#endif  // EXAMPLES_ANALYTICAL_APPS_GOR_MATCH_GOR_MATCH_H_