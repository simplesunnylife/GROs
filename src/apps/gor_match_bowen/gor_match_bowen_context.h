#ifndef EXAMPLES_ANALYTICAL_APPS_GOR_MATCH_GOR_MATCH_CONTEXT_H_
#define EXAMPLES_ANALYTICAL_APPS_GOR_MATCH_GOR_MATCH_CONTEXT_H_

#include <grape/app/context_base.h>
#include <grape/grape.h>

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <queue>

#include "gor/gor.h"

#include "gundam/graph_type/small_graph.h"
#include "gundam/graph_type/large_graph.h"
#include "gundam/graph_type/graph.h"

namespace grape {

/**
 * @brief Context for the parallel version of GorMatch.
 *
 * @tparam FRAG_T
 */
template <typename FRAG_T>
class GorMatchContext : public VertexDataContext<FRAG_T, int64_t> {
 public:
  using oid_t = typename FRAG_T::oid_t;
  using vid_t = typename FRAG_T::vid_t;

  explicit GorMatchContext(const FRAG_T& fragment)
       : VertexDataContext<FRAG_T, int64_t>(fragment, true){
    return;
  }

  void Init(ParallelMessageManager& messages, std::string yaml_file,
            int frag_num) {
    this->yaml_file_ = yaml_file;
    this->fid_ = this->fragment().fid();
    this->frag_num_ = frag_num;
    // #######################
    // #  program parameter  #
    // #######################
    this->ResetParameter();
    this->current_processing_gor_idx_ = 0;
    this->total_match_time_ = 0;
#ifdef PROFILING
    preprocess_time = 0;
    exec_time = 0;
    postprocess_time = 0;
#endif
  }

  void Output(std::ostream& os) {
#ifdef PROFILING
    VLOG(2) << "preprocess_time: " << preprocess_time << "s.";
    VLOG(2) << "exec_time: " << exec_time << "s.";
    VLOG(2) << "postprocess_time: " << postprocess_time << "s.";
#endif
  }

#ifdef PROFILING
  double preprocess_time = 0;
  double exec_time = 0;
  double postprocess_time = 0;
#endif

  using VertexIDType = uint32_t;
  using   EdgeIDType = uint32_t;

  using VertexLabelType = uint32_t;
  using   EdgeLabelType = uint32_t;

  // using Pattern   = GUNDAM::LargeGraph2<VertexIDType, VertexLabelType, std::string,
  //                                         EdgeIDType,   EdgeLabelType, std::string>;
  
  using   Pattern = GUNDAM::SmallGraph<VertexIDType, VertexLabelType, 
                                         EdgeIDType,   EdgeLabelType>;

  // using DataGraph = GUNDAM::LargeGraph<VertexIDType, VertexLabelType, std::string,
  //                                        EdgeIDType,   EdgeLabelType, std::string>;

  // using Pattern = GUNDAM::SmallGraph<VertexIDType, VertexLabelType, 
  //                                      EdgeIDType,   EdgeLabelType>;

  using DataGraph = GUNDAM::Graph<
      GUNDAM::SetVertexIDType<VertexIDType>,
      GUNDAM::SetVertexAttributeStoreType<GUNDAM::AttributeType::kGrouped>,
      // GUNDAM::SetVertexAttributeStoreType<GUNDAM::AttributeType::kSeparated>,
      GUNDAM::SetVertexLabelType<VertexLabelType>,
      GUNDAM::SetVertexLabelContainerType<GUNDAM::ContainerType::Map>,
      GUNDAM::SetVertexIDContainerType<GUNDAM::ContainerType::Map>,
      GUNDAM::SetVertexPtrContainerType<GUNDAM::ContainerType::Vector>,
      GUNDAM::SetEdgeLabelContainerType<GUNDAM::ContainerType::Vector>,
      GUNDAM::SetVertexAttributeKeyType<std::string>,
      GUNDAM::SetEdgeIDType<EdgeIDType>,
      GUNDAM::SetEdgeAttributeStoreType<GUNDAM::AttributeType::kGrouped>,
      // GUNDAM::SetEdgeAttributeStoreType<GUNDAM::AttributeType::kSeparated>,
      GUNDAM::SetEdgeLabelType<EdgeLabelType>,
      GUNDAM::SetEdgeAttributeKeyType<std::string>>;

  using GorType = gor::GraphOracleRule<Pattern, DataGraph>;

  using CandidateSet = std::map<typename GUNDAM::VertexHandle< Pattern >::type, 
                    std::vector<typename GUNDAM::VertexHandle<DataGraph>::type>>;

  inline void ToNextGor() {
    assert(this->current_processing_gor_idx_ 
         < this->gor_set_.size());
    this->current_processing_gor_idx_++;
    return;
  }

  inline bool HasGorToProcess() const {
    return this->current_processing_gor_idx_ 
         < this->gor_set_.size();
  }

  inline bool HasProcessedAllGor() const {
    return !this->HasGorToProcess();
  }

  inline GorType& CurrentGor() {
    assert(this->current_processing_gor_idx_ >= 0 
        && this->current_processing_gor_idx_ < this->gor_set_.size());
    return this->gor_set_[this->current_processing_gor_idx_].first;
  }

  inline std::string CurrentGorName() const {
    return "# gor # vfile: " + std::get<0>(this->gor_path_set_[this->current_processing_gor_idx_])
              + " # efile: " + std::get<1>(this->gor_path_set_[this->current_processing_gor_idx_])
              + " # xfile: " + std::get<2>(this->gor_path_set_[this->current_processing_gor_idx_])
              + " # yfile: " + std::get<3>(this->gor_path_set_[this->current_processing_gor_idx_])
              + " # pivot: " + GUNDAM::ToString(this->CurrentPivot()->id());
              + " #";
  }

  inline const typename GUNDAM::VertexHandle<Pattern>::type& CurrentPivot() const {
    assert(this->current_processing_gor_idx_ >= 0 
        && this->current_processing_gor_idx_ < this->gor_set_.size());
    assert(this->gor_set_[this->current_processing_gor_idx_].second);
    assert(this->gor_set_[this->current_processing_gor_idx_].second
        == this->gor_set_[this->current_processing_gor_idx_]
                .first.pattern().FindVertex(this->gor_set_[0].second->id()));
    return this->gor_set_[this->current_processing_gor_idx_].second;
  }

  inline const size_t GorSetSize() const {
    return this->gor_set_.size();
  }

  inline void ResetParameter() {
    this->candidate_set_initialized_ = false;
    this->remained_pivot_size_.clear();
    this->remained_pivot_size_.resize(this->frag_num_, 0);
    this->current_allocated_pivot_set_idx_ = 0;
    this->final_result_ = 0;
    this->candidate_set_.clear();
    if (this->work_load_balance_) {
      // in this case, the data_graph_pivot_range_ is allocated
      // from the coordinator, needs to be reset in each round for
      // different gars
      assert(this->data_graph_pivot_range_.empty());
      this->data_graph_pivot_range_.clear();
    }
    return;
  }

  inline void AddMatchTime(double time) {
    this->total_match_time_ += time;
    return;
  }

  inline double TotalMatchTime() const {
    return this->total_match_time_;
  }

  // inline auto PartitionDeltaVertexBegin() {
  //   return ;
  // }

  // inline auto PartitionDeltaVertexEnd() {
  //   return ;
  // }

  // inline auto PartitionVertexBegin() {
  //   return this->;
  // }

  // inline auto PartitionVertexEnd() {
  //   return ;
  // }

  std::string yaml_file_;
  int fid_;
  int frag_num_;

  DataGraph data_graph_;

  // partition the pivot range in data graph
  // for the incremental, partition the original graph 
  std::vector<typename GUNDAM::VertexHandle<DataGraph>::type> 
    data_graph_pivot_range_;

  CandidateSet candidate_set_;

  // also need to be stored, to output the name of current processing gor
  std::vector<std::tuple<std::string,
                         std::string,
                         std::string,
                         std::string,
                         std::string>> gor_path_set_;

  std::vector<std::pair<GorType, typename GUNDAM::VertexHandle<Pattern>::type>> 
    gor_set_;

  bool candidate_set_initialized_;

  size_t final_result_;

  std::ofstream time_log_file_;

  double time_limit_;
  double time_limit_per_supp_;

  // ############################################
  // ##   parameter for incremental matching   ##
  // ############################################
  bool is_incremental_;
  std::vector<typename GUNDAM::VertexHandle<DataGraph>::type> 
  // holds all delta vertexes, to update the match in existed pivot set
              delta_vertex_handle_set_,
  // holds delta vertexes in this partition, to update the match in the delta vertex set
    partition_delta_vertex_handle_set_;

  // ########################################
  // ##   parameter for workload balance   ##
  // ########################################
  bool work_load_balance_;
  int32_t balance_period_ms_;
  int32_t pivoted_candidate_balance_size_;

  // to be stored in each worker, for more efficient transform
  std::vector<typename GUNDAM::VertexHandle<DataGraph>::type> 
    sorted_vertex_handle_set_;

  // maintained in coordinator(kTimmerFragID), to monitor
  // how many pivot are remained in each worker
  std::vector<size_t> remained_pivot_size_;

  // maintained in coordinator(kTimmerFragID), to mark how many
  // pivots are allocated to workers
  size_t current_allocated_pivot_set_idx_;

 private:
  size_t current_processing_gor_idx_;

  // maintained in coordinator(kTimmerFragID), calc the total time 
  // for multi gors
  double total_match_time_;
};

}  // namespace grape

#endif  // EXAMPLES_ANALYTICAL_APPS_GOR_MATCH_GOR_MATCH_CONTEXT_H_
