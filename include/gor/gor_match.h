#ifndef GOR_GOR_MATCH_H_
#define GOR_GOR_MATCH_H_

#include "gor.h"

#include "include/gundam/type_getter/edge_handle.h"
#include "include/gundam/type_getter/vertex_handle.h"

#include "include/gundam/algorithm/simulation.h"

#include "include/util/log.h"

#include "gar/gar_supp.h"

#include <vector>
#include <map>

namespace gor {

template <typename GraphPatternType, 
          typename    DataGraphType>
size_t GORMatch(GraphOracleRule<GraphPatternType, 
                                   DataGraphType>& gor,
  const typename 
    GUNDAM::VertexHandle<GraphPatternType>::type& pivot,
                                   DataGraphType& data_graph,
    const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
       std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& candidate_set_for_pivot_vertex,
    const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
       std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& candidate_set_for_all_vertexes,
       std::vector<
          std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
       std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>>& maximal_match_set,
            double time_limit = -1.0) {
  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using CandidateSet = std::map<GraphPatternVertexHandle,
                       std::vector<DataGraphVertexHandle>>;

  assert(pivot);
  assert(pivot == gor.pattern().FindVertex(pivot->id()));

  auto begin_time = std::time(NULL);

  assert( candidate_set_for_pivot_vertex.find(pivot)
       != candidate_set_for_pivot_vertex.end());
  assert(!candidate_set_for_pivot_vertex.find(pivot)->second.empty() );

  const std::vector<DataGraphVertexHandle>& 
      pivot_candidate = candidate_set_for_pivot_vertex.find(pivot)->second;

  assert(!pivot_candidate.empty());

  GraphPatternType rhs_pattern(
    gar::_gar_supp::LiteralPattern(gor,  true));

  std::vector<GraphPatternVertexHandle> rhs_vertex_handle_set;
  for (auto vertex_it = rhs_pattern.VertexBegin();
           !vertex_it.IsDone();
            vertex_it++) {
    rhs_vertex_handle_set.emplace_back(gor.pattern().FindVertex(vertex_it->id()));
    assert(rhs_vertex_handle_set.back());
  }
  std::sort(rhs_vertex_handle_set.begin(),
            rhs_vertex_handle_set.end());

  size_t count = 0;

  for (const auto& pivot_candidate_handle 
                 : pivot_candidate) {

    if (time_limit > 0 
     && time_limit < (std::time(NULL) - begin_time)) {
      // has reached time limit
      break;
    }
    
    CandidateSet pivoted_candidate_set(
                         candidate_set_for_all_vertexes);
    assert(pivoted_candidate_set.find(pivot)
        != pivoted_candidate_set.end());
    pivoted_candidate_set[pivot].clear();
    assert(pivoted_candidate_set.find(pivot)->second.empty());
    pivoted_candidate_set[pivot].emplace_back(pivot_candidate_handle);

    if (!GUNDAM::_dp_iso_using_match
               ::RefineCandidateSet(gor.pattern(),
                                       data_graph, 
                            pivoted_candidate_set, pivot)) {
      continue;
    }
    assert(pivoted_candidate_set.size() == gor.pattern().CountVertex());

    size_t pivoted_count = 1;
    for (const auto& [query_vertex_handle,
                      dst_candidate_set] : pivoted_candidate_set) {
      if (!std::binary_search(rhs_vertex_handle_set.begin(),
                              rhs_vertex_handle_set.end(),
                            query_vertex_handle)) {
        continue;
      }
      // util::Debug("wenzhi here!");
      pivoted_count *= dst_candidate_set.size();
    }

    maximal_match_set.emplace_back(std::move(pivoted_candidate_set));
    count += pivoted_count;
  }

  return count;
}

template <typename GraphPatternType, 
          typename    DataGraphType>
size_t GORMatch(GraphOracleRule<GraphPatternType, 
                                   DataGraphType>& gor,
  const typename 
    GUNDAM::VertexHandle<GraphPatternType>::type& pivot,
                                   DataGraphType& data_graph,
    const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
       std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& candidate_set_for_pivot_vertex,
       std::vector<
          std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
       std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>>& maximal_match_set,
            double time_limit = -1.0) {
  return GORMatch(gor, pivot, data_graph,
          candidate_set_for_pivot_vertex,
          candidate_set_for_pivot_vertex, maximal_match_set, time_limit);
}

template <typename GraphPatternType, 
          typename    DataGraphType>
size_t GORMatch(GraphOracleRule<GraphPatternType, 
                                   DataGraphType>& gor,
  typename GUNDAM::VertexHandle<GraphPatternType>::type& pivot,
                                   DataGraphType& data_graph,
          std::vector<
          std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
          std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type>>>& maximal_match_set,
            double time_limit = -1.0) {

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using CandidateSet = std::map<GraphPatternVertexHandle,
                       std::vector<DataGraphVertexHandle>>;  

  CandidateSet candidate_set;

  auto res = GUNDAM::Simulation<GUNDAM::MatchSemantics::kDualSimulation>(
                    gor.pattern(), 
                       data_graph,
                    candidate_set);

  if (res == 0){
    return 0;
  }
  
  return GORMatch(gor, pivot, data_graph, candidate_set, maximal_match_set, time_limit);
}

template <typename GraphPatternType, 
          typename    DataGraphType>
size_t GORMatch(GraphOracleRule<GraphPatternType, 
                                   DataGraphType>& gor,
  typename GUNDAM::VertexHandle<GraphPatternType>::type& pivot,
                                   DataGraphType& data_graph,
            double time_limit = -1.0) {

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  std::vector<std::map<GraphPatternVertexHandle,
              std::vector<DataGraphVertexHandle>>> maximal_match_set;

  return GORMatch(gor, pivot, data_graph, maximal_match_set, time_limit);
}

template <typename GraphPatternType, 
          typename    DataGraphType>
size_t IncrementalGORMatch(    GraphOracleRule<GraphPatternType, 
                                                  DataGraphType>& gor,
           const typename GUNDAM::VertexHandle<GraphPatternType>::type& pivot,
                                                  DataGraphType& data_graph,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& candidate_set,
  const 
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>& delta_target_graph,
     std::vector<
        std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>>& maximal_match_set) {

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using CandidateSet = std::map<GraphPatternVertexHandle,
                       std::vector<DataGraphVertexHandle>>;

  assert(pivot && pivot == gor.pattern().FindVertex(pivot->id()));

  // util::Debug("##wenzhi here in GORMatch!##");

  assert( candidate_set.find(pivot)
       != candidate_set.end());
  assert(!candidate_set.find(pivot)->second.empty() );

  const std::vector<DataGraphVertexHandle>& 
      pivot_candidate = candidate_set.find(pivot)->second;

  assert(!pivot_candidate.empty());

  GraphPatternType rhs_pattern(
    gar::_gar_supp::LiteralPattern(gor,  true));

  std::vector<GraphPatternVertexHandle> rhs_vertex_handle_set;
  for (auto vertex_it = rhs_pattern.VertexBegin();
           !vertex_it.IsDone();
            vertex_it++) {
    rhs_vertex_handle_set.emplace_back(gor.pattern().FindVertex(vertex_it->id()));
    assert(rhs_vertex_handle_set.back());
  }
  std::sort(rhs_vertex_handle_set.begin(),
            rhs_vertex_handle_set.end());

  assert(std::is_sorted(delta_target_graph.begin(),
                        delta_target_graph. end ()));

  std::vector<std::tuple<GraphPatternVertexHandle,
                std::vector<DataGraphVertexHandle>, // target_list_in_delta_graph
                std::vector<DataGraphVertexHandle>> // target_list_in_original_graph
             > has_delta_target_graph_pattern_vertex;
  
  for (auto &[query_handle, target_list] : candidate_set) {
    assert(std::is_sorted(target_list.begin(),
                          target_list. end ()));

    if (pivot == query_handle) {
      continue;
    }

    std::vector<DataGraphVertexHandle> target_list_in_delta_graph;
    std::set_intersection(target_list .begin(),
                          target_list . end (),
                    delta_target_graph.begin(),
                    delta_target_graph. end (),
       std::back_inserter(target_list_in_delta_graph));

    if (target_list_in_delta_graph.empty()) {
      // does not have candidate contained in the delta graph
      continue;
    }
    
    std::vector<DataGraphVertexHandle> target_list_in_original_graph;
    std::set_difference(target_list .begin(),
                        target_list . end (),
                  delta_target_graph.begin(),
                  delta_target_graph. end (),
     std::back_inserter(target_list_in_original_graph));

    assert(target_list_in_original_graph.size()
         + target_list_in_delta_graph.size()
        == target_list.size());

    has_delta_target_graph_pattern_vertex.emplace_back();

    std::get<0>(has_delta_target_graph_pattern_vertex.back()) = query_handle;
    std::get<1>(has_delta_target_graph_pattern_vertex.back()).swap(target_list_in_delta_graph);
    std::get<2>(has_delta_target_graph_pattern_vertex.back()).swap(target_list_in_original_graph);
  }
  std::sort(has_delta_target_graph_pattern_vertex.begin(),
            has_delta_target_graph_pattern_vertex.end());
          
  const int kTotalMask = (1 << (has_delta_target_graph_pattern_vertex.size()));
  // pre-process all candidate set
  std::vector<CandidateSet> total_masked_pivoted_candidate_set;
  total_masked_pivoted_candidate_set.reserve(kTotalMask - 1);
  for (int mask = 1; mask < kTotalMask; mask++) {
    CandidateSet pivoted_candidate_set;
    bool illegal_candidate_set = false;
    for (int bit_pos = has_delta_target_graph_pattern_vertex.size() - 1;
             bit_pos >= 0 ; 
             bit_pos--) {
      const auto&        query_vertex_handle    = std::get<0>(has_delta_target_graph_pattern_vertex[bit_pos]);
      const auto& target_list_in_delta_graph    = std::get<1>(has_delta_target_graph_pattern_vertex[bit_pos]);
      const auto& target_list_in_original_graph = std::get<2>(has_delta_target_graph_pattern_vertex[bit_pos]);
      assert(pivoted_candidate_set.find(query_vertex_handle)
          == pivoted_candidate_set.end());
      assert(!target_list_in_delta_graph.empty());
      if (mask & (1 << bit_pos)) {
        // select target_list_in_delta_graph for this vertex
        auto [ pivoted_candidate_set_it,
               pivoted_candidate_set_ret ]
             = pivoted_candidate_set.emplace(query_vertex_handle,
                                             target_list_in_delta_graph);
        // should have added successfully
        assert(pivoted_candidate_set_ret);
        continue;
      }
      // select target_list_in_original_graph for this vertex
      if (target_list_in_original_graph.empty()) {
        illegal_candidate_set = true;
        break;
      }
      auto [ pivoted_candidate_set_it,
             pivoted_candidate_set_ret ]
           = pivoted_candidate_set.emplace(query_vertex_handle,
                                           target_list_in_original_graph);
      // should have added successfully
      assert(pivoted_candidate_set_ret);
    }

    if (illegal_candidate_set) {
      assert(pivoted_candidate_set.size() < has_delta_target_graph_pattern_vertex.size());
      continue;
    }
    assert(pivoted_candidate_set.size() == has_delta_target_graph_pattern_vertex.size());
    for (const auto& [query_vertex,
                      target_set] : candidate_set) {
      if (query_vertex == pivot) {
        continue;
      }
      auto pivoted_candidate_set_it 
         = pivoted_candidate_set.find(query_vertex);
      if (pivoted_candidate_set_it
       != pivoted_candidate_set.end()) {
        // has already added in the pivoted_candidate_set
        continue;
      }
      pivoted_candidate_set.emplace_hint(pivoted_candidate_set_it,
                                         query_vertex, target_set);
    }
    assert(pivoted_candidate_set.size()
        == gor.pattern().CountVertex() - 1);
    total_masked_pivoted_candidate_set.emplace_back(
       std::move(pivoted_candidate_set));
  }

  size_t count = 0;

  for (const auto& pivot_candidate_handle 
                 : pivot_candidate) {  
    std::vector<std::vector<DataGraphVertexHandle>> rhs_vertex_candidate_set;
    rhs_vertex_candidate_set.resize(rhs_vertex_handle_set.size());

    util::Debug("################");
    
    for (const auto& masked_pivoted_candidate_set
             : total_masked_pivoted_candidate_set) {
      CandidateSet pivoted_candidate_set(masked_pivoted_candidate_set);
      assert(pivoted_candidate_set.find(pivot)
          == pivoted_candidate_set.end());
      auto [ pivoted_candidate_set_it,
             pivoted_candidate_set_ret ]
           = pivoted_candidate_set.emplace(pivot, 
        std::vector<DataGraphVertexHandle>{pivot_candidate_handle});
      assert(pivoted_candidate_set_ret);
      assert(pivoted_candidate_set.size() == gor.pattern().CountVertex());
      
      bool illegal = false;
      for (size_t rhs_vertex_idx = 0;
                  rhs_vertex_idx < rhs_vertex_handle_set.size();
                  rhs_vertex_idx++) {
        const auto& rhs_vertex_handle = rhs_vertex_handle_set[rhs_vertex_idx];
        assert(rhs_vertex_idx < rhs_vertex_candidate_set.size());
        auto& rhs_vertex_candidate = rhs_vertex_candidate_set[rhs_vertex_idx];
        assert(std::is_sorted(rhs_vertex_candidate.begin(),
                              rhs_vertex_candidate.end()));
        assert( pivoted_candidate_set.find(rhs_vertex_handle)
             != pivoted_candidate_set.end() );
        assert(!pivoted_candidate_set.find(rhs_vertex_handle)->second.empty());
        std::vector<DataGraphVertexHandle> difference_candidate_set;
        std::set_difference(pivoted_candidate_set[rhs_vertex_handle].begin(), 
                            pivoted_candidate_set[rhs_vertex_handle].end(),
                            rhs_vertex_candidate.begin(), 
                            rhs_vertex_candidate.end(),
                            std::back_inserter(difference_candidate_set));
        if (difference_candidate_set.empty()) {
          illegal = true;
          break;
        }
        pivoted_candidate_set[rhs_vertex_handle].swap(difference_candidate_set);
      }

      if (illegal) {
        continue;
      }
      
      if (!GUNDAM::_dp_iso_using_match
                 ::RefineCandidateSet(gor.pattern(),
                                         data_graph, 
                              pivoted_candidate_set, pivot)) {
        continue;
      }

      for (size_t rhs_vertex_idx = 0;
                  rhs_vertex_idx < rhs_vertex_handle_set.size();
                  rhs_vertex_idx++) {
        const auto& rhs_vertex_handle = rhs_vertex_handle_set[rhs_vertex_idx];
        assert(rhs_vertex_idx < rhs_vertex_candidate_set.size());
        auto& rhs_vertex_candidate = rhs_vertex_candidate_set[rhs_vertex_idx];
        assert(std::is_sorted(rhs_vertex_candidate.begin(),
                              rhs_vertex_candidate.end()));
        assert( pivoted_candidate_set.find(rhs_vertex_handle)
             != pivoted_candidate_set.end() );
        assert(!pivoted_candidate_set.find(rhs_vertex_handle)->second.empty());
        std::vector<DataGraphVertexHandle> merged_candidate_set;
        std::set_union(rhs_vertex_candidate.begin(), 
                       rhs_vertex_candidate.end(),
                       pivoted_candidate_set[rhs_vertex_handle].begin(), 
                       pivoted_candidate_set[rhs_vertex_handle].end(),
                       std::back_inserter(merged_candidate_set));
        rhs_vertex_candidate.swap(merged_candidate_set);
      }
    }
      
    size_t pivoted_count = 1;
    for (const auto& dst_candidate_set : rhs_vertex_candidate_set) {
      pivoted_count *= dst_candidate_set.size();
    }

    maximal_match_set.emplace_back();
    for (size_t rhs_vertex_idx = 0;
                rhs_vertex_idx < rhs_vertex_handle_set.size();
                rhs_vertex_idx++) {
      auto [ maximal_match_set_back_it,
             maximal_match_set_back_ret ]
           = maximal_match_set.back().emplace(std::move(rhs_vertex_handle_set   [rhs_vertex_idx]),
                                              std::move(rhs_vertex_candidate_set[rhs_vertex_idx]));
      assert(maximal_match_set_back_ret);
    }
    count += pivoted_count;
  }

  return count;
}


template <typename GraphPatternType, 
          typename    DataGraphType>
size_t IncrementalGORMatch(    GraphOracleRule<GraphPatternType, 
                                                  DataGraphType>& gor,
           const typename GUNDAM::VertexHandle<GraphPatternType>::type& pivot,
                                                  DataGraphType& data_graph,
  const 
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>& delta_target_graph,
     std::vector<
        std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>>& maximal_match_set) {

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using CandidateSet = std::map<GraphPatternVertexHandle,
                       std::vector<DataGraphVertexHandle>>;  

  assert(maximal_match_set.empty());
  maximal_match_set.clear();

  CandidateSet candidate_set;

  auto res = GUNDAM::Simulation<GUNDAM::MatchSemantics::kDualSimulation>(
                    gor.pattern(), 
                       data_graph,
                    candidate_set);

  if (res == 0){
    return 0;
  }  

  return IncrementalGORMatch(gor, pivot, data_graph, 
                                      candidate_set, 
                                 delta_target_graph, 
                                  maximal_match_set);
}

template <typename GraphPatternType, 
          typename    DataGraphType>
size_t IncrementalGORMatch(    GraphOracleRule<GraphPatternType, 
                                                  DataGraphType>& gor,
           const typename GUNDAM::VertexHandle<GraphPatternType>::type& pivot,
                                                  DataGraphType& data_graph,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& candidate_set,
  const 
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>& delta_target_graph) {

  std::vector<
     std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
  std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>> maximal_match_set;

  return IncrementalGORMatch(gor, pivot, data_graph, 
                                      candidate_set, 
                                 delta_target_graph, 
                                  maximal_match_set);
}

}      // namespace gor

#endif  // GOR_GOR_MATCH_H_