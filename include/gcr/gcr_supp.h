#ifndef GCR_SUPP_H_
#define GCR_SUPP_H_

#include <map>
#include <vector>

#include "include/gundam/algorithm/match_semantics.h"
#include "include/gundam/algorithm/match_using_match.h"

#include "gar/gar.h"
#include "gar/gar_supp.h"

namespace gcr {

// return std::pair<x_supp, xy_supp>
template <typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GcrSupp(gar::GraphAssociationRule<GraphPatternType,
                                       DataGraphType>& gcr,
                                       DataGraphType& data_graph,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
        std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type>>& gcr_pattern_candidate_set,
    //  std::function<bool(const GUNDAM::Match<GraphPatternType,
    //                                            DataGraphType>&)>  satisfy_x_supp_callback,
    //  std::function<bool(const GUNDAM::Match<GraphPatternType,
    //                                            DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0) {

  using GraphPatternVertexHandleType = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandleType = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using CandidateSetContainer = std::map<GraphPatternVertexHandleType,
                                std::vector<DataGraphVertexHandleType>>;

  assert(time_limit          == -1.0);
  assert(time_limit_per_supp == -1.0);

  auto& pattern = gcr.pattern();

  // when has only one non-constant literal, can evaluate
  // it through the simulation
  size_t lhs_constant_literal_counter = 0,
         lhs_variable_literal_counter = 0;
  for (const auto& x_literal_ptr : gcr.x_literal_set()) {
    assert(x_literal_ptr->type() != gar::LiteralType::kNoneLiteral);
    /* ################################## *
     * ##  statistic the lhs literals  ## *
     * ################################## */
    switch (x_literal_ptr->type()) {
     case gar::LiteralType::kConstantLiteral:
      lhs_constant_literal_counter++;
      continue;
     case gar::LiteralType::kVariableLiteral:
      lhs_variable_literal_counter++;
      continue;
     default:
      // unsupported literal type
      assert(false);
      continue;
    }
  }
  size_t lhs_non_constant_literal_counter 
                  = gcr.x_literal_set().Count() 
           - lhs_constant_literal_counter;

  if ( lhs_non_constant_literal_counter != 0) {
    // cannot utilize simulation
    return gar::GarSupp(gcr, data_graph, 
                        gcr_pattern_candidate_set, time_limit,
                                                   time_limit_per_supp);
  }
  /* ############################## *
   * ##  can utilize simulation  ## *
   * ############################## */

  gar::LiteralInfo<GraphPatternType,
                      DataGraphType> y_literal_info;
  assert(y_literal_info.literal_type() == gar::LiteralType::kNoneLiteral);
  y_literal_info = (*gcr.y_literal_set().begin())->info();

  const CandidateSetContainer* filted_constant_candidate_set_ptr
                  = std::addressof(gcr_pattern_candidate_set);

  CandidateSetContainer candidate_set_satisfy_constant_rules;
  if (lhs_constant_literal_counter > 0) {
    /* ######################################## *
     * ##  filt all constant literal in lhs  ## *
     * ######################################## */
    // has constant_literal in lhs, filt 
    candidate_set_satisfy_constant_rules = gcr_pattern_candidate_set;
    for (const auto& x_literal_ptr : gcr.x_literal_set()) {
      assert(
          x_literal_ptr->type() != gar::LiteralType::kNoneLiteral);
      if (x_literal_ptr->type() != gar::LiteralType::kConstantLiteral) {
        continue;
      }

      gar::LiteralInfo<GraphPatternType,
                          DataGraphType> literal_info = x_literal_ptr->info();

      // filt out all candidate vertexes not satisfy the constant rules
      const auto x_handle     = pattern.FindVertex(literal_info.x_id());
      const auto x_attr_key   = literal_info.x_attr_key();
      const auto x_value_str  = literal_info.c_str();
      const auto x_value_type = literal_info.data_type();
        auto  x_candidate_set_it  = candidate_set_satisfy_constant_rules.find(x_handle);
      assert(x_candidate_set_it != candidate_set_satisfy_constant_rules.end());
      auto& x_candidate_set 
          = x_candidate_set_it->second;
      std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type> satisfy_literal_x_candidate_set;
      satisfy_literal_x_candidate_set.reserve(
                      x_candidate_set.size());
      for (const auto& dst_handle : x_candidate_set) {
        const auto x_dst_attr_handle = dst_handle->FindAttribute(x_attr_key);
        if (!x_dst_attr_handle) {
          // does not has attribute, not satisfy this literal
          continue;
        }
        // has this attribute
        if (x_dst_attr_handle->value_type() != x_value_type
          || x_dst_attr_handle->value_str()  != x_value_str) {
          // has different type or value, not satisfy the literal
          continue;
        }
        // satisfy literal
        satisfy_literal_x_candidate_set.emplace_back(dst_handle);
      }
      if (satisfy_literal_x_candidate_set.size() == 0) {
        // illegal pattern, cannot satisfy all literals in lhs
        return std::pair(0,0);
      }
      x_candidate_set.swap(satisfy_literal_x_candidate_set);
    }
    GUNDAM::_dp_iso_using_match
          ::RefineCandidateSet(pattern, data_graph, 
                               candidate_set_satisfy_constant_rules);
    filted_constant_candidate_set_ptr = std::addressof(candidate_set_satisfy_constant_rules);
  }
  const CandidateSetContainer& filted_constant_candidate_set
                            = *filted_constant_candidate_set_ptr;
  
  if (y_literal_info.literal_type() == gar::LiteralType::kConstantLiteral) {
    // all literals are constant literal
    assert(y_literal_info.literal_type() == gar::LiteralType::kConstantLiteral);
    const auto x_handle = pattern.FindVertex(y_literal_info.x_id());
    assert(x_handle);
     auto  x_candidate_set_it  = filted_constant_candidate_set.find(x_handle);
    assert(x_candidate_set_it != filted_constant_candidate_set.end());
     auto& x_candidate_set 
         = x_candidate_set_it->second;
    if (x_candidate_set.empty()) {
      return std::pair(0, 0);
    }
     auto  x_supp_counter = x_candidate_set.size();
    assert(x_supp_counter > 0);
    size_t xy_supp_counter = 0;
    const auto x_attr_key   = y_literal_info.x_attr_key();
    const auto x_value_str  = y_literal_info.c_str();
    const auto x_value_type = y_literal_info.data_type();
    for (const auto& dst_handle : x_candidate_set) {
      const auto x_dst_attr_handle = dst_handle->FindAttribute(x_attr_key);
      if (!x_dst_attr_handle) {
        // does not has attribute, not satisfy this literal
        continue;
      }
      if (x_dst_attr_handle->value_type() != x_value_type
       || x_dst_attr_handle->value_str()  != x_value_str) {
        // has different type or value, not satisfy the y_literal_info
        continue;
      }
      // satisfy the y_literal_info
      xy_supp_counter++;
    }
    assert(x_supp_counter > 0);
    return std::pair(xy_supp_counter, ((float) xy_supp_counter) 
                                    / ((float)  x_supp_counter));
  }
  switch (y_literal_info.literal_type()) {
   case gar::LiteralType::kVariableLiteral: {
     auto  x_handle = pattern.FindVertex(y_literal_info.x_id());
    assert(x_handle);
     auto  y_handle = pattern.FindVertex(y_literal_info.y_id());
    assert(y_handle);

    const auto& x_attr_key = y_literal_info.x_attr_key();
    const auto& y_attr_key = y_literal_info.y_attr_key();

    assert( filted_constant_candidate_set.find(x_handle)
         != filted_constant_candidate_set.end() );
    const auto&  x_candidate_set 
          = filted_constant_candidate_set.find(x_handle)->second;

    assert( filted_constant_candidate_set.find(y_handle)
         != filted_constant_candidate_set.end() );
    const auto&  y_candidate_set 
          = filted_constant_candidate_set.find(y_handle)->second;

    uint64_t x_supp_counter = x_candidate_set.size()
                            * y_candidate_set.size();

    if (x_supp_counter == 0) {
      return std::pair(0, 0);
    }

    std::map<std::string,  // value_str
                  typename DataGraphType
                      ::VertexCounterType> x_histogram,
                                           y_histogram;

    for (const auto& x_candidate : x_candidate_set) {
      const auto x_attr_handle = x_candidate->FindAttribute(x_attr_key);
      if (!x_attr_handle) {
        // does not have this literal
        continue;
      }
      x_histogram[x_attr_handle->value_str()]++;
    }
    for (const auto& y_candidate : y_candidate_set) {
      const auto y_attr_handle = y_candidate->FindAttribute(y_attr_key);
      if (!y_attr_handle) {
        // does not have this literal
        continue;
      }
      y_histogram[y_attr_handle->value_str()]++;
    }
    /* ##################################### *
     * ## wenzhi: optimize me             ## *
     * ##     compare the two histograms  ## *
     * ##     in the merge sort manner    ## *
     * ##################################### */
    uint64_t xy_supp_counter = 0;
    for (const auto& [x_str, x_counter] : x_histogram) {
      auto y_histogram_it =  y_histogram.find(x_str);
      if ( y_histogram_it == y_histogram.end() ) {
        continue;
      }
      xy_supp_counter += x_counter * y_histogram_it->second;
    }
    assert(x_supp_counter >  0);
    assert(x_supp_counter >= xy_supp_counter);
    return std::pair(x_supp_counter,
                    xy_supp_counter);
    }
   case gar::LiteralType::kEdgeLiteral: {
     auto  x_handle = pattern.FindVertex(y_literal_info.x_id());
    assert(x_handle);
     auto  y_handle = pattern.FindVertex(y_literal_info.y_id());
    assert(y_handle);

    assert( filted_constant_candidate_set.find(x_handle)
         != filted_constant_candidate_set.end() );
    const auto& x_vertex_x_supp_candidate_set 
              = filted_constant_candidate_set.find(x_handle)->second;

    assert( filted_constant_candidate_set.find(y_handle)
          != filted_constant_candidate_set.end() );
    const auto& y_vertex_x_supp_candidate_set 
              = filted_constant_candidate_set.find(y_handle)->second;

    uint64_t x_supp_counter = 0,
            xy_supp_counter = 0;

    for (const auto& vertex_handle : x_vertex_x_supp_candidate_set) {
      for (auto out_edge_it = vertex_handle->OutEdgeBegin();
               !out_edge_it.IsDone();
                out_edge_it++) {
        if (out_edge_it->label() != y_literal_info.edge_label()) {
          continue;
        }
        if (out_edge_it->dst_handle()->label() != y_handle->label()) {
          continue;
        }
        if (!std::binary_search(y_vertex_x_supp_candidate_set.begin(),
                                y_vertex_x_supp_candidate_set.end(),
                                out_edge_it->dst_handle())) {
          continue;
        }
        xy_supp_counter++;
      }
    }
    assert(x_supp_counter > 0);
    return std::pair(x_supp_counter,
                    xy_supp_counter);
    }
   default:
    assert(false);
    break;
  }
  assert(false);
  return std::pair(0, 0);
}

template <typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GcrSupp(gar::GraphAssociationRule<GraphPatternType,
                                       DataGraphType>& gcr,
                                       DataGraphType& data_graph) {

  using GraphPatternVertexHandleType = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandleType = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using CandidateSetContainer = std::map<GraphPatternVertexHandleType,
                                std::vector<DataGraphVertexHandleType>>;

  CandidateSetContainer gcr_pattern_candidate_set;
  if (!GUNDAM::_dp_iso_using_match::InitCandidateSet<
       GUNDAM::MatchSemantics::kHomomorphism>(gcr.pattern(),
                                              data_graph,
                                              gcr_pattern_candidate_set)) {
    return std::pair(0, 0);
  }
  if (!GUNDAM::_dp_iso_using_match
             ::RefineCandidateSet(gcr.pattern(),
                                  data_graph,
                                  gcr_pattern_candidate_set)) {
    return std::pair(0, 0);
  }
  return GcrSupp(gcr, data_graph, gcr_pattern_candidate_set);
}

}; // gcr

#endif // GCR_SUPP_H_