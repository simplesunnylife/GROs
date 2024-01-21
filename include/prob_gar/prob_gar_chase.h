#ifndef _PROB_GAR_PROB_GAR_CHASE_H
#define _PROB_GAR_PROB_GAR_CHASE_H

#include "prob_gar.h"

#include "include/gundam/algorithm/match_using_match.h"

#include "include/gundam/match/match.h"

#include "include/gundam/type_getter/edge_label.h"
#include "include/gundam/type_getter/vertex_id.h"

#include "include/gundam/tool/map_edge_to.h"
#include "include/gundam/tool/has_edge.h"

#include "include/gar/literal_info.h"

#include <vector>
#include <stack>
#include <map>
#include <set>
#include <tuple>

#include "gundam/type_getter/edge_id.h"

namespace prob_gar {

namespace _prob_gar_chase {

template<typename GraphPatternType,
         typename    DataGraphType>
float Dfs(const std::vector<ProbGraphAssociationRule<GraphPatternType, 
                                                        DataGraphType>> &prob_gar_list,
          // to store how an edge can be chased by different rules
          const std::map<typename GUNDAM::EdgeID<DataGraphType>::type,
                         std::vector<std::tuple<size_t, size_t>>
                        >& edge_id_chased_by_gars,
          // to store how many chased edges are utilized for the chasing of 
          // the gars
          const std::map<std::tuple<size_t, size_t>,
                         std::vector<typename GUNDAM::EdgeID<DataGraphType>::type>
                        >& gar_used_chased_edge_id_set,
          const typename GUNDAM::EdgeID<DataGraphType>::type& edge_id,
           std::set<std::tuple<size_t, size_t>>& match_in_path) {

  // this edge should have already contained in the edge_id_chased_by_gars
  assert(edge_id_chased_by_gars.find(edge_id)
      != edge_id_chased_by_gars. end());

  const auto& match_set_for_edge_id = edge_id_chased_by_gars.find(edge_id)->second;

  std::vector<float> prob_set;

  for (const auto& match_idx : match_set_for_edge_id) {
    if (match_in_path.find(match_idx) 
     != match_in_path.end()) {
      // this edge has already been considered
      continue;
    }
    const auto& prob_gar_idx = std::get<0>(match_idx);
    assert(prob_gar_idx >= 0
        && prob_gar_idx <  prob_gar_list.size());
    // this match have not been considered before 
    // add this match_idx to match_in_path
    auto [ match_in_path_it,
           match_in_path_ret ] = match_in_path.emplace(match_idx);
    assert(match_in_path_ret);

    // std::cout << "##\tgar_idx:" << prob_gar_idx 
    //           << "\tmatch_idx:" << std::get<1>(match_idx)
    //           << std::endl;

    auto gar_used_chased_edge_id_set_it 
       = gar_used_chased_edge_id_set.find(match_idx);
    assert(gar_used_chased_edge_id_set_it
        != gar_used_chased_edge_id_set.end());

    const auto& used_chased_edge_id
          = gar_used_chased_edge_id_set_it->second;

    float prob = 1.0;

    // this match have already been added 
    for (const auto& edge_id_for_chase : used_chased_edge_id) {
      prob *= Dfs<GraphPatternType,
                     DataGraphType>(prob_gar_list,
                                    edge_id_chased_by_gars,
                                  gar_used_chased_edge_id_set,
                                    edge_id_for_chase,
                                    match_in_path);
    } 

    prob_set.emplace_back(prob * prob_gar_list[prob_gar_idx].prob());

    // erase this match
    match_in_path.erase(match_in_path_it);
  }

  float neg_prob = 1.0;

  for (const auto& prob : prob_set) {
    neg_prob *= (1.0 - prob);
  }

  return 1.0 - neg_prob;
}

};

// forward inference
template <typename GraphPatternType, 
          typename    DataGraphType>
int ProbGARChase(
    const std::vector<ProbGraphAssociationRule<GraphPatternType, 
                                                  DataGraphType>> &prob_gar_list, 
    DataGraphType &data_graph,
    std::map<typename GUNDAM::EdgeID<DataGraphType>::type, float>& diff_edge_set) {

  // std::cout << "prob_gar_list.size(): "
  //           <<  prob_gar_list.size() << std::endl;

  using GraphPatternVertexHandleType = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandleType = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using GraphPatternEdgeHandleType = typename GUNDAM::EdgeHandle<GraphPatternType>::type;
  using    DataGraphEdgeHandleType = typename GUNDAM::EdgeHandle<   DataGraphType>::type;

  using  VertexIDType = typename GUNDAM:: VertexID<GraphPatternType>::type;
  using EdgeLabelType = typename GUNDAM::EdgeLabel<GraphPatternType>::type;

  using DataGraphVertexIDType = typename GUNDAM::VertexID<DataGraphType>::type;

  using DataGraphEdgeIDType    = typename GUNDAM::EdgeID   <DataGraphType>::type;
  using DataGraphEdgeLabelType = typename GUNDAM::EdgeLabel<DataGraphType>::type;

  using ProbGarType = ProbGraphAssociationRule<GraphPatternType,
                                                  DataGraphType>;

  using MatchType = GUNDAM::Match<const GraphPatternType,
                                           DataGraphType>;

  using MatchSetType = GUNDAM::MatchSet<const GraphPatternType,
                                                 DataGraphType>;

  #ifndef NDEBUG
  for (const auto& prob_gar 
                 : prob_gar_list) {
    if (prob_gar.y_literal_set().Count() == 1
   && (*prob_gar.y_literal_set().begin())->type() == gar::LiteralType::kEdgeLiteral) {
      // can only process prob_gar with edge literal rhs
      continue;
    }
    assert(false);
    return -1;
  }
  #endif // NDEBUG

  DataGraphEdgeIDType max_edge_id = 0;

  std::vector<DataGraphEdgeIDType> new_edge_id_set,
                              original_edge_id_set;
   
  for (auto vertex_it = data_graph.VertexBegin();
           !vertex_it.IsDone();
            vertex_it++) {
    for (auto out_edge_it = vertex_it->OutEdgeBegin();
             !out_edge_it.IsDone();
              out_edge_it++) {
      max_edge_id = max_edge_id > out_edge_it->id()?
                    max_edge_id : out_edge_it->id();
           new_edge_id_set.emplace_back(out_edge_it->id());
      original_edge_id_set.emplace_back(out_edge_it->id());
    }
  }

  std::sort(new_edge_id_set.begin(), 
            new_edge_id_set.end());

  std::sort(original_edge_id_set.begin(), 
            original_edge_id_set.end());

  DataGraphEdgeIDType edge_id_allocator = max_edge_id;
  edge_id_allocator++;

  std::vector<MatchSetType> prob_gar_match_set_list;
  prob_gar_match_set_list.resize(prob_gar_list.size());

  // to store how an edge can be chased by different rules
  std::map<DataGraphEdgeIDType,
           std::vector<std::tuple<size_t, // prob_gar idx 
                                  size_t> // match idx 
                      >
          > edge_id_chased_by_gars;

  // to store how many chased edges are utilized for the chasing of 
  // the gars
  std::map<std::tuple<size_t,  // prob_gar idx 
                      size_t>, // match idx 
           std::vector<DataGraphEdgeIDType>
          > gar_used_chased_edge_id_set;

  bool is_first_round = true;
  // keep matching if there are new edges added in the
  // last round
  while (!new_edge_id_set.empty()) {
    // since adding an edge to pattern can change the edge handle
    // as well as the vertex handle, store all edges that need
    // to be added and add them together after the matching
    std::map<std::tuple<DataGraphVertexIDType,
                        DataGraphEdgeLabelType,
                        DataGraphVertexIDType>, 
             DataGraphEdgeIDType> edge_to_add;

    for (size_t prob_gar_idx = 0;
                prob_gar_idx < prob_gar_list.size();
                prob_gar_idx++) {
      auto& prob_gar           = prob_gar_list          [prob_gar_idx];
      auto& prob_gar_match_set = prob_gar_match_set_list[prob_gar_idx];
      const auto rhs_literal_info = (*prob_gar.y_literal_set().begin())->info();
      assert(rhs_literal_info.literal_type() == gar::LiteralType::kEdgeLiteral);

      std::function<bool(const MatchType&)> prob_gar_prune_callback 
        = [](const MatchType& match) -> bool {
        // prune nothing
        return false;
      };

      std::function<bool(const MatchType&)> prob_gar_match_callback 
       = [&prob_gar_idx,
          &prob_gar,
          &prob_gar_match_set,
             &new_edge_id_set,
          #ifndef NDEBUG
          &data_graph,
          #endif // NDEBUG
          &is_first_round,
          &original_edge_id_set,
          &rhs_literal_info,
          &edge_id_chased_by_gars,
          &edge_id_allocator,
          &edge_to_add,
          &gar_used_chased_edge_id_set](const MatchType& match) -> bool {

        // set of the edge ids that are used in this match
        std::vector<DataGraphEdgeIDType> used_chased_edge_id_set;

        bool has_edge_added_in_last_round = false;

        for (auto edge_it = prob_gar.pattern().EdgeBegin();
                 !edge_it.IsDone();
                  edge_it++) {

          DataGraphEdgeHandleType dst_edge_handle = GUNDAM::MapEdgeTo(match, edge_it);
          // should have corresponding edge of this edge
          assert(dst_edge_handle);

          // to find whether this edge has already been contained
          // in the original graph
          if (!std::binary_search(original_edge_id_set.begin(),
                                  original_edge_id_set. end (),
                                       dst_edge_handle->id())) {
            // this edge has not been contained in the original data graph
            used_chased_edge_id_set.emplace_back(dst_edge_handle->id());
          }

          if (has_edge_added_in_last_round) {
            // does not need to find it in edge_it_set again
            continue;
          }
          if (!std::binary_search(new_edge_id_set.begin(),
                                  new_edge_id_set. end (),
                                  dst_edge_handle->id())) {
            // this edge is not added in last round
            continue;
          }
          // this edge is added in last round
          has_edge_added_in_last_round = true;
        }

        if (!is_first_round
         && !has_edge_added_in_last_round) {
          // does not match to new added edges, does not need to add 
          // this match in the match set
          // should have been 
          return true;
        }

        // this match contains new edges added in this round

        // to find whether the generated rules are contained 
        // in the orignal data graph or added in previous rounds
        DataGraphVertexHandleType x_handle = match.MapTo(rhs_literal_info.x_id());
        assert(x_handle);
        DataGraphVertexHandleType y_handle = match.MapTo(rhs_literal_info.y_id());
        assert(y_handle);

        for (auto out_edge_it = x_handle->OutEdgeBegin();
                 !out_edge_it.IsDone();
                  out_edge_it++) {
          if (out_edge_it->label() != rhs_literal_info.edge_label()) {
            continue;
          }
          if (out_edge_it->dst_handle() != y_handle) {
            continue;
          }
          // this edge exists
          // to find whether the chased edge contained in the original data graph
          if (std::binary_search(original_edge_id_set.begin(),
                                 original_edge_id_set. end (),
                                      out_edge_it->id())) {
            // this edge is contained in the orignal data graph
            #ifndef NDEBUG
            // all edge should not be contained in the original graph
            for (auto out_edge_it = x_handle->OutEdgeBegin();
                     !out_edge_it.IsDone();
                      out_edge_it++) {
              if (out_edge_it->label() != rhs_literal_info.edge_label()) {
                continue;
              }
              if (out_edge_it->dst_handle() != y_handle) {
                continue;
              }
              // this edge should be contained in the original graph
              assert(std::binary_search(original_edge_id_set.begin(),
                                        original_edge_id_set. end (),
                                            out_edge_it->id()));
            }
            if constexpr (GUNDAM::GraphParameter<DataGraphType>::vertex_level_count_edge
                       && GUNDAM::GraphParameter<DataGraphType>::vertex_level_edge_label_index
                       && GUNDAM::GraphParameter<DataGraphType>::duplicate_edge) {
              assert(x_handle->CountOutEdge(rhs_literal_info.edge_label(),
                                            y_handle) > 1);
            }
            #endif // NDEBUG
            return true;
          }
          // has not been contained in original data graph
          #ifndef NDEBUG
          // all edge should not be contained in the original graph
          for (auto out_edge_it = x_handle->OutEdgeBegin();
                   !out_edge_it.IsDone();
                    out_edge_it++) {
            if (out_edge_it->label() != rhs_literal_info.edge_label()) {
              continue;
            }
            if (out_edge_it->dst_handle() != y_handle) {
              continue;
            }
            // this edge should not be contained in the original graph
            assert(!std::binary_search(original_edge_id_set.begin(),
                                       original_edge_id_set. end (),
                                            out_edge_it->id()));
          }
          if constexpr (GUNDAM::GraphParameter<DataGraphType>::vertex_level_count_edge
                     && GUNDAM::GraphParameter<DataGraphType>::vertex_level_edge_label_index
                     && GUNDAM::GraphParameter<DataGraphType>::duplicate_edge) {
            // this edge is not contained in the original data graph, 
            // should have been chased only once
            assert(x_handle->CountOutEdge(rhs_literal_info.edge_label(),
                                          y_handle) == 1);
          }
          #endif // NDEBUG

          // the chased edge have already existed
          // but not contained in the original graph

          // add to the match set
          prob_gar_match_set.AddMatch(match);
          // this edge should have been chased by other gars
          assert(edge_id_chased_by_gars.find(out_edge_it->id())
              != edge_id_chased_by_gars.end());
          // add to the match result
          edge_id_chased_by_gars[out_edge_it->id()].emplace_back(
                       prob_gar_idx,
                       prob_gar_match_set.size() - 1);
          assert(gar_used_chased_edge_id_set.find(std::tuple(prob_gar_idx,
                                                             prob_gar_match_set.size() - 1))
              == gar_used_chased_edge_id_set.end());
          auto [ gar_used_chased_edge_id_set_it,
                 gar_used_chased_edge_id_set_ret ] 
               = gar_used_chased_edge_id_set.emplace(std::tuple(prob_gar_idx,
                                                                prob_gar_match_set.size() - 1),
                                                     std::move(used_chased_edge_id_set));
          // should have been added successfully
          assert(gar_used_chased_edge_id_set_ret);
          // continue matching
          return true;
        }
        // does not have such an edge, needs to add it 

        // to find whether this edge have already been added
        // before in this round
        auto edge_to_add_it = edge_to_add.find(std::tuple(x_handle->id(),
                                           rhs_literal_info.edge_label(),
                                                          y_handle->id()));

        #ifndef NDEBUG
        if (edge_to_add_it != edge_to_add.end()) {
          // this edge has not been added before in this round

          // should have already been added before
          assert(edge_id_chased_by_gars.find(edge_to_add_it->second) 
              != edge_id_chased_by_gars.end());

          // should not have been contained in data graph
          assert(!data_graph.FindEdge(edge_to_add_it->second));
        }
        #endif // NDEBUG

        if (edge_to_add_it == edge_to_add.end()) {
          // this edge has not been added before in this round yet
          edge_to_add_it = edge_to_add.emplace_hint(edge_to_add_it,
                                               std::tuple(x_handle->id(),
                                           rhs_literal_info.edge_label(),
                                                          y_handle->id()),
                                                    edge_id_allocator++);
          // should have been added successfully
          assert(edge_to_add_it != edge_to_add.end());
          assert(edge_to_add_it->second + 1 == edge_id_allocator);
        }

        assert(edge_to_add_it != edge_to_add.end());
        prob_gar_match_set.AddMatch(match);

        edge_id_chased_by_gars[edge_to_add_it->second].emplace_back(prob_gar_idx,
                                                                    prob_gar_match_set.size() - 1);
        assert(edge_id_chased_by_gars.find(edge_to_add_it->second)
            != edge_id_chased_by_gars. end());
        // >= 1 since this edge might also be chased by other gars in this round
        assert(edge_id_chased_by_gars.find(edge_to_add_it->second)->second.size() >= 1);
        assert(std::get<0>(edge_id_chased_by_gars.find(edge_to_add_it->second)->second.back())
                == prob_gar_idx);
        assert(std::get<1>(edge_id_chased_by_gars.find(edge_to_add_it->second)->second.back())
                == prob_gar_match_set.size() - 1);

        auto [ gar_used_chased_edge_id_set_it,
               gar_used_chased_edge_id_set_ret ] 
             = gar_used_chased_edge_id_set.emplace(std::tuple(prob_gar_idx,
                                                              prob_gar_match_set.size() - 1),
                                                    std::move(used_chased_edge_id_set));
        // should have been added successfully
        assert(gar_used_chased_edge_id_set_ret);
        // continue matching
        return true;
      };
      
      GUNDAM::MatchUsingMatch(prob_gar.pattern(),
                              data_graph,
                              prob_gar_prune_callback,
                              prob_gar_match_callback);
    }

    new_edge_id_set.clear();
    for (const auto& [edge_type, edge_id] : edge_to_add) {
      auto [edge_handle,
            edge_ret] = data_graph.AddEdge(std::get<0>(edge_type),
                                           std::get<2>(edge_type),
                                           std::get<1>(edge_type),
                                                       edge_id);
      assert(edge_ret);
      new_edge_id_set.emplace_back(edge_id);
    }
    std::sort(new_edge_id_set.begin(), 
              new_edge_id_set.end());
    // match set no-longer works
    is_first_round = false;

    #ifndef NDEBUG
    for (const auto& [edge_id, match_set] : edge_id_chased_by_gars) {
      // edge should have been added in the data graph
      assert(data_graph.FindEdge(edge_id));
    }
    #endif // NDEBUG
  }

  #ifndef NDEBUG
  for (const auto& [idx, used_chased_edge_id_set]
                   : gar_used_chased_edge_id_set) {
    const auto& prob_gar_idx = std::get<0>(idx);
    assert(used_chased_edge_id_set.size() >= 0);
    assert(prob_gar_list[prob_gar_idx].pattern().CountEdge()
        >= used_chased_edge_id_set.size());
  }
  #endif // NDEBUG

  assert(original_edge_id_set.size() 
                + edge_id_chased_by_gars.size() == data_graph.CountEdge());


  // std::cout << "#############" << std::endl;
  // for (const auto& [edge_id,
  //                  match_idx_set] : edge_id_chased_by_gars) {
  //   std::cout << "edge_id: " << edge_id << std::endl;
  //   for (const auto& [gar_idx, match_idx] : match_idx_set) {
  //     std::cout << "\tgar_idx:" << gar_idx
  //               << "\tmatch_idx:" << match_idx
  //               << std::endl;
  //   }
  // }

  // std::cout << "#############" << std::endl;
  // for (const auto& [match_idx,
  //                    edge_id_set] : gar_used_chased_edge_id_set) {
  //   std::cout << "gar_idx: " << std::get<0>(match_idx)
  //         << "\tmatch_idx: " << std::get<1>(match_idx) << std::endl;
  //   for (const auto& edge_id : edge_id_set) {
  //     std::cout << "\tedge_id:" << edge_id << std::endl;
  //   }
  // }

  // dfs searching through the match set
  for (const auto& [edge_id, match_set] : edge_id_chased_by_gars) {

    std::set<std::tuple<size_t, size_t>> match_in_path;

    diff_edge_set.emplace(edge_id, _prob_gar_chase::Dfs<GraphPatternType,
                                                           DataGraphType>(
                                                        prob_gar_list,
                                                        edge_id_chased_by_gars,
                                                       gar_used_chased_edge_id_set,
                                                        edge_id,
                                                        match_in_path));
  }
  
  return diff_edge_set.size();
}

template <class GraphPatternType, 
          class    DataGraphType>
int ProbGARChase(
    const ProbGraphAssociationRule<GraphPatternType, 
                                      DataGraphType> &prob_gar, 
    DataGraphType &data_graph,
    std::map<typename GUNDAM::EdgeID<DataGraphType>::type, float>& diff_edge_set) {
  std::vector<ProbGraphAssociationRule<GraphPatternType, 
                                          DataGraphType>> prob_gar_list;
  prob_gar_list.emplace_back(prob_gar);
  return ProbGARChase(prob_gar_list, data_graph, diff_edge_set);
}

}; // namespace prob_gar

#endif  // _PROB_GAR_PROB_GAR_CHASE_H