#ifndef GAR_SUPP_H_
#define GAR_SUPP_H_

#include <map>

#include "gar/gar.h"
#include "gar/same_gar.h"

#include "gundam/algorithm/match_using_match.h"

#include "gundam/match/match.h"

#include "gundam/type_getter/vertex_handle.h"
#include "gundam/type_getter/edge_handle.h"
#include "gundam/tool/operator/preserve.h"
#include "gundam/tool/sub_graph_of.h"

#include "gundam/io/csvgraph.h"

namespace gar{
namespace _gar_supp{

template <typename LiteralContainerType,
          typename          VertexIDType>
void AddVertexIDTo(LiteralContainerType& literal_container,
                  std::set<VertexIDType>& vertex_id_set){
  for (const auto& literal_ptr : literal_container) {
    auto literal_info = literal_ptr->info();
    switch (literal_info.literal_type()){
    case gar::LiteralType::kAttrValueLiteral:
      // x.A
    case gar::LiteralType::kConstantLiteral:
      // x.A == c
      vertex_id_set.emplace(literal_info.x_id());
      continue;

    case gar::LiteralType::kEdgeLiteral:
      // x_id_ - edge_label_ -> y_id_
    case gar::LiteralType::kVariableLiteral:
      // x.A == y.B
    case gar::LiteralType::kMlLiteral:
      // x_id_ - (module_url_: edge_label_) -> y_id_
      vertex_id_set.emplace(literal_info.x_id());
      vertex_id_set.emplace(literal_info.y_id());
      continue;
    case gar::LiteralType::kKeyLiteral:
      // x_id_ == y_id_
      vertex_id_set.emplace(literal_info.x_id());
      vertex_id_set.emplace(literal_info.y_id());
      continue;
    default:
      // unknown literal type
      std::cout << "unknown literal type" << std::endl;
      assert(false);
      break;
    }
  }
  return;
}

/* ############################################################# *
 * ##  this function is utilized to grep the sub-pattern      ## *
 * ##  that are contained in the literals and user-specified  ## *
 * ##  pivot set                                              ## *     
 * ##  for example, consider the following example gar:       ## *
 * ##     pattern:  0 -> 1 -> 2 -> 3                          ## *
 * ##           X:  0.A = a0                                  ## *
 * ##           Y:  2.A = a1                                  ## *
 * ##   when pivoted_vertex_set is empty                      ## *
 * ##     rhs_only == true                                    ## *
 * ##       return 2  (isolated vertex)                       ## *
 * ##     rhs_only == false                                   ## *
 * ##       return 0   2  (two isolated vertexex)             ## *
 * ##                                                         ## *
 * ##   when pivoted_vertex_set == {1}                        ## *
 * ##     rhs_only == true                                    ## *
 * ##       return 1 -> 2                                     ## *
 * ##     rhs_only == false                                   ## *
 * ##       return 0 -> 1 -> 2                                ## *
 * ############################################################# */
template <typename GraphPatternType,
          typename    DataGraphType>
inline GraphPatternType LiteralPattern(
  GraphAssociationRule<GraphPatternType,
                          DataGraphType>& gar, bool rhs_only,
  const std::vector<typename 
  GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set
      = std::vector<typename 
  GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  using VertexIDType = typename GUNDAM::VertexID<GraphPatternType>::type;

  std::set<VertexIDType> preserve_vertex_id_set;

  AddVertexIDTo(gar.y_literal_set(), preserve_vertex_id_set);
  if (!rhs_only) {
    AddVertexIDTo(gar.x_literal_set(), preserve_vertex_id_set);
  }
  
  for (const auto& vertex_handle : pivoted_vertex_set) {
    preserve_vertex_id_set.emplace(vertex_handle->id());
  }

  std::set<VertexIDType> erase_vertex_id_set;
  for (auto vertex_it = gar.pattern().VertexBegin();
           !vertex_it.IsDone();
            vertex_it++) {
    if (preserve_vertex_id_set.find(vertex_it->id())
     != preserve_vertex_id_set.end()){
      // contained in preserve_vertex_id_set
      continue;
    }
    erase_vertex_id_set.emplace(vertex_it->id());
  }
  auto literal_pattern(gar.pattern());
  for (const auto& erase_vertex_id : erase_vertex_id_set){
    literal_pattern.EraseVertex(erase_vertex_id);
  }
  return literal_pattern;
}

template <typename GraphPatternType,
          typename    DataGraphType>
inline GraphAssociationRule<GraphPatternType,
                               DataGraphType> 
RhsGar(GraphAssociationRule<GraphPatternType,
                               DataGraphType>& gar,
       const std::vector<typename 
       GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set
           = std::vector<typename 
       GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
    pivoted_rhs_gar(_gar_supp::LiteralPattern(gar, true, pivoted_vertex_set));
  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    gar_pattern_to_rhs_pattern(gar.pattern(), 
                   pivoted_rhs_gar.pattern(), "same_id_map");
  assert(gar_pattern_to_rhs_pattern.size()
             == pivoted_rhs_gar.pattern().CountVertex());
  for (const auto& y_literal_ptr : gar.y_literal_set()) {
    #ifndef NDEBUG
    std::vector<GraphPatternVertexHandle> y_literal_vertex_handle_set;
    y_literal_ptr->CalPivot(y_literal_vertex_handle_set);
    assert(y_literal_vertex_handle_set.size() == 1
        || y_literal_vertex_handle_set.size() == 2);
    for (const auto& y_handle : y_literal_vertex_handle_set) {
      assert(gar_pattern_to_rhs_pattern.HasMap(y_handle));
    }
    #endif // NDEBUG
    pivoted_rhs_gar.AddY(y_literal_ptr->info());
  }
  for (const auto& x_literal_ptr : gar.x_literal_set()) {
    std::vector<GraphPatternVertexHandle> x_literal_vertex_handle_set;
    x_literal_ptr->CalPivot(x_literal_vertex_handle_set);
    assert(x_literal_vertex_handle_set.size() == 1
        || x_literal_vertex_handle_set.size() == 2);
    bool x_literal_mapped = true;
    for (const auto& x_handle : x_literal_vertex_handle_set) {
      if (!gar_pattern_to_rhs_pattern.HasMap(x_handle)) {
        x_literal_mapped = false;
        break;
      }
    }
    if (!x_literal_mapped) {
      continue;
    }
    // add all x literals that are contained 
    pivoted_rhs_gar.AddX(x_literal_ptr->info());
  }
  return pivoted_rhs_gar;
}

template <typename GraphPatternType,
          typename    DataGraphType>
inline GraphAssociationRule<GraphPatternType,
                               DataGraphType> 
LiteralGar(GraphAssociationRule<GraphPatternType,
                                   DataGraphType>& gar,
           const std::vector<typename 
           GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
               = std::vector<typename 
           GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
    literal_gar(_gar_supp::LiteralPattern(gar, false, pivoted_vertex_set));
  GUNDAM::Match<GraphPatternType, 
                GraphPatternType> 
    gar_pattern_to_literal_pattern(gar.pattern(), 
                           literal_gar.pattern(),
                           "same_id_map");
  assert(gar_pattern_to_literal_pattern.size()
        == literal_gar.pattern().CountVertex());
  for (const auto& y_literal_ptr   
             : gar.y_literal_set()){
    #ifndef NDEBUG
    std::vector<GraphPatternVertexHandle> y_literal_vertex_handle_set;
    y_literal_ptr->CalPivot(y_literal_vertex_handle_set);
    assert(y_literal_vertex_handle_set.size() == 1
        || y_literal_vertex_handle_set.size() == 2);
    for (const auto& y_handle : y_literal_vertex_handle_set) {
      assert(gar_pattern_to_literal_pattern.HasMap(y_handle));
    }
    #endif // NDEBUG
    literal_gar.AddY(y_literal_ptr->info());
  }
  for (const auto& x_literal_ptr 
             : gar.x_literal_set()){
    std::vector<GraphPatternVertexHandle> x_literal_vertex_handle_set;
    x_literal_ptr->CalPivot(x_literal_vertex_handle_set);
    assert(x_literal_vertex_handle_set.size() == 1
        || x_literal_vertex_handle_set.size() == 2);
    bool x_literal_mapped = true;
    for (const auto& x_handle : x_literal_vertex_handle_set) {
      if (!gar_pattern_to_literal_pattern.HasMap(x_handle)){
        x_literal_mapped = false;
        break;
      }
    }
    if (!x_literal_mapped){
      continue;
    }
    literal_gar.AddX(x_literal_ptr->info());
  }
  return std::move(literal_gar);
}

} // namespace _gar_supp

/*  ################################################  *
 *  ##                                            ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##          pivoted_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##      pivoted_rhs_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##  pivoted_literal_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##              gar_pattern_match_callback    ##  *
 *  ##                                            ##  *
 *  ################################################  */

// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
  const GUNDAM::Match<GraphPatternType,
                         DataGraphType>& gar_pattern_to_data_graph_partial_match,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_contained_in_literals,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  assert(time_limit == -1.0
      || time_limit > 0);

  assert(time_limit_per_supp == -1.0
      || time_limit_per_supp > 0);

  assert(pivoted_vertex_set.empty() 
      || pivoted_vertex_set.size() > gar_pattern_to_data_graph_partial_match.size());

  if (// time_limit is set but time_limit_per_supp is not set
      (time_limit_per_supp == -1.0 && time_limit != -1.0)
      // both time_limit and time_limit_per_supp are set
      // but there are time_limit_per_supp > time_limit
   || (time_limit_per_supp != -1.0 && time_limit != -1.0
    && time_limit_per_supp > time_limit)) {
    time_limit_per_supp = time_limit;
  }

  auto& gar_pattern = gar.pattern();

  using MatchType = GUNDAM::Match<GraphPatternType, 
                                     DataGraphType>;

  using CandidateSetContainer = std::map<GraphPatternVertexHandle, 
                                std::vector<DataGraphVertexHandle>>;

  /* ############################################################## *
   * ##   obtain pattern contained in all literals and pivotes   ## *
   * ############################################################## */
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
        pivoted_literal_gar = _gar_supp::LiteralGar(gar, pivoted_vertex_set);
  auto& pivoted_literal_pattern = pivoted_literal_gar.pattern();

  /* ############################################################## *
   * ##   obtain pattern contained in rhs literals and pivotes   ## *
   * ############################################################## */
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
         pivoted_rhs_gar = _gar_supp::RhsGar(gar, pivoted_vertex_set);
  auto&  pivoted_rhs_pattern = pivoted_rhs_gar.pattern();
  assert(pivoted_literal_pattern.CountVertex()
          >= pivoted_rhs_pattern.CountVertex());
  assert(GUNDAM::SubGraphOf(pivoted_rhs_pattern,
                        pivoted_literal_pattern));

  /* ############################################# *
   * ##   obtain pattern contained in pivotes   ## *
   * ############################################# */
  GraphPatternType pivoted_pattern = pivoted_vertex_set.empty()?
                                     pivoted_rhs_pattern
                                   : GUNDAM::PreserveVertexSet(gar.pattern(),
                                     pivoted_vertex_set);
  assert(pivoted_rhs_pattern.CountVertex() 
          >= pivoted_pattern.CountVertex());
  assert(pivoted_vertex_set.empty()
      || pivoted_pattern.CountVertex() 
      == pivoted_vertex_set.size());
  assert(GUNDAM::SubGraphOf(pivoted_pattern,
                        pivoted_rhs_pattern));
  #ifndef NDEBUG
  for (auto map_it = gar_pattern_to_data_graph_partial_match.MapBegin();
           !map_it.IsDone();
            map_it++) {
    assert(pivoted_pattern.FindVertex(map_it->src_handle()->id()));
  }
  #endif // NDEBUG

  /* ########################################## *
   * ##  parameters to be used in callbacks  ## *
   * ########################################## */
  const bool pivoted_pattern_same_as_pivoted_rhs_pattern 
               = GUNDAM::SamePattern(pivoted_pattern,
                                     pivoted_rhs_pattern);

  const bool pivoted_rhs_gar_same_as_pivoted_literal_gar 
                           = SameGar(pivoted_rhs_gar,
                                     pivoted_literal_gar);

  const bool pivoted_literal_pattern_same_as_gar_pattern
           = GUNDAM::SamePattern(pivoted_literal_pattern,
                                             gar_pattern);

  /* ############################ *
   * ##  from pivoted_pattern  ## *
   * ############################ */
  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_rhs_pattern_to_pivoted_pattern_partial_match(
    pivoted_rhs_pattern,   pivoted_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_literal_pattern_to_pivoted_pattern_partial_match(
    pivoted_literal_pattern,   pivoted_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    gar_pattern_to_pivoted_pattern_partial_match(
    gar_pattern,   pivoted_pattern, "same_id_map");

  /* ################################ *
   * ##  from pivoted_rhs_pattern  ## *
   * ################################ */
  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match(
    pivoted_literal_pattern,   pivoted_rhs_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    gar_pattern_to_pivoted_rhs_pattern_partial_match(
    gar_pattern,   pivoted_rhs_pattern, "same_id_map");

  /* #################################### *
   * ##  from pivoted_literal_pattern  ## *
   * #################################### */
  GUNDAM::Match<GraphPatternType, GraphPatternType>
    gar_pattern_to_pivoted_literal_pattern_partial_match(
    gar_pattern,   pivoted_literal_pattern, "same_id_map");

  CandidateSetContainer pivoted_pattern_candidate_set,
                    pivoted_rhs_pattern_candidate_set,
                pivoted_literal_pattern_candidate_set;

  // process pivoted_literal_pattern_candidate_set
  //         and pivoted_rhs_pattern_candidate_set
  for (const auto& [gar_pattern_candidate_ptr,
                    gar_pattern_candidate]
                  : gar_pattern_candidate_set_for_vertexes_contained_in_literals) {
    // should not be null
    assert(gar_pattern_candidate_ptr);
    const auto kPivotedLiteralPatternPtr 
              = pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
    if (kPivotedLiteralPatternPtr) {
      // this vertex is contained in pivoted_literal_pattern
      auto [ pivoted_literal_pattern_candidate_set_it,
             pivoted_literal_pattern_candidate_set_ret ]
           = pivoted_literal_pattern_candidate_set.emplace(
            kPivotedLiteralPatternPtr,
            gar_pattern_candidate);
      // should added successfully
      assert(pivoted_literal_pattern_candidate_set_ret);
      const auto kPivotedRhsPatternPtr
                = pivoted_rhs_pattern.FindVertex(gar_pattern_candidate_ptr->id());
      if (kPivotedRhsPatternPtr){
        // this vertex is contained in pivoted_rhs_pattern
        auto [ pivoted_rhs_pattern_candidate_set_it,
               pivoted_rhs_pattern_candidate_set_ret ]
             = pivoted_rhs_pattern_candidate_set.emplace(
              kPivotedRhsPatternPtr,
                gar_pattern_candidate);
        // should added successfully
        assert(pivoted_rhs_pattern_candidate_set_ret);
        const auto kPivotedPatternPtr
                  = pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        if (kPivotedPatternPtr){
          // this vertex is contained in pivoted_rhs_pattern
          auto [ pivoted_pattern_candidate_set_it,
                 pivoted_pattern_candidate_set_ret ]
               = pivoted_pattern_candidate_set.emplace(
                 kPivotedPatternPtr,
                  gar_pattern_candidate);
          // should added successfully
          assert(pivoted_pattern_candidate_set_ret);
        }
      }
      #ifndef NDEBUG
      else {
        // should not be contained in pivoted_rhs_pattern
        assert(!pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id()));
      }
      #endif // NDEBUG
    }
    #ifndef NDEBUG
    else {
      // should not be contained in pivoted_rhs_pattern and pivoted_pattern
      assert(!pivoted_rhs_pattern.FindVertex(gar_pattern_candidate_ptr->id()));
      assert(    !pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id()));
    }
    #endif // NDEBUG
  }

  bool has_satisfy_x       = false,
       has_satisfy_x_not_y = false,
       has_satisfy_x_and_y = false;

  typename DataGraphType::VertexCounterType x_supp = 0,
                                           xy_supp = 0;

  /* ################################################ *
   * ##  callbacks from gar_pattern to data graph  ## *
   * ################################################ */
  std::function<bool(const MatchType&)> 
    gar_pattern_prune_callback
    = [](const MatchType& 
    gar_pattern_to_data_graph_match) -> bool {
    // prune nothing, continue the matching
    return false;
  };
  // based on the partial match from literal_kPattern to data graph
  // complete the match of the entire gar_pattern 
  std::function<bool(const MatchType&)> 
    gar_pattern_match_callback 
    = [&has_satisfy_x](const MatchType& 
    gar_pattern_to_data_graph_match) -> bool {
    // only need to complete the match, find one match that 
    // satisfy all literal is enough
    has_satisfy_x = true;
    // terminate matching
    return false;
  };

  /* ############################################################ *
   * ##  callbacks from pivoted_literal_pattern to data graph  ## *
   * ############################################################ */
  // if the current pivoted_literal_pattern_to_data_graph_match
  // has already violate a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    pivoted_literal_pattern_prune_callback
    = [&](const MatchType& 
    pivoted_literal_pattern_to_data_graph_match) -> bool {
    for (const auto& x_literal_ptr : pivoted_literal_gar.x_literal_set()) {
      if (!x_literal_ptr->MappedBy(pivoted_literal_pattern_to_data_graph_match)){
        // this literal is not considered in the current
        // pivoted_literal_pattern_to_data_graph_match,
        // this match cannot be pruned
        //
        // move to the next literal
        continue;
      }
      if (!x_literal_ptr->Satisfy(pivoted_literal_pattern_to_data_graph_match)) {
        // does not satisfy all literals
        // does not need further matching
        return true;
      }
    }
    // satisfy all literals
    // continue matching
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_literal_pattern_match_callback
    = [&](const MatchType& 
    pivoted_literal_pattern_to_data_graph_match) -> bool {
    if (!pivoted_literal_pattern_same_as_gar_pattern) {
      // pivoted_rhs_gar is not the same as pivoted_literal_gar
      // first match the pivoted_rhs_literal, then match
      // the pivoted_literal_gar
      GUNDAM::Match<GraphPatternType, DataGraphType>  
                      gar_pattern_to_data_graph_partial_match
        = pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_pivoted_literal_pattern_partial_match);
      assert(!has_satisfy_x);
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kMerge>(
                              gar_pattern,   data_graph,
                              gar_pattern_to_data_graph_partial_match,
                              gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                              gar_pattern_prune_callback,
                              gar_pattern_match_callback,
                              time_limit_per_supp);
      if (has_satisfy_x) {
        // has found satisfy X, terminate matching
        return false;
      }
      // continue matching
      return true;
    }
    assert(!has_satisfy_x);
    // gar_pattern of literal_gar is the same as the
    // entire gar_pattern, does not need further match
    has_satisfy_x = true;
    // has found satisfy X, terminate matching
    return false;
  };

  /* ######################################################## *
   * ##  callbacks from pivoted_rhs_pattern to data graph  ## *
   * ######################################################## */
  // if the current match_state has already violate
  // a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    pivoted_rhs_pattern_prune_callback
    = [&](const MatchType&
    pivoted_rhs_pattern_to_data_graph_match) -> bool {
    for (const auto& x_literal_ptr : pivoted_rhs_gar.x_literal_set()) {
      if (!x_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_match)) {
        // this literal is not considered in the current
        // pivoted_rhs_pattern_to_data_graph_match, 
        // this match cannot be pruned move to the next literal
        continue;
      }
      if (!x_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_match)) {
        // does not satisfy all literals
        // does not need further matching
        return true;
      }
    }
    // satisfy all literals
    // continue matching
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_rhs_pattern_match_callback
    = [&](const MatchType& 
    pivoted_rhs_pattern_to_data_graph_match) -> bool {
    assert(!has_satisfy_x);
    assert(!has_satisfy_x_not_y);
    /* ############################################## *
     * ##  first check whether satisfy y literals  ## *
     * ############################################## */
    bool satisify_y_literals = true;
    for (const auto& y_literal_ptr : pivoted_rhs_gar.y_literal_set()) {
      assert(y_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_match));
      if (y_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_match)) {
        continue;
      }
      satisify_y_literals = false;
      break;
    }

    if (!pivoted_rhs_gar_same_as_pivoted_literal_gar) {
      // pivoted_rhs_gar is not the same as literal_gar
      // first match the rhs literal, then match
      // the literal_gar

      /* ####################################################### *
       * ##  get the paritial match from pivoted_rhs_pattern  ## *
       * ##  to data_graph                                    ## *
       * ####################################################### */
      GUNDAM::Match<GraphPatternType, DataGraphType>
      pivoted_literal_pattern_to_data_graph_partial_match
        = pivoted_rhs_pattern_to_data_graph_match(
      pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match);
      
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              pivoted_literal_pattern,   data_graph, 
                              pivoted_literal_pattern_to_data_graph_partial_match,
                              pivoted_literal_pattern_candidate_set,
                              pivoted_literal_pattern_prune_callback,
                              pivoted_literal_pattern_match_callback,
                              time_limit_per_supp);
    }
    else {
      // the pivoted_rhs_gar is the same as pivoted_literal_gar
      // directly match literal_gar to data_graph
      if (!pivoted_literal_pattern_same_as_gar_pattern) {
        /* ####################################################### *
         * ##  get the paritial match from pivoted_rhs_pattern  ## *
         * ##  to data_graph                                    ## *
         * ####################################################### */
        GUNDAM::Match<GraphPatternType, DataGraphType>
                  gar_pattern_to_data_graph_partial_match
        = pivoted_rhs_pattern_to_data_graph_match(
                  gar_pattern_to_pivoted_rhs_pattern_partial_match);

        GUNDAM::MatchUsingMatch<match_semantics,
                                GUNDAM::MatchAlgorithm::kDagDp,
                                GUNDAM::MergeNecConfig::kNotMerge>(
                                gar_pattern,   data_graph, 
                                gar_pattern_to_data_graph_partial_match,
                                gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                                gar_pattern_prune_callback,
                                gar_pattern_match_callback,
                                time_limit_per_supp);
      }
      else {
        // pivoted_rhs_gar_pattern is the same as gar_pattern
        has_satisfy_x = true;
      }
    }
    assert(!has_satisfy_x_not_y);
    assert(!has_satisfy_x_and_y);
    if (!has_satisfy_x) {
      // continue matching
      return true;
    }
    // satisfy x
    // reset parameter
    has_satisfy_x = false;
    if (satisify_y_literals) {
      // satisfy x and y
      assert(!has_satisfy_x_not_y);
      has_satisfy_x_and_y = true;
      // continue matching
      return true;
    }
    // satisfy x but not y
    has_satisfy_x_and_y = false;
    has_satisfy_x_not_y =  true;
    // no longer needs matching
    return false;
  };

  /* #################################################### *
   * ##  callbacks from pivoted_pattern to data graph  ## *
   * #################################################### */
  std::function<bool(const MatchType&)>
    pivoted_pattern_prune_callback
    = [](const MatchType& 
    pivoted_pattern_to_data_graph_match) -> bool {
    // prune nothing
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_pattern_match_callback
    = [&](const MatchType& 
    pivoted_pattern_to_data_graph_match) -> bool {
    /* ######################################## *
     * ##  get a match for pivot vertex set  ## *
     * ######################################## */
    has_satisfy_x       = false;
    has_satisfy_x_not_y = false;
    has_satisfy_x_and_y = false;
    
    /* ####################################################### *
     * ##  get the paritial match from pivoted_rhs_pattern  ## *
     * ##  to data_graph                                    ## *
     * ####################################################### */
    GUNDAM::Match<GraphPatternType, DataGraphType>
    pivoted_rhs_pattern_to_data_graph_partial_match
      = pivoted_pattern_to_data_graph_match(
    pivoted_rhs_pattern_to_pivoted_pattern_partial_match);

    if (!pivoted_pattern_same_as_pivoted_rhs_pattern) {
      // pivoted_rhs_gar is not the same as literal_gar
      // first match the rhs literal, then match
      // the literal_gar
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              pivoted_rhs_pattern,   data_graph, 
                              pivoted_rhs_pattern_to_data_graph_partial_match,
                              pivoted_rhs_pattern_candidate_set,
                              pivoted_rhs_pattern_prune_callback,
                              pivoted_rhs_pattern_match_callback,
                              time_limit_per_supp);
      assert(!has_satisfy_x);
      assert(!has_satisfy_x_and_y 
          || !has_satisfy_x_not_y); // cannot be both true
      has_satisfy_x = has_satisfy_x_and_y 
                   || has_satisfy_x_not_y;
      if (!has_satisfy_x) {
        return true;
      }
      x_supp++;
      const bool kSatisfyXSuppCallbackRet
                = satisfy_x_supp_callback(pivoted_pattern_to_data_graph_match);
      if (has_satisfy_x_not_y) {
        return kSatisfyXSuppCallbackRet;
      }
      assert(has_satisfy_x_and_y);
      xy_supp++;
      const bool kSatisfyXYSuppCallbackRet
                = satisfy_xy_supp_callback(pivoted_pattern_to_data_graph_match);
      return kSatisfyXSuppCallbackRet
         && kSatisfyXYSuppCallbackRet;
    }
    /* ########################################################### *
     * ##  pivoted_pattern is the same as pivoted_rhs_pattern,  ## *
     * ##  directly whether the match satisfy y literals here   ## *
     * ########################################################### */
    assert(pivoted_rhs_pattern_to_data_graph_partial_match.size()
            == pivoted_pattern_to_data_graph_match.size());

    /* ########################################################### *
     * ##  get the paritial match from pivoted_literal_pattern  ## *
     * ##  to data_graph                                        ## *
     * ########################################################### */
    GUNDAM::Match<GraphPatternType, DataGraphType>
    pivoted_literal_pattern_to_data_graph_partial_match
          = pivoted_pattern_to_data_graph_match(
    pivoted_literal_pattern_to_pivoted_pattern_partial_match);
    if (!pivoted_rhs_gar_same_as_pivoted_literal_gar) {
      assert( GUNDAM::SamePattern(pivoted_rhs_pattern, pivoted_pattern));
      assert(!        SameGar    (pivoted_rhs_gar,     pivoted_literal_gar));
      assert(!GUNDAM::SamePattern(pivoted_rhs_pattern, pivoted_literal_pattern));
      
      assert(!has_satisfy_x);
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              pivoted_literal_pattern,   data_graph, 
                              pivoted_literal_pattern_to_data_graph_partial_match,
                              pivoted_literal_pattern_candidate_set,
                              pivoted_literal_pattern_prune_callback,
                              pivoted_literal_pattern_match_callback,
                              time_limit_per_supp);
    }
    else {
      // pivoted_rhs_pattern is the same as pivoted_literal_pattern
      // directly match to pivoted_literal_pattern
      /* ############################################### *
       * ##  get the paritial match from gar_pattern  ## *
       * ##  to data_graph                            ## *
       * ############################################### */
      assert(        SameGar    (pivoted_rhs_gar, pivoted_literal_gar));
      assert(GUNDAM::SamePattern(pivoted_pattern, pivoted_literal_pattern));
      /* ######################################################## *
       * ##  chech whether the match satisfy lhs literal here  ## *
       * ######################################################## */
      bool satisify_x_literals = true;
      for (const auto& x_literal_ptr : pivoted_literal_gar.x_literal_set()) {
        assert(x_literal_ptr->MappedBy(pivoted_literal_pattern_to_data_graph_partial_match));
        if (x_literal_ptr->Satisfy(pivoted_literal_pattern_to_data_graph_partial_match)) {
          continue;
        }
        satisify_x_literals = false;
        break;
      }
      if (!satisify_x_literals) {
        // continue matching
        assert(!has_satisfy_x);
        return true;
      }
      if (!pivoted_literal_pattern_same_as_gar_pattern) {
        GUNDAM::Match<GraphPatternType, DataGraphType>
              gar_pattern_to_data_graph_partial_match
        = pivoted_pattern_to_data_graph_match(
              gar_pattern_to_pivoted_pattern_partial_match);
        assert(!has_satisfy_x);
        GUNDAM::MatchUsingMatch<match_semantics,
                                GUNDAM::MatchAlgorithm::kDagDp,
                                GUNDAM::MergeNecConfig::kNotMerge>(
                                gar_pattern,   data_graph, 
                                gar_pattern_to_data_graph_partial_match,
                                gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                                gar_pattern_prune_callback,
                                gar_pattern_match_callback,
                                time_limit_per_supp);
      }
      else {
        assert(GUNDAM::SamePattern(pivoted_pattern, gar_pattern));
        // pivoted_pattern is the same as gar_pattern, no longer need to match
        assert(!has_satisfy_x);
        has_satisfy_x = true;
      }
    }
    assert(!has_satisfy_x_and_y);
    assert(!has_satisfy_x_not_y);

    if (!has_satisfy_x) {
      return true;
    }
    x_supp++;

    GUNDAM::Match<GraphPatternType, DataGraphType>
          gar_pattern_to_data_graph_partial_match
    = pivoted_pattern_to_data_graph_match(
          gar_pattern_to_pivoted_pattern_partial_match);

    const bool kSatisfyXSuppCallbackRet
              = satisfy_x_supp_callback(gar_pattern_to_data_graph_partial_match);

    /* ######################################################## *
    * ##  chech whether the match satisfy rhs literal here  ## *
    * ######################################################## */
    bool satisify_y_literals = true;
    for (const auto& y_literal_ptr : pivoted_rhs_gar.y_literal_set()) {
      assert(y_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_partial_match));
      if (y_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_partial_match)) {
        continue;
      }
      satisify_y_literals = false;
      break;
    }
    if (!satisify_y_literals) {
      return kSatisfyXSuppCallbackRet;
    }
    // satisfy y literals
    xy_supp++;
    // continue matching
    const bool kSatisfyXYSuppCallbackRet
              = satisfy_xy_supp_callback(gar_pattern_to_data_graph_partial_match);
    return kSatisfyXSuppCallbackRet
       && kSatisfyXYSuppCallbackRet;
  };

  MatchType pivoted_pattern_partial_match;
  for (auto pivoted_pattern_vertex_it = pivoted_pattern.VertexBegin();
           !pivoted_pattern_vertex_it.IsDone();
            pivoted_pattern_vertex_it++) {
    auto target_handle = gar_pattern_to_data_graph_partial_match.MapTo(pivoted_pattern_vertex_it->id());
    if (!target_handle) {
      continue;
    }
    pivoted_pattern_partial_match.AddMap(pivoted_pattern_vertex_it,
                                         target_handle);
  }

  // first match pivoted gar_pattern, then match rhs gar_pattern
  GUNDAM::MatchUsingMatch<match_semantics,
                          GUNDAM::MatchAlgorithm::kDagDp,
                          GUNDAM::MergeNecConfig::kNotMerge>(
                          pivoted_pattern, data_graph,
                          pivoted_pattern_partial_match,
                          pivoted_pattern_candidate_set,
                          pivoted_pattern_prune_callback,
                          pivoted_pattern_match_callback,
                          time_limit);
  assert(x_supp >= xy_supp);
  return std::pair(x_supp, xy_supp);
}

// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
                                  DataGraphType&  data_graph,
           const GUNDAM::Match<GraphPatternType,
                                  DataGraphType>& pattern_to_data_graph_partial_match,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& pattern_candidate_set,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  return GarSupp<match_semantics>
                (gar, data_graph, 
                 pattern_to_data_graph_partial_match,
                 pattern_candidate_set,
                 pattern_candidate_set,
                 satisfy_x_supp_callback,
                 satisfy_xy_supp_callback,
                 time_limit,
                 time_limit_per_supp,
                 pivoted_vertex_set);
}

// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
                                  DataGraphType&  data_graph,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& pattern_candidate_set,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  const GUNDAM::Match<GraphPatternType,
                         DataGraphType> match_state;
  return GarSupp<match_semantics>
                (gar, data_graph, 
                 match_state,
                 pattern_candidate_set,
                 pattern_candidate_set,
                 satisfy_x_supp_callback,
                 satisfy_xy_supp_callback,
                 time_limit,
                 time_limit_per_supp,
                 pivoted_vertex_set);
}

// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
                                  DataGraphType& data_graph,
           const GUNDAM::Match<GraphPatternType,
                                  DataGraphType>& pattern_to_data_graph_partial_match,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& pattern_candidate_set,
            int64_t  x_supp_limit = -1,
            int64_t xy_supp_limit = -1,
            double     time_limit = -1.0,
            double     time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()){

  using MatchType = GUNDAM::Match<GraphPatternType,
                                     DataGraphType>;

  assert(x_supp_limit == -1 
     || xy_supp_limit == -1 
     || (x_supp_limit >= xy_supp_limit));

  int64_t x_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, DataGraphType>&)>
    satisfy_x_supp_callback 
        = [&x_supp_counter,
           &x_supp_limit](const MatchType& match) {
    if (x_supp_limit == -1){
      // x_supp_limit is not specified 
      // continue matching
      return true;
    }
    assert(x_supp_counter < x_supp_limit);
    x_supp_counter++;
    if (x_supp_counter == x_supp_limit){
      // has reached x_supp_limit
      // stop matching
      return false;
    }
    return true;
  };

  int64_t xy_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, DataGraphType>&)> 
    satisfy_xy_supp_callback 
        = [&xy_supp_counter,
           &xy_supp_limit](const MatchType& match) {
    if (xy_supp_limit == -1) {
      // xy_supp_limit is not specified 
      // continue matching
      return true;
    }
    assert(xy_supp_counter < xy_supp_limit);
    xy_supp_counter++;
    if (xy_supp_counter == xy_supp_limit) {
      // has reached xy_supp_limit
      // stop matching
      return false;
    }
    return true;
  };
  
  auto [x_supp, xy_supp] = GarSupp<match_semantics>(
                                   gar, data_graph,
                             pattern_to_data_graph_partial_match,
                             pattern_candidate_set,
                             satisfy_x_supp_callback,
                            satisfy_xy_supp_callback,
                                  time_limit,
                                  time_limit_per_supp,
                                  pivoted_vertex_set);

  assert((( x_supp_limit == -1) || ( x_supp ==  x_supp_counter) && ( x_supp <=  x_supp_limit))
      && ((xy_supp_limit == -1) || (xy_supp == xy_supp_counter) && (xy_supp <= xy_supp_limit)));

  return std::pair(x_supp, xy_supp);
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& pattern_candidate_set,
     std::function<bool(const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type, 
                                       typename GUNDAM::VertexHandle<   DataGraphType>::type>&)> satisfy_x_supp_callback,
     std::function<bool(const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type, 
                                       typename GUNDAM::VertexHandle<   DataGraphType>::type>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()){

  const GUNDAM::Match<GraphPatternType,
                         DataGraphType> pattern_to_data_graph_partial_match;

  return GarSupp<match_semantics>(
                 gar, data_graph,
                 pattern_to_data_graph_partial_match,
                 pattern_candidate_set,
                 satisfy_x_supp_callback,
                satisfy_xy_supp_callback,
                              time_limit,
                              time_limit_per_supp,
                           pivoted_vertex_set );
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
  const GUNDAM::Match<GraphPatternType,
                         DataGraphType>& pattern_to_data_graph_partial_match,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()){
  
  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using CandidateSetType = std::map<GraphPatternVertexHandle, 
                        std::vector<   DataGraphVertexHandle>>;

  CandidateSetType pattern_candidate_set;
  if (!GUNDAM::_dp_iso::InitCandidateSet<match_semantics>(
                gar.pattern(),
                data_graph,
                pattern_candidate_set)) {
    return std::pair(0, 0);
  }
  if (!GUNDAM::_dp_iso::RefineCandidateSet(
                gar.pattern(),
                data_graph,
                pattern_candidate_set)) {
    return std::pair(0, 0);
  }

  return GarSupp<match_semantics>(
                 gar, data_graph,
                 pattern_to_data_graph_partial_match,
                 pattern_candidate_set,
                 satisfy_x_supp_callback,
                satisfy_xy_supp_callback,
                              time_limit,
                              time_limit_per_supp,
                           pivoted_vertex_set );
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()){

  const GUNDAM::Match<GraphPatternType,
                         DataGraphType> pattern_to_data_graph_partial_match;

  return GarSupp<match_semantics>(
                 gar, data_graph,
                 pattern_to_data_graph_partial_match,
                 satisfy_x_supp_callback,
                satisfy_xy_supp_callback,
                              time_limit,
                              time_limit_per_supp,
                           pivoted_vertex_set );
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& pattern_candidate_set,
            int64_t  x_supp_limit = -1,
            int64_t xy_supp_limit = -1,
            double     time_limit = -1.0,
            double     time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  GUNDAM::Match<GraphPatternType,
                   DataGraphType> pattern_to_data_graph_partial_match;
  return GarSupp<match_semantics>(
                 gar, data_graph, pattern_to_data_graph_partial_match, 
                                  pattern_candidate_set, 
                                           x_supp_limit, 
                                          xy_supp_limit,
                                             time_limit,
                                             time_limit_per_supp,
                                          pivoted_vertex_set);
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
            int64_t  x_supp_limit = -1,
            int64_t xy_supp_limit = -1,
            double     time_limit = -1.0,
            double     time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using MatchType = GUNDAM::Match<GraphPatternType, 
                                     DataGraphType>;

  using CandidateSetContainer = std::map<GraphPatternVertexHandle, 
                                std::vector<DataGraphVertexHandle>>;

  auto& pattern  = gar.pattern();

  CandidateSetContainer pattern_candidate_set;
  if (!GUNDAM::_dp_iso::InitCandidateSet<match_semantics>(
                pattern ,
                data_graph,
                pattern_candidate_set)) {
    return std::make_pair(0, 0);
  }
  if (!GUNDAM::_dp_iso::RefineCandidateSet(
                pattern , 
                data_graph, 
                pattern_candidate_set)) {
    return std::make_pair(0, 0);
  }
  return GarSupp<match_semantics>(
                 gar, data_graph, pattern_candidate_set, 
                 x_supp_limit, 
                xy_supp_limit,
                   time_limit,
                   time_limit_per_supp,
                pivoted_vertex_set );
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
inline std::pair<typename DataGraphType::VertexCounterType,
                 typename DataGraphType::VertexCounterType> 
  IncrementalGarSupp(
     GraphAssociationRule<GraphPatternType,
                             DataGraphType>& gar,
                             DataGraphType&  data_graph,
    std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type> &delta_target_graph,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
        double     time_limit = -1.0,
        double     time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;
  using CandidateSet = std::map<GraphPatternVertexHandle, 
                       std::vector<DataGraphVertexHandle>>;
  CandidateSet candidate_set;
  if (!GUNDAM::_dp_iso_using_match::InitCandidateSet<match_semantics>(
           gar.pattern(), 
              data_graph, 
           candidate_set)) {
    return std::pair(0, 0);
  }

  if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(
           gar.pattern(), 
              data_graph, 
           candidate_set)) {
    return std::pair(0, 0);
  }

  std::pair<typename DataGraphType::VertexCounterType,
            typename DataGraphType::VertexCounterType> gar_match_result;

  std::sort(delta_target_graph.begin(), 
            delta_target_graph.end());

  std::vector<GraphPatternVertexHandle> has_delta_target_graph_pattern_vertex;

  for (auto &[query_handle, target_list] : candidate_set) {
    bool found_new_vertex = false;
    for (auto &target_handle : target_list) {
      if (std::binary_search(delta_target_graph.begin(),
                             delta_target_graph. end (),
                                   target_handle)) {
        found_new_vertex = true;
        break;
      }
    }
    if (!found_new_vertex) {
      continue;
    }
    has_delta_target_graph_pattern_vertex.emplace_back(query_handle);
  }
  std::sort(has_delta_target_graph_pattern_vertex.begin(),
            has_delta_target_graph_pattern_vertex.end());
  int total_mask = (1 << (has_delta_target_graph_pattern_vertex.size()));
  for (int mask = 1; mask < total_mask; mask++) {
    std::vector<GraphPatternVertexHandle> this_mask_vertex;
    for (int bit_pos = 0;
             bit_pos < has_delta_target_graph_pattern_vertex.size(); 
             bit_pos++) {
      if (mask & (1 << bit_pos)) {
        this_mask_vertex.emplace_back(
            has_delta_target_graph_pattern_vertex[bit_pos]);
      }
    }
    CandidateSet copy_candidate_set{candidate_set};
    for (auto &[query_handle, target_list] : copy_candidate_set) {
      if (std::binary_search(this_mask_vertex.begin(), 
                             this_mask_vertex.end(),
                             query_handle)) {
        std::vector<DataGraphVertexHandle> this_vertex_target_list;
        for (auto &target_handle : target_list) {
          if (std::binary_search(delta_target_graph.begin(),
                                 delta_target_graph.end(), 
                                       target_handle)) {
            this_vertex_target_list.emplace_back(target_handle);
          }
        }
        std::swap(this_vertex_target_list, 
                              target_list);
        continue;
      }
      std::vector<DataGraphVertexHandle> this_vertex_target_list;
      for (auto &target_handle : target_list) {
        if (!std::binary_search(delta_target_graph.begin(),
                                delta_target_graph.end(), 
                                      target_handle)) {
          this_vertex_target_list.emplace_back(target_handle);
        }
      }
      std::swap(this_vertex_target_list, 
                            target_list);
    }

    auto partial_result = GarSupp<match_semantics>(
                                  gar, data_graph, copy_candidate_set, 
                                  satisfy_x_supp_callback,
                                 satisfy_xy_supp_callback,
                                    time_limit,
                                    time_limit_per_supp,
                                 pivoted_vertex_set );

    gar_match_result.first  += partial_result.first;
    gar_match_result.second += partial_result.second;
  }
  return gar_match_result;
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
inline std::pair<typename DataGraphType::VertexCounterType,
                 typename DataGraphType::VertexCounterType> 
  IncrementalGarSupp(
     GraphAssociationRule<GraphPatternType,
                             DataGraphType>& gar,
                             DataGraphType&  data_graph,
    std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type> &delta_target_graph,
        int64_t  x_supp_limit = -1,
        int64_t xy_supp_limit = -1,
        double     time_limit = -1.0,
        double     time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  using MatchType = GUNDAM::Match<GraphPatternType,
                                     DataGraphType>;

  assert(x_supp_limit == -1 
     || xy_supp_limit == -1 
     || (x_supp_limit >= xy_supp_limit));

  int64_t x_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, DataGraphType>&)>
    satisfy_x_supp_callback 
        = [&x_supp_counter,
           &x_supp_limit](const MatchType& match) {
    if (x_supp_limit == -1){
      // x_supp_limit is not specified 
      // continue matching
      return true;
    }
    assert(x_supp_counter < x_supp_limit);
    x_supp_counter++;
    if (x_supp_counter == x_supp_limit){
      // has reached x_supp_limit
      // stop matching
      return false;
    }
    return true;
  };

  int64_t xy_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, DataGraphType>&)> 
    satisfy_xy_supp_callback 
        = [&xy_supp_counter,
           &xy_supp_limit](const MatchType& match) {
    if (xy_supp_limit == -1) {
      // xy_supp_limit is not specified 
      // continue matching
      return true;
    }
    assert(xy_supp_counter < xy_supp_limit);
    xy_supp_counter++;
    if (xy_supp_counter == xy_supp_limit) {
      // has reached xy_supp_limit
      // stop matching
      return false;
    }
    return true;
  };
  
  return IncrementalGarSupp(gar, data_graph, delta_target_graph,
                            satisfy_x_supp_callback,
                            satisfy_xy_supp_callback,
                            time_limit,
                            time_limit_per_supp,
                            pivoted_vertex_set);
}

} // namespace gar

#endif // GAR_SUPP_H_
