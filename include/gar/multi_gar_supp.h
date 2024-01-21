#ifndef MULTI_GAR_SUPP_H_
#define MULTI_GAR_SUPP_H_

#include <map>

#include "gar/gar.h"
#include "gar/same_gar.h"
#include "gar/gar_supp.h"

#include "gundam/algorithm/match_using_match.h"

#include "gundam/match/match.h"

#include "gundam/type_getter/vertex_handle.h"
#include "gundam/type_getter/edge_handle.h"

#include "gundam/io/csvgraph.h"

namespace gar{

namespace _multi_gar_supp{

template <typename GraphPatternType,
          typename    DataGraphType>
class GarPatternGetter{
 private:
  using GarType = GraphAssociationRule<GraphPatternType,
                                          DataGraphType>;

 public:
  using size_type = typename std::vector<GarType>::size_type;

  GarPatternGetter(std::vector<GarType>*   gar_list_ptr,
             const std::vector<size_type>& vailable_gar_idx_set)
                            :gar_list_ptr_(gar_list_ptr),
                     vailable_gar_idx_set_(vailable_gar_idx_set){
    assert( this->gar_list_ptr_);
    assert(!this->vailable_gar_idx_set_.empty());
    return;
  }
  
  inline GraphPatternType& operator[]( size_type pos ){
    assert(pos >= 0 && pos < this->vailable_gar_idx_set_.size());
    auto vailable_gar_idx  = this->vailable_gar_idx_set_[pos];
    assert(vailable_gar_idx >= 0
        && vailable_gar_idx < this->gar_list_ptr_->size());
    return this->gar_list_ptr_[vailable_gar_idx].pattern();
  }
   
  inline const GraphPatternType& operator[]( size_type pos ) const {
    assert(pos >= 0 && pos < this->vailable_gar_idx_set_.size());
    auto vailable_gar_idx  = this->vailable_gar_idx_set_[pos];
    assert(vailable_gar_idx >= 0
        && vailable_gar_idx < this->gar_list_ptr_->size());
    return this->gar_list_ptr_[vailable_gar_idx].pattern();
  }

  inline size_type gar_idx_at( size_type pos ) const {
    assert(pos >= 0 && pos < this->vailable_gar_idx_set_.size());
    auto vailable_gar_idx  = this->vailable_gar_idx_set_[pos];
    return vailable_gar_idx;
  }

  inline size_type size() const {
    return vailable_gar_idx_set_.size();
  }

 private:
  std::vector<GarType>* gar_list_ptr_;

  std::vector<size_type> vailable_gar_idx_set_;
};

template <typename GraphPatternType,
          typename    DataGraphType>
class LiteralPatternGetter{
 private:
  using GarType = GraphAssociationRule<GraphPatternType,
                                          DataGraphType>;

 public:
  using size_type = typename std::vector<GarType>::size_type;

  GarPatternGetter(std::vector<std::pair<GarType, 
                                         std::vector<size_type>>>* literal_gar_list_ptr)
                                           : literal_gar_list_ptr_(literal_gar_list_ptr){
    assert( this->literal_gar_list_ptr_);
    assert(!this->vailable_gar_idx_set_.empty());
    return;
  }
  
  inline GraphPatternType& operator[]( size_type pos ){
    assert(pos >= 0 && pos < this->vailable_gar_idx_set_.size());
    auto vailable_gar_idx  = this->vailable_gar_idx_set_[pos];
    assert(vailable_gar_idx >= 0
        && vailable_gar_idx < this->gar_list_ptr_->size());
    return this->gar_list_ptr_[vailable_gar_idx].pattern();
  }
   
  inline const GraphPatternType& operator[]( size_type pos ) const {
    assert(pos >= 0 && pos < this->vailable_gar_idx_set_.size());
    auto vailable_gar_idx  = this->vailable_gar_idx_set_[pos];
    assert(vailable_gar_idx >= 0
        && vailable_gar_idx < this->gar_list_ptr_->size());
    return this->gar_list_ptr_[vailable_gar_idx].pattern();
  }

  inline size_type size() const {
    return vailable_gar_idx_set_.size();
  }

 private:
  std::vector<std::pair<GarType, std::vector<size_type>>>* literal_gar_list_ptr_;
};

};

template <typename GraphPatternType,
          typename    DataGraphType>
void MultiGarSupp(
        std::vector<GraphAssociationRule<GraphPatternType,
                                            DataGraphType>>& gar_list,
             DataGraphType& data_graph,
     std::function<bool(int,
                        const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type, 
                                       typename GUNDAM::VertexHandle<   DataGraphType>::type>&)> satisfy_x_supp_callback,
     std::function<bool(int,
                        const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type, 
                                       typename GUNDAM::VertexHandle<   DataGraphType>::type>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0){
  // gar_list should have the same rhs
  using MatchMap = std::map<GraphPatternVertexHandle, 
                               DataGraphVertexHandle>;

  using CandidateSetContainer = std::map<GraphPatternVertexHandle, 
                                std::vector<DataGraphVertexHandle>>;

  using GarType = GraphAssociationRule<GraphPatternType,
                                          DataGraphType>;

  assert(!gar_list.empty());

  GarType rhs_gar = _gar_supp::RhsGar(*gar_list.begin());
  #ifndef NDEBUG
  // to verify that all gar have the same rhs pattern
  for (const auto& gar : gar_list) {
    GarType current_rhs_gar = _gar_supp::RhsGar(gar);
    assert(SameGar(current_rhs_gar, rhs_gar));
  }
  #endif // NDEBUG

  constexpr int kLiteralGarIdx = 0,
                 kGarIdxSetIdx = 1,
             kPatternGetterIdx = 2;

  std::vector<std::tuple<GarType, 
                         std::vector<size_t>,
                         PatternGetterType>> literal_gar_tuple_list;
  literal_gar_tuple_list.reserve(gar_list.size());

  // get the literal graph for the pattern of each gar
  for (size_t gar_idx = 0; gar_idx < gar_list.size(); gar_idx++){
    GarType literal_gar = _gar_supp::LiteralGar(gar_list[gar_idx]);

    bool has_same_literal_gar = false;
    for (const auto& existed_literal_gar
                           : literal_gar_tuple_list) {
      if (SameGar(literal_gar, existed_literal_gar.first)) {
        // if there exist same literal gar, then add this idx to it 
        auto [ gar_idx_it,
               gar_idx_ret ] = existed_literal_gar.second.emplace_back(gar_idx);
        assert(gar_idx_ret);
        has_same_literal_gar = true;
        break;
      }
    }
    // does not contain the same literal pattern before
    if (!has_same_literal_gar) {
      literal_gar_tuple_list.emplace_back(std::move(literal_gar),
                                    std::vector<size_t>{gar_idx},
                                    PatternGetterType());
    }
  }
  // obtain a list of different literal gars and the set of gars
  // that corresponds to them
  literal_gar_tuple_list.shrink_to_fit();

  #ifndef NDEBUG
  // the same gar_idx should not appear two times
  // in the set of gar_idx
  std::set<size_t> gar_idx_set;
  for (const auto& literal_gar : literal_gar_tuple_list) {
    for (const auto& gar_idx : std::get<kGarIdxSetIdx>(literal_gar)){
      auto [ gar_idx_set_it,
             gar_idx_set_ret ] = gar_idx_set.emplace(gar_idx);
      // should added successfully
      assert(gar_idx_set_ret);
    }
  }
  // all gar should have been added
  assert(gar_idx_set.size() == gar_list.size());
  #endif // NDEBUG

  // construct a pattern getter for the set of gars that it contains
  for (auto& literal_gar_tuple : literal_gar_tuple_list) {
    const auto& gar_idx_set = std::get<kGarIdxSetIdx>(literal_gar_tuple);
    // construct the pattern getter based on the
    // gathered set of gar_idx that share this literal_gar
    PatternGetterType pattern_getter(&gar_list,
                                      gar_idx_set);
    std::get<kPatternGetterIdx>(literal_gar_tuple) = std::move(pattern_getter);
  }

  // mark whether each gar needs to be continue matching
  std::vector<bool> reach_end_condition(gar_list.size(), false);

  std::vector<size_t> x_supp_list(gar_list.size(), 0),
                     xy_supp_list(gar_list.size(), 0);

  auto rhs_pattern_match_callback
    = [](const MatchMap& match) {
    if (prune_rhs_pattern_not_satisfy_x_literal_set_callback(match)){
      // this match does not satisfy all literals
      // continue matching
      return true;
    }

    // mark whether each gar needs to be continue matching
    std::vector<bool> has_literal_pattern_match(gar_list.size(), false);

    std::vector<bool> match_set;

    // to mark whether has been completed
    std::vector<bool> has_complete_match(gar_list.size(), false);

    for (auto& literal_gar_tuple 
             : literal_gar_tuple_list) {
      auto& literal_gar        = std::get<   kLiteralGarIdx>(literal_gar_tuple);
      auto& gar_idx_set        = std::get<    kGarIdxSetIdx>(literal_gar_tuple);
      auto& gar_pattern_getter = std::get<kPatternGetterIdx>(literal_gar_tuple);
      
      assert(gar_idx_set.size()
          == gar_pattern_getter.size());

      #ifndef NDEBUG
      size_t match_gar_counter = 0;
      const size_t kTotalGarToMatch = gar_pattern_getter.size();
      #endif // NDEBUG

      auto prune_literal_pattern_not_satisfy_x_literal_set_callback
             = [&literal_gar,
                &has_complete_match](const MatchMap& match_state){
        for (const auto& x_literal_ptr 
           : literal_gar.x_literal_set()) {
          if (!x_literal_ptr->MappedBy(match_state)){
            // this literal is not considered in the current
            // match_state, this match cannot be pruned
            // move to the next literal
            continue;
          }
          if (!x_literal_ptr->Satisfy(match_state)) {
            // does not satisfy all literals
            // does not need further matching
            return true;
          }
        }
        // satisfy all literals
        // continue matching
        return false;
      };

      // match for each literal gar based on the 
      // current match from rhs gar to data graph
      auto literal_pattern_match_callback 
       = [&literal_pattern](const MatchMap& match) -> bool {
        // has got a match from literal gar to data graph
        // try to complete this match for the gars
        // that shares this literal gar
        auto literal_gar_idx_set = literal_gar.second;

        auto gar_prune_callback = [&has_complete_match](
                        int gar_getter_idx,
                        const MatchMap& match){
          assert(gar_getter_idx >= 0
              && gar_getter_idx < gar_pattern_getter.size());
          auto gar_idx = gar_pattern_getter.gar_idx_at(size);
          assert(gar_idx >= 0
              && gar_idx < has_complete_match.size());
          // needs to be pruned since the match_callback
          // is no longer needs to be called
          if (!has_complete_match[gar_idx]){
            return true;
          }
          return false;
        };

        auto gar_match_callback = [&has_complete_match
                                  #ifndef NDEBUG
                                  // parameter to be used in debug
                                  ,&match_gar_counter
                                  ,&kTotalGarToMatch
                                  #endif // NDEBUG
                                  ](
                        int gar_getter_idx,
                        const MatchMap& match) {
          #ifndef NDEBUG
          match_gar_counter++;
          assert(match_gar_counter <= kTotalGarToMatch);
          #endif // NDEBUG
          assert(gar_getter_idx >= 0
              && gar_getter_idx < gar_pattern_getter.size());
          auto gar_idx = gar_pattern_getter.gar_idx_at(size);
          assert(gar_idx >= 0
              && gar_idx < has_complete_match.size());
          // mark this pattern has been completed match
          has_complete_match[gar_idx] = true;
          // no longer need to match
          return false;
        };

        GUNDAM::MultiQueryDpiso(gar_pattern_getter, 
                                     data_graph,
                                gar_prune_callback,
                                gar_match_callback);
      };
      
      auto& literal_pattern = literal_gar.first.pattern();

      if (!SameAs(literal_gar,
                      rhs_gar)){
        GUNDAM::MatchUsingMatch(literal_pattern,   data_graph,
                                literal_pattern_to_data_graph_partial_match,
                                literal_pattern_candidate_set,
              prune_literal_pattern_not_satisfy_x_literal_set_callback,
                                        literal_pattern_match_callback,
                                        time_limit_per_supp);
        continue;
      }

      if () {
        
      }
      
      // is literal pattern is the same as rhs pattern
      // does not need to further matching
      GUNDAM::MultiQueryDpiso(gar_pattern_getter,   data_graph,
                              gar_pattern_getter_to_common_pattern_partial_match,
                              common_pattern_to_data_graph_match,
                              gar_pattern_candidate_set,
                              gar_prune_callback,
                              gar_match_callback,
                              time_limit_per_supp);
    }

    for (int gar_idx = 0; 
             gar_idx < has_complete_match.size(); 
             gar_idx++){
      if (!has_complete_match[gar_idx]) {
        // this support does not satisfy X for gar_idx'th pattern
        continue;
      }
      // this support satisfy X for gar_idx'th pattern
      x_supp_list[gar_idx]++;
      GUNDAM::Match<GraphPatternType, 
                       DataGraphType>
            pattern_to_data_graph_partial_match
      = rhs_pattern_to_data_graph_match(
            pattern_to_rhs_pattern_partial_match[gar_idx]);
      assert(pattern_to_data_graph_partial_match.size()
          == pattern_to_rhs_pattern_partial_match[gar_idx].size());
      assert(pattern_to_data_graph_partial_match.size()
          == rhs_pattern_to_data_graph_match.size());

      MatchMap pattern_to_data_graph_match;
      for (auto map_it = pattern_to_data_graph_partial_match.MapBegin();
               !map_it.IsDone();
                map_it++){
        pattern_to_data_graph_match.emplace(map_it->src_handle(),
                                            map_it->dst_handle());
      }
      assert(pattern_to_data_graph_match.size()
                                == match.size());
      
      const bool kSatisfyXSuppCallbackRet
                = satisfy_x_supp_callback(gar_idx, pattern_to_data_graph_match);

      bool satisify_y_literals = true;
      for (const auto& y_literal_ptr 
             : rhs_gar.y_literal_set()){
        assert(y_literal_ptr->MappedBy(match));
        if (y_literal_ptr->Satisfy(match)){
          continue;
        }
        satisify_y_literals = false;
        break;
      }
      if (satisify_y_literals) {
        xy_supp_list[gar_idx]++;
        const bool kSatisfyXYSuppCallbackRet
                  = satisfy_xy_supp_callback(gar_idx, pattern_to_data_graph_match);
        if (!kSatisfyXSuppCallbackRet
        || !kSatisfyXYSuppCallbackRet){
          // reached ending condition, end matching
          assert(0 <= gar_idx && gar_idx < reach_end_condition.size());
          reach_end_condition[gar_idx] = true;
          continue;
        }
        assert(kSatisfyXSuppCallbackRet
            && kSatisfyXYSuppCallbackRet);
        // does not reach ending condition, continue matching
        return true;
      }
      if (!kSatisfyXSuppCallbackRet){
        // reached ending condition, end matching
        reach_end_condition[gar_idx] = true;
        continue;
      }
    }

    bool satisify_y_literals = true;
    for (const auto& y_literal_ptr 
            : rhs_gar.y_literal_set()){
      assert(y_literal_ptr->MappedBy(match));
      if (y_literal_ptr->Satisfy(match)){
        continue;
      }
      satisify_y_literals = false;
      break;
    }
    if (satisify_y_literals){
      xy_supp++;
      const bool kSatisfyXYSuppCallbackRet
                = satisfy_xy_supp_callback(pattern_to_data_graph_match);
      if (!kSatisfyXSuppCallbackRet
       || !kSatisfyXYSuppCallbackRet){
        // reached ending condition, end matching
        reach_end_condition[gar_idx] = true;
        continue;
      }
      assert(kSatisfyXSuppCallbackRet
          && kSatisfyXYSuppCallbackRet);
      // does not reach ending condition, continue matching
      continue;
    }
    if (!kSatisfyXSuppCallbackRet){
      // reached ending condition, end matching
      reach_end_condition[gar_idx] = true;
      continue;
    }
    // does not reach ending condition, continue matching
    return true;
  };

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    rhs_pattern_to_pattern_match
   (rhs_pattern,   pattern , "same_id_map");

  GUNDAM::Match<GraphPatternType, DataGraphType> 
    rhs_pattern_to_data_graph_partial_match
      = pattern_to_data_graph_partial_match(
    rhs_pattern_to_pattern_match);

  GUNDAM::MatchUsingMatch(rhs_pattern,   data_graph,
                          rhs_pattern_to_data_graph_partial_match,
                          rhs_pattern_candidate_set,
                    prune_rhs_pattern_not_satisfy_x_literal_set_callback,
                                              rhs_pattern_match_callback,
                          time_limit);
  return std::make_pair(x_supp, xy_supp);
}

} // namespace gar

#endif // MULTI_GAR_SUPP_H_
