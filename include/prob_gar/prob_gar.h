#ifndef _PROB_GAR_PROB_GAR_H
#define _PROB_GAR_PROB_GAR_H

#include "gar/gar.h"

namespace prob_gar {

template <typename   Pattern, 
          typename DataGraph>
class ProbGraphAssociationRule 
 : public gar::GraphAssociationRule<Pattern, DataGraph> {
 private:
  using GraphAssociationRuleType 
 = gar::GraphAssociationRule<Pattern, DataGraph>;

  using PatternVertexHandleType = typename GUNDAM::VertexHandle<Pattern>::type;

 public:
  ProbGraphAssociationRule()
    : GraphAssociationRuleType(),
      prob_(1.0){
    this->pivot_vertex_set_.clear();
    for (auto vertex_it = this->pattern().VertexBegin();
             !vertex_it.IsDone();
              vertex_it++) {
      this->pivot_vertex_set_.emplace_back(vertex_it);
    }
    return;
  }

  ProbGraphAssociationRule(const Pattern& pattern)
    : GraphAssociationRuleType(pattern),
      prob_(1.0){
    this->pivot_vertex_set_.clear();
    for (auto vertex_it = this->pattern().VertexBegin();
             !vertex_it.IsDone();
              vertex_it++) {
      this->pivot_vertex_set_.emplace_back(vertex_it);
    }
    return;
  }

  ProbGraphAssociationRule(const Pattern& pattern,
                                    float prob)
    : GraphAssociationRuleType(pattern),
      prob_(prob){
    this->pivot_vertex_set_.clear();
    for (auto vertex_it = this->pattern().VertexBegin();
             !vertex_it.IsDone();
              vertex_it++) {
      this->pivot_vertex_set_.emplace_back(vertex_it);
    }
    return;
  }

  ProbGraphAssociationRule(
    const GraphAssociationRuleType& gar)
         :GraphAssociationRuleType( gar),
          prob_(1.0){
    this->pivot_vertex_set_.clear();
    for (auto vertex_it = this->pattern().VertexBegin();
             !vertex_it.IsDone();
              vertex_it++) {
      this->pivot_vertex_set_.emplace_back(vertex_it);
    }
    return;
  }

  ProbGraphAssociationRule(
    const GraphAssociationRuleType& gar,
                              float prob)
         :GraphAssociationRuleType( gar),
          prob_(prob){
    this->pivot_vertex_set_.clear();
    for (auto vertex_it = this->pattern().VertexBegin();
             !vertex_it.IsDone();
              vertex_it++) {
      this->pivot_vertex_set_.emplace_back(vertex_it);
    }
    return;
  }

  inline void Reset() {
    GraphAssociationRuleType::Reset();
    this->prob_ = 1.0;
    return;
  }

  inline const auto& prob() const {
    return this->prob_;
  }

  inline void set_prob(float prob) {
    assert(prob >= 0 && prob <= 1.0);
    this->prob_ = prob;
    return;
  }

 private:
  float prob_;

  std::vector<PatternVertexHandleType> pivot_vertex_set_;
};

}; // namespace prob_gar

#endif // _PROB_GAR_PROB_GAR_H