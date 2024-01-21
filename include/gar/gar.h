#ifndef GAR_H_
#define GAR_H_

#include "gundam/component/container2.h"
#include "gundam/component/generator.h"

#include "gundam/io/rapidcsv.h"
#include "gundam/algorithm/vf2.h"

#include "literal.h"
#include "ml_literal.h"

namespace gar {

template <class Pattern, class DataGraph>
class GraphAssociationRule {
 public:
  using PatternType = Pattern;
  using DataGraphType = DataGraph;

 public:
  using PatternVertexType   = typename Pattern::VertexType;
  using PatternVertexIDType = typename Pattern::VertexType::IDType;

  using DataGraphVertexType   = typename DataGraph::VertexType;
  using DataGraphVertexIDType = typename DataGraph::VertexType::IDType;
  
  using DataGraphEdgeType      = typename DataGraph::EdgeType;
  using DataGraphEdgeLabelType = typename DataGraphEdgeType::LabelType;
  using DataGraphVertexAttributeKeyType =
      typename DataGraph::VertexType::AttributeKeyType;

  using VertexAttributePair =
      std::pair<PatternVertexIDType, DataGraphVertexAttributeKeyType>;

  using LiteralSetType = GUNDAM::PointerVector<Literal<Pattern, DataGraph>>;

  void AddLiteral(LiteralSetType& literal_set,
            const LiteralInfo<Pattern, 
                            DataGraph>& literal_info){
    switch (literal_info.literal_type()) {
    case LiteralType::kVariableLiteral:
      literal_set.template Add<VariableLiteral<Pattern, 
                                             DataGraph>>(this->pattern_, literal_info);
      break;
    case LiteralType::kAttrValueLiteral:
      literal_set.template Add<AttributeLiteral<Pattern, 
                                              DataGraph>>(this->pattern_, literal_info);
      break;
    case LiteralType::kConstantLiteral:
      literal_set.template Add<ConstantLiteral<Pattern, 
                                             DataGraph>>(this->pattern_, literal_info);
      break;
    case LiteralType::kEdgeLiteral:
      literal_set.template Add<EdgeLiteral<Pattern, 
                                         DataGraph>>(this->pattern_, literal_info);
      break;
    case LiteralType::kMlLiteral:
      literal_set.template Add<MlLiteral<Pattern, 
                                       DataGraph>>(this->pattern_, literal_info);
      break;
    case LiteralType::kKeyLiteral:
      literal_set.template Add<KeyLiteral<Pattern, 
                                        DataGraph>>(this->pattern_, literal_info);
      break;
    default:
      assert(false);
      break;
    }
    return;
  }

 public:
  GraphAssociationRule(){};

  //template <class UPattern>
  //GraphAssociationRule(UPattern &&pattern)
  //    : pattern_(std::forward<UPattern>(pattern)){};

  GraphAssociationRule(const Pattern &pattern) 
                           : pattern_(pattern){
    return;
  };

  GraphAssociationRule &operator=(const GraphAssociationRule &b) {
    this->pattern_ = b.pattern_;
    this->x_literal_set_.Clear();
    this->y_literal_set_.Clear();
    this->x_literal_set_ = b.x_literal_set_;
    this->y_literal_set_ = b.y_literal_set_;
    for (auto &ptr : x_literal_set_) {
      ptr = ptr->Copy(pattern_);
    }
    for (auto &ptr : y_literal_set_) {
      ptr = ptr->Copy(pattern_);
    }
    return *this;
  }

  GraphAssociationRule(const GraphAssociationRule &b)
                     : pattern_(b.pattern_) {
    x_literal_set_ = b.x_literal_set_;
    y_literal_set_ = b.y_literal_set_;
    for (auto &ptr : x_literal_set_) {
      ptr = ptr->Copy(pattern_);
    }
    for (auto &ptr : y_literal_set_) {
      ptr = ptr->Copy(pattern_);
    }
    return;
  }

  GraphAssociationRule(GraphAssociationRule &&) = default;

  GraphAssociationRule &operator=(GraphAssociationRule &&) = default;

  ~GraphAssociationRule() {}

  void AddX(const LiteralInfo<Pattern, 
                              DataGraph>& literal_info) {
    this->AddLiteral(this->x_literal_set_,
                             literal_info);
    return;
  }

  void AddY(const LiteralInfo<Pattern, 
                              DataGraph>& literal_info) {
    this->AddLiteral(this->y_literal_set_,
                             literal_info);
    return;
  }

  template <class T, typename... Args>
  void AddX(Args... args) {
    x_literal_set_.template Add<T>(pattern_, args...);
    return;
  }

  template <class T, typename... Args>
  void AddY(Args... args) {
    y_literal_set_.template Add<T>(pattern_, args...);
    return;
  }

  void Reset() {
    x_literal_set_.Clear();
    y_literal_set_.Clear();
    pattern_.Clear();
  }

  // removed by wenzhi, literals would no-longer work 
  // after the pattern is set
  // template <class UPattern>
  // void SetPattern(UPattern &&pattern) {
  //   pattern_ = std::forward<UPattern>(pattern);
  // }

  // marked as should be removed by wenzhi, pattern cannot be modified
  // in order to avoid the literals no longer work
  // however, the program cannot be compiled once it is removed
  Pattern &pattern() { return pattern_; }
  const Pattern &pattern() const { return pattern_; }

  LiteralSetType &x_literal_set() { return x_literal_set_; }
  const LiteralSetType &x_literal_set() const { return x_literal_set_; }

  LiteralSetType &y_literal_set() { return y_literal_set_; }
  const LiteralSetType &y_literal_set() const { return y_literal_set_; }

  inline void ClearXLiteralSet() {
    this->x_literal_set_.Clear();
    return;
  }

  inline void ClearYLiteralSet() {
    this->y_literal_set_.Clear();
    return;
  }

 private:
  PatternType pattern_;
  LiteralSetType x_literal_set_, 
                 y_literal_set_;
};

}  // namespace gar

#endif