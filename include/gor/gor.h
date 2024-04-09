#ifndef GOR_GOR_H_
#define GOR_GOR_H_

#include "gar/gar.h"
#include "oracle_literal.h"

namespace gor {

// template <typename Pattern, 
//           typename DataGraph>
// class GraphOracleRule : public gar::GraphAssociationRule<Pattern,
//                                                        DataGraph> {
//  public:
//    using GraphAssociationRuleType
//   = gar::GraphAssociationRule<Pattern, DataGraph>;
//  public:
//   using GraphAssociationRuleType::GraphAssociationRuleType;
  
// };


// does not support additional literal now
// add by bowen 
template <typename Pattern, 
          typename DataGraph>
class GraphOracleRule : public gar::GraphAssociationRule<Pattern,
                                                       DataGraph> {
 public:
   using GraphAssociationRuleType
  = gar::GraphAssociationRule<Pattern, DataGraph>;
 public:
  using GraphAssociationRuleType::GraphAssociationRuleType;
  using DataGraphType = DataGraph;
  using PatternType = Pattern;
  using DataGraphVertexType   = typename DataGraph::VertexType;
  using UnaryPredicate = std::function<bool(const DataGraphVertexType&)>;
  using BinaryPredicate = std::function<bool(const DataGraphVertexType&, const DataGraphVertexType&)>; 
  // 添加新的谓词相关类型定义
  
  struct PredicateInfo {
    enum class Type {
      Unary,
      Binary,
      Aggregation,
      MachineLearning
    } type;
    
    UnaryPredicate unary_predicate; //一元谓词
    BinaryPredicate binary_predicate; //二元谓词
  };

    // 添加处理谓词的新方法
  void AddUnaryPredicate(const typename Pattern::VertexType::IDType& pattern_vertex_id, UnaryPredicate predicate) {
    predicate_infos_[pattern_vertex_id] = {PredicateInfo::Type::Unary, predicate, {}};
  }

  void AddBinaryPredicate(const typename Pattern::VertexType::IDType& pattern_vertex_id1, 
                          const typename Pattern::VertexType::IDType& pattern_vertex_id2, 
                          BinaryPredicate predicate) {
    binary_predicate_infos_[std::make_pair(pattern_vertex_id1, pattern_vertex_id2)] = predicate;
  }

  private:
    PatternType pattern_;
    // LiteralSetType x_literal_set_, y_literal_set_; //gar有相关定义了，只需要添加谓词结构应该就可以

    // 添加存储谓词信息的数据结构
    std::map<typename Pattern::VertexType::IDType, PredicateInfo> predicate_infos_;
    std::map<std::pair<typename Pattern::VertexType::IDType, typename Pattern::VertexType::IDType>, BinaryPredicate> binary_predicate_infos_;

};

}  // namespace gor

#endif