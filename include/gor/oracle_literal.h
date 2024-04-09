#ifndef EXAMPLES_ANALYTICAL_ORACLE_LITERAL_H_
#define EXAMPLES_ANALYTICAL_ORACLE_LITERAL_H_

#include "gar/literal.h"
#include "gar/literal_info.h"


// 记录oracle literal
// 根据gar的ML literal模板仿写的oracle模板
namespace gor {

template <class Pattern, class DataGraph>
class OracleLiteral : public gar::Literal<Pattern, DataGraph> {
 private:
  using BaseLiteralType = gar::Literal<Pattern, DataGraph>;

  using LiteralInfoType = typename BaseLiteralType::LiteralInfoType;

  using LiteralStandAloneInfoType =
      typename BaseLiteralType::LiteralStandAloneInfoType;

  using PatternVertexHandle = typename BaseLiteralType::PatternVertexHandle;
  using PatternVertexIDType = typename BaseLiteralType::PatternVertexIDType;

  using DataGraphVertexIDType = typename BaseLiteralType::DataGraphVertexIDType;
  using DataGraphVertexHandle = typename BaseLiteralType::DataGraphVertexHandle;
  
  using DataGraphEdgeIDType = typename BaseLiteralType::DataGraphEdgeIDType;
  using DataGraphEdgeHandle = typename BaseLiteralType::DataGraphEdgeHandle;

  using DataGraphEdgeLabelType = typename BaseLiteralType::DataGraphEdgeLabelType;

  //一阶函数和二阶函数，定义说明：
  // 定义 Oracle 函数的类型
  using UnaryOracleFunction = std::function<bool(const DataGraphVertexHandle&)>;
  using BinaryOracleFunction = std::function<bool(const DataGraphVertexHandle&, const DataGraphVertexHandle&)>;


 public:
  //   一元函数
  OracleLiteral(Pattern &pattern, const PatternVertexIDType &src_id,
                UnaryOracleFunction unary_function)
      : src_handle_{pattern.FindVertex(src_id)},
        unary_oracle_function_{unary_function},
        is_unary_{true} {
    assert(src_handle_);
  }

  // 二元函数
  OracleLiteral(Pattern &pattern, const PatternVertexIDType &src_id,
                const PatternVertexIDType &dst_id,
                BinaryOracleFunction binary_function)
      : src_handle_{pattern.FindVertex(src_id)},
        dst_handle_{pattern.FindVertex(dst_id)},
        binary_oracle_function_{binary_function},
        is_unary_{false} {
    assert(src_handle_);
    assert(dst_handle_);
  }

  ~OracleLiteral() {
    return;
  }


  virtual bool Satisfy(const std::map<PatternVertexHandle, DataGraphVertexHandle> &match_result) const override {
    if (is_unary_) {
      // 一元 Oracle 函数的处理
      auto it = match_result.find(src_handle_);
      if (it == match_result.end()) return false;
      return unary_oracle_function_(it->second);
    } else {
      // 二元 Oracle 函数的处理
      auto it_src = match_result.find(src_handle_);
      auto it_dst = match_result.find(dst_handle_);
      if (it_src == match_result.end() || it_dst == match_result.end()) return false;
      return binary_oracle_function_(it_src->second, it_dst->second);
    }
  }


  virtual bool Update(
      const std::map<PatternVertexHandle, DataGraphVertexHandle> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<typename BaseLiteralType::DataGraphEdgeIDType> &edge_id_gen,
      std::set<typename BaseLiteralType::DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<typename BaseLiteralType::DataGraphEdgeIDType> *diff_edge_set = nullptr) const override {
    // OracleLiteral may not directly modify the data graph but might influence other operations
    // The update logic depends on the specific Oracle function
    return false; // Returning false as default, meaning no direct update is performed
  }


  // Other necessary member functions...

 private:
  PatternVertexHandle src_handle_, dst_handle_;
  bool is_unary_;  // 标记这是一元还是二元 Oracle
  UnaryOracleFunction unary_oracle_function_;
  BinaryOracleFunction binary_oracle_function_;

};

}  // namespace gar

#endif  // EXAMPLES_ANALYTICAL_ORACLE_LITERAL_H_
