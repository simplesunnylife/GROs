#ifndef EXAMPLES_ANALYTICAL_LITERAL_H_
#define EXAMPLES_ANALYTICAL_LITERAL_H_

#include <cstdint>
#include <map>
#include <set>

// #include "gar/literal_info.h"
#include "gar/literal_stand_alone_info.h"
#include "gundam/component/generator.h"
#include "gundam/data_type/datatype.h"
#include "gundam/match/match.h"

#include "gundam/type_getter/edge_attribute_handle.h"
#include "gundam/type_getter/edge_handle.h"
#include "gundam/type_getter/edge_label.h"
#include "gundam/type_getter/edge_id.h"

#include "gundam/type_getter/vertex_attribute_handle.h"
#include "gundam/type_getter/vertex_handle.h"
#include "gundam/type_getter/vertex_label.h"
#include "gundam/type_getter/vertex_id.h"

namespace gar {


template <class Pattern, class DataGraph>
class LiteralInfo;

static constexpr char CsvSeparator = ',';

template <typename T>
std::string ProtectedSeparatorVal(T &&val) {
  std::stringstream s_stream;
  s_stream << val;
  std::string stream_str = s_stream.str();
  if (stream_str.find(CsvSeparator) != stream_str.npos &&
      !(stream_str[0] == '"' && stream_str.back() == '"')) {
    return (std::string) "\"" + s_stream.str() + (std::string) "\"";
  }
  return s_stream.str();
}
template <class Pattern, class DataGraph>
class Literal {
 public:
  using LiteralInfoType 
      = LiteralInfo<Pattern, DataGraph>;

  using LiteralStandAloneInfoType 
      = LiteralStandAloneInfo<Pattern, DataGraph>;

  using PatternVertexIDType = typename GUNDAM::VertexID    <Pattern>::type;
  using PatternVertexHandle = typename GUNDAM::VertexHandle<Pattern>::type;

  using DataGraphVertexIDType = typename GUNDAM::VertexID    <DataGraph>::type;
  using DataGraphVertexHandle = typename GUNDAM::VertexHandle<DataGraph>::type;

  using DataGraphEdgeIDType    = typename GUNDAM::EdgeID    <DataGraph>::type;
  using DataGraphEdgeHandle    = typename GUNDAM::EdgeHandle<DataGraph>::type;
  using DataGraphEdgeLabelType = typename GUNDAM::EdgeLabel <DataGraph>::type;

  using DataGraphVertexAttributeKeyType =
      typename DataGraph::VertexType::AttributeKeyType;

  virtual ~Literal() {}
  virtual bool Satisfy(
      const std::map<PatternVertexHandle, 
                   DataGraphVertexHandle> &) const = 0;

  virtual bool Satisfy(const GUNDAM::Match<Pattern, DataGraph>&) const = 0;

  virtual bool Update(
      const std::map<PatternVertexHandle, 
                   DataGraphVertexHandle> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<DataGraphEdgeIDType> &edge_id_gen,
      std::set<DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<DataGraphEdgeIDType>   *diff_edge_set        = nullptr) const = 0;

  virtual bool Update(
      const GUNDAM::Match<Pattern, DataGraph> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<DataGraphEdgeIDType> &edge_id_gen,
      std::set<DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<DataGraphEdgeIDType>   *diff_edge_set        = nullptr) const = 0;

  virtual Literal<Pattern, DataGraph> *Copy(Pattern &pattern) const = 0;

  virtual void Print(std::ostream &out,
                     std::string gar_name = std::string()) const = 0;

  virtual void CalPivot(std::set<PatternVertexHandle> &pivot) const = 0;

  virtual void CalPivot(std::vector<PatternVertexHandle> &pivot) const = 0;

  virtual LiteralInfoType info() const = 0;
  virtual LiteralStandAloneInfoType stand_alone_info() const = 0;
  virtual enum LiteralType type() const = 0;
  virtual bool MappedBy(
      const std::map<PatternVertexHandle, 
                   DataGraphVertexHandle> &) const = 0;

  virtual bool MappedBy(const GUNDAM::Match<Pattern, DataGraph> &) const = 0;
};

// x.A
template <class Pattern, class DataGraph>
class AttributeLiteral : public Literal<Pattern, DataGraph> {
 private:
  friend class LiteralInfo<Pattern, DataGraph>;

  using BaseLiteralType = Literal<Pattern, DataGraph>;

  using LiteralInfoType = typename BaseLiteralType::LiteralInfoType;

  using LiteralStandAloneInfoType =
      typename BaseLiteralType::LiteralStandAloneInfoType;

  using PatternVertexHandle = typename BaseLiteralType::PatternVertexHandle;
  using PatternVertexIDType = typename BaseLiteralType::PatternVertexIDType;

  using DataGraphVertexHandle = typename BaseLiteralType::DataGraphVertexHandle;
  using DataGraphVertexIDType = typename BaseLiteralType::DataGraphVertexIDType;

  using DataGraphEdgeHandle = typename BaseLiteralType::DataGraphEdgeHandle;
  using DataGraphEdgeIDType = typename BaseLiteralType::DataGraphEdgeIDType;

  using DataGraphVertexAttributeKeyType =
      typename BaseLiteralType::DataGraphVertexAttributeKeyType;

  PatternVertexHandle vertex_handle_;
  DataGraphVertexAttributeKeyType attr_key_;

  inline static bool Satisfy(const DataGraphVertexHandle& dst_handle,
                             const DataGraphVertexAttributeKeyType& attr_key) {
    auto vertex_attr_handle = dst_handle->FindAttribute(attr_key);
    if (!vertex_attr_handle) {
      // does not have attribute for this key, does not satisfy
      // this literal return false
      return false;
    }
    // contain this literal
    return true;
  }

  inline bool Satisfy(const DataGraphVertexHandle& dst_handle) const {
    return AttributeLiteral::Satisfy(dst_handle, this->attr_key_);
  }

  inline bool InnerUpdate(const DataGraphVertexHandle&  data_graph_vertex_handle,
                       std::set<DataGraphVertexIDType> *diff_vertex_attr_set) const {
    assert(data_graph_vertex_handle);
    auto [add_attr_handle, 
          add_attr_ret] = data_graph_vertex_handle->AddAttribute(this->attr_key_,
                                                                 std::string{"#"});
    if (!add_attr_ret) {
      // add fail, already have this attribute
      return false;
    }
    // added success, which means the does not have this attribute before
    if (diff_vertex_attr_set != nullptr) {
      diff_vertex_attr_set->insert(data_graph_vertex_handle->id());
    }
    return true;
  }

 public:
  AttributeLiteral(Pattern &pattern, PatternVertexIDType vertex_id,
                   DataGraphVertexAttributeKeyType attr_key)
      : vertex_handle_(pattern.FindVertex(vertex_id)), attr_key_(attr_key) {
    assert(this->vertex_handle_);
    return;
  }

  AttributeLiteral(Pattern &pattern, const LiteralInfoType &literal_info)
      : vertex_handle_(pattern.FindVertex(literal_info.x_id())),
        attr_key_(literal_info.x_attr_key()) {
    assert(literal_info.literal_type() == LiteralType::kAttrValueLiteral);
    assert(this->vertex_handle_);
    return;
  }

  ~AttributeLiteral() { return; }

  virtual bool Satisfy(
      const std::map<PatternVertexHandle, 
                   DataGraphVertexHandle> &match_result) const override {
    assert(this->MappedBy(match_result));
    DataGraphVertexHandle dst_handle =
        match_result.find(this->vertex_handle_)->second;
    assert(dst_handle);
    return this->Satisfy(dst_handle);
  }

  virtual bool Satisfy(const GUNDAM::Match<Pattern, DataGraph>& match) const override {
    assert(this->MappedBy(match));
    DataGraphVertexHandle dst_handle = match.MapTo(this->vertex_handle_);
    assert(dst_handle);
    return this->Satisfy(dst_handle);
  }

  virtual bool MappedBy(
      const std::map<PatternVertexHandle, DataGraphVertexHandle> &match_result)
      const override {
    assert(this->vertex_handle_);
    return match_result.find(this->vertex_handle_) != match_result.end();
  }

  virtual bool MappedBy(
      const GUNDAM::Match<Pattern, DataGraph> &match) const override {
    assert(this->vertex_handle_);
    return match.HasMap(this->vertex_handle_);
  }

  virtual bool Update(
      const std::map<PatternVertexHandle, DataGraphVertexHandle> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<DataGraphEdgeIDType> &edge_id_gen,
      std::set<DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<DataGraphEdgeIDType>   *diff_edge_set        = nullptr) const override {
    assert(this->MappedBy(match_result));
    assert(!this->Satisfy(match_result));

    assert(this->vertex_handle_);

    DataGraphVertexHandle data_graph_vertex_handle 
           = match_result.find(this->vertex_handle_)->second;
    assert(data_graph_vertex_handle);

    return this->InnerUpdate(data_graph_vertex_handle,
                             diff_vertex_attr_set);
  }

  virtual bool Update(
      const GUNDAM::Match<Pattern, DataGraph> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<DataGraphEdgeIDType> &edge_id_gen,
      std::set<DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<DataGraphEdgeIDType>   *diff_edge_set        = nullptr) const override {
    assert(this->MappedBy(match_result));
    assert(!this->Satisfy(match_result));

    assert(this->vertex_handle_);

    DataGraphVertexHandle data_graph_vertex_handle =
        match_result.MapTo(this->vertex_handle_);
    assert(data_graph_vertex_handle);

    return this->InnerUpdate(data_graph_vertex_handle,
                             diff_vertex_attr_set);
  }

  virtual Literal<Pattern, DataGraph> *Copy(Pattern &pattern) const override {
    assert(this->vertex_handle_);
    PatternVertexIDType vertex_id = this->vertex_handle_->id();
    return new AttributeLiteral<Pattern, DataGraph>(pattern, vertex_id,
                                                    this->attr_key_);
  }

  virtual void Print(std::ostream &out,
                     std::string gar_name = std::string()) const override {
    assert(this->vertex_handle_);
    out << "Attribute," << ProtectedSeparatorVal(this->vertex_handle_->id())
        << "," << ProtectedSeparatorVal(this->attr_key_) << ",,,,";
    if (gar_name != "") {
      out << "," << gar_name;
    }
    out << std::endl;
    return;
  }

  virtual void CalPivot(std::set<PatternVertexHandle> &pivot) const override {
    assert(this->vertex_handle_);
    pivot.insert(this->vertex_handle_);
    return;
  }

  virtual void CalPivot(
      std::vector<PatternVertexHandle> &pivot) const override {
    assert(this->vertex_handle_);
    pivot.emplace_back(this->vertex_handle_);
    return;
  }

  inline LiteralInfoType info() const override {
    assert(this->vertex_handle_);
    return LiteralInfoType(this->vertex_handle_->id(), this->attr_key_);
  }

  inline LiteralStandAloneInfoType stand_alone_info() const override {
    assert(this->vertex_handle_);
    return LiteralStandAloneInfoType(this->vertex_handle_->label(),
                                     this->attr_key_);
  }

  virtual enum LiteralType type() const override {
    return LiteralType::kAttrValueLiteral;
  }
};

// x -l-> y
template <class Pattern, class DataGraph>
class EdgeLiteral : public Literal<Pattern, DataGraph> {
 private:
  friend class LiteralInfo<Pattern, DataGraph>;

  using BaseLiteralType = Literal<Pattern, DataGraph>;

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

  inline static bool Satisfy(const DataGraphVertexHandle  &match_src_handle,
                             const DataGraphVertexHandle  &match_dst_handle,
                             const DataGraphEdgeLabelType &edge_label) {
    assert(match_src_handle);
    assert(match_dst_handle);
    for (auto it = match_src_handle->OutEdgeBegin(edge_label);
             !it.IsDone(); 
              it++) {
      if (it->dst_handle() == match_dst_handle) 
        return true;
    }
    return false;
  }

  inline bool Satisfy(const DataGraphVertexHandle &match_src_handle,
                      const DataGraphVertexHandle &match_dst_handle) const {

    return EdgeLiteral::Satisfy(match_src_handle,
                                match_dst_handle, this->edge_label_);
  }

  inline bool InnerUpdate(const DataGraphVertexHandle &match_src_handle,
                          const DataGraphVertexHandle &match_dst_handle,
                                            DataGraph &data_graph,
                    GUNDAM::SimpleArithmeticIDGenerator<
                                DataGraphEdgeIDType>  &edge_id_gen,
                       std::set<DataGraphEdgeIDType>  *diff_edge_set) const {
    assert(match_src_handle);
    assert(match_dst_handle);

    auto [edge_handle, 
          edge_ret] = data_graph.AddEdge(match_src_handle->id(), 
                                         match_dst_handle->id(),
                                         this->edge_label_, 
                                         edge_id_gen.GetID());
    assert(edge_handle);
    assert(edge_ret);
    if (!edge_ret) {
      return false;
    }
    if (diff_edge_set) {
      diff_edge_set->emplace(edge_handle->id());
    }
    return true;
  }

 public:
  EdgeLiteral(Pattern &pattern, PatternVertexIDType src_id,
              PatternVertexIDType dst_id, DataGraphEdgeLabelType edge_label)
      : src_handle_{pattern.FindVertex(src_id)},
        dst_handle_{pattern.FindVertex(dst_id)},
        edge_label_{edge_label} {
    assert(src_handle_);
    assert(dst_handle_);
    return;
  }

  EdgeLiteral(Pattern &pattern, const LiteralInfoType &literal_info)
      : src_handle_{pattern.FindVertex(literal_info.x_id())},
        dst_handle_{pattern.FindVertex(literal_info.y_id())},
        edge_label_{literal_info.edge_label()} {
    assert(literal_info.literal_type() == LiteralType::kEdgeLiteral);
    assert(src_handle_);
    assert(dst_handle_);
    return;
  }

  virtual bool Satisfy(
      const std::map<PatternVertexHandle, 
                   DataGraphVertexHandle> &match_result) const override {
    assert(this->MappedBy(match_result));
    auto match_src_handle = match_result.find(this->src_handle_)->second,
         match_dst_handle = match_result.find(this->dst_handle_)->second;
    assert(match_src_handle);
    assert(match_dst_handle);
    return this->Satisfy(match_src_handle,
                         match_dst_handle);
  }

  virtual bool Satisfy(
      const GUNDAM::Match<Pattern, DataGraph> &match) const override {
    assert(this->MappedBy(match));
    auto match_src_handle = match.MapTo(this->src_handle_),
         match_dst_handle = match.MapTo(this->dst_handle_);
    assert(match_src_handle);
    assert(match_dst_handle);
    return this->Satisfy(match_src_handle,
                         match_dst_handle);
  }

  virtual bool MappedBy(
      const std::map<PatternVertexHandle, 
                   DataGraphVertexHandle> &match_result)
      const override {
    return (match_result.find(this->src_handle_) != match_result.end()) &&
           (match_result.find(this->dst_handle_) != match_result.end());
  }

  virtual bool MappedBy(
      const GUNDAM::Match<Pattern, DataGraph> &match) const override {
    return match.HasMap(this->src_handle_)
        && match.HasMap(this->dst_handle_);
  }

  virtual bool Update(
      const std::map<PatternVertexHandle, DataGraphVertexHandle> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<DataGraphEdgeIDType> &edge_id_gen,
      std::set<DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<DataGraphEdgeIDType>   *diff_edge_set        = nullptr) const override {
    assert(this->MappedBy(match_result));
    assert(!this->Satisfy(match_result));

    auto match_src_handle = match_result.find(this->src_handle_)->second,
         match_dst_handle = match_result.find(this->dst_handle_)->second;
    assert(match_src_handle);
    assert(match_dst_handle);

    return this->InnerUpdate(match_src_handle,
                             match_dst_handle,
                             data_graph,
                             edge_id_gen,
                             diff_edge_set);
  }

  virtual bool Update(
      const GUNDAM::Match<Pattern, DataGraph> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<DataGraphEdgeIDType> &edge_id_gen,
      std::set<DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<DataGraphEdgeIDType>   *diff_edge_set        = nullptr) const override {
    assert(this->MappedBy(match_result));
    assert(!this->Satisfy(match_result));

    auto match_src_handle = match_result.MapTo(this->src_handle_),
         match_dst_handle = match_result.MapTo(this->dst_handle_);
    assert(match_src_handle);
    assert(match_dst_handle);

    return this->InnerUpdate(match_src_handle,
                             match_dst_handle,
                             data_graph,
                             edge_id_gen,
                             diff_edge_set);
  }

  virtual Literal<Pattern, DataGraph> *Copy(Pattern &pattern) const override {
    assert(this->src_handle_);
    assert(this->dst_handle_);
    PatternVertexIDType src_id = this->src_handle_->id();
    PatternVertexIDType dst_id = this->dst_handle_->id();
    return new EdgeLiteral<Pattern, DataGraph>(pattern, src_id, dst_id,
                                               this->edge_label_);
  }

  virtual void Print(std::ostream &out,
                     std::string gar_name = std::string()) const override {
    out << "Edge," << ProtectedSeparatorVal(this->src_handle_->id()) << ",,"
        << ProtectedSeparatorVal(this->dst_handle_->id()) << ",,"
        << ProtectedSeparatorVal(this->edge_label_) << ",";
    if (gar_name != "") {
      out << "," << gar_name;
    }
    out << std::endl;
  }

  virtual void CalPivot(std::set<PatternVertexHandle> &pivot) const override {
    pivot.insert(this->src_handle_);
    pivot.insert(this->dst_handle_);
    return;
  }

  virtual void CalPivot(
      std::vector<PatternVertexHandle> &pivot) const override {
    pivot.emplace_back(this->src_handle_);
    pivot.emplace_back(this->dst_handle_);
    return;
  }

  inline LiteralInfoType info() const override {
    return LiteralInfoType(this->src_handle_->id(), this->dst_handle_->id(),
                           this->edge_label_);
  }

  inline LiteralStandAloneInfoType stand_alone_info() const override {
    return LiteralStandAloneInfoType(this->src_handle_->label(),
                                     this->dst_handle_->label(),
                                     this->edge_label_);
  }

  virtual enum LiteralType type() const override {
    return LiteralType::kEdgeLiteral;
  }

 private:
  PatternVertexHandle src_handle_, dst_handle_;

  DataGraphEdgeLabelType edge_label_;
};

// x.A = y.B
template <class Pattern, class DataGraph>
class VariableLiteral : public Literal<Pattern, DataGraph> {
 private:
  friend class LiteralInfo<Pattern, DataGraph>;

  using LiteralInfoType = typename Literal<Pattern, DataGraph>::LiteralInfoType;

  using LiteralStandAloneInfoType = LiteralStandAloneInfo<Pattern, DataGraph>;

  using BaseLiteralType = Literal<Pattern, DataGraph>;
  using PatternVertexHandle = typename BaseLiteralType::PatternVertexHandle;
  using PatternVertexIDType = typename BaseLiteralType::PatternVertexIDType;

  using DataGraphVertexHandle = typename BaseLiteralType::DataGraphVertexHandle;
  using DataGraphVertexIDType = typename BaseLiteralType::DataGraphVertexIDType;
  using DataGraphVertexAttributeKeyType =
      typename BaseLiteralType::DataGraphVertexAttributeKeyType;

  using DataGraphEdgeIDType = typename BaseLiteralType::DataGraphEdgeIDType;
  using DataGraphEdgeHandle = typename BaseLiteralType::DataGraphEdgeHandle;

  std::pair<PatternVertexHandle, DataGraphVertexAttributeKeyType> x_, y_;

  PatternVertexHandle x_handle_, y_handle_;

  DataGraphVertexAttributeKeyType x_attr_key_, y_attr_key_;

  inline static bool Satisfy(const DataGraphVertexHandle  &match_x_handle,
                             const DataGraphVertexHandle  &match_y_handle,
                             const DataGraphVertexAttributeKeyType &x_attr_key,
                             const DataGraphVertexAttributeKeyType &y_attr_key) {
    assert(match_x_handle);
    assert(match_y_handle);
    auto x_attr_handle = match_x_handle->FindAttribute(x_attr_key);
    auto y_attr_handle = match_y_handle->FindAttribute(y_attr_key);
    if (x_attr_handle.IsNull() || y_attr_handle.IsNull()) {
      // one of the vertex does not contain the specified attribute
      return false;
    }
    assert(x_attr_handle && y_attr_handle);
    GUNDAM::BasicDataType x_value_type = x_attr_handle->value_type();
    GUNDAM::BasicDataType y_value_type = y_attr_handle->value_type();
    if (x_value_type != y_value_type) {
      // these two attributes have different types
      // therefore, these two attribute are not the same
      return false;
    }
    switch (x_value_type) {
      case GUNDAM::BasicDataType::kTypeInt:
        if (x_attr_handle->template const_value<int>() !=
            y_attr_handle->template const_value<int>()) {
          return false;
        }
        return true;
      case GUNDAM::BasicDataType::kTypeDouble:
        if (x_attr_handle->template const_value<double>() !=
            y_attr_handle->template const_value<double>()) {
          return false;
        }
        return true;
      case GUNDAM::BasicDataType::kTypeString:
        if (x_attr_handle->template const_value<std::string>() !=
            y_attr_handle->template const_value<std::string>()) {
          return false;
        }
        return true;
      default:
        // unknown data type
        assert(false);
        break;
    }
    return false;
  }

  inline bool Satisfy(
      const DataGraphVertexHandle &match_x_handle,
      const DataGraphVertexHandle &match_y_handle) const {
      return VariableLiteral::Satisfy(match_x_handle, 
                                      match_y_handle,
                                      this->x_attr_key_,
                                      this->y_attr_key_);
  }

  inline bool InnerUpdate(const DataGraphVertexHandle &x_handle,
                          const DataGraphVertexHandle &y_handle,
                       std::set<DataGraphVertexIDType> *diff_vertex_attr_set) const {
    assert(x_handle);
    assert(y_handle);

    auto x_attr_handle = x_handle->FindAttribute(this->x_attr_key_);
    auto y_attr_handle = y_handle->FindAttribute(this->y_attr_key_);
    if (x_attr_handle.IsNull()) {
      bool x_attr_ret = false;
      std::tie(x_attr_handle, x_attr_ret) =
          x_handle->AddAttribute(this->x_attr_key_, std::string{"#"});
      assert(x_attr_handle);
      assert(x_attr_ret);
    }
    if (y_attr_handle.IsNull()) {
      bool y_attr_ret = false;
      std::tie(y_attr_handle, y_attr_ret) =
          y_handle->AddAttribute(this->y_attr_key_, std::string{"#"});
      assert(y_attr_handle);
      assert(y_attr_ret);
    }
    assert(x_attr_handle);
    assert(y_attr_handle);

    bool x_empty_flag =
        (x_attr_handle->value_type() == GUNDAM::BasicDataType::kTypeString) &&
        (x_attr_handle->template value<std::string>() == "#");
    bool y_empty_flag =
        (y_attr_handle->value_type() == GUNDAM::BasicDataType::kTypeString) &&
        (y_attr_handle->template value<std::string>() == "#");

    if (!x_empty_flag && !y_empty_flag) {
      // can do nothing
      // whethere there are x.A == y.B or x.A != y.B
      return false;
    }
    if (x_empty_flag && y_empty_flag) {
      // both evaluation method are empty
      // assert(false);
      return true;
    }
    auto src_attr_handle = x_empty_flag ? y_attr_handle : x_attr_handle;
    auto dst_handle = x_empty_flag ? x_handle : y_handle;
    auto &dst_attr_key =
        x_empty_flag ? x_attr_handle->key() : y_attr_handle->key();

    auto value_type = src_attr_handle->value_type();

    switch (value_type) {
      case GUNDAM::BasicDataType::kTypeInt:
        dst_handle->SetAttribute(dst_attr_key,
                                 src_attr_handle->template value<int>());
        break;
      case GUNDAM::BasicDataType::kTypeDouble:
        dst_handle->SetAttribute(dst_attr_key,
                                 src_attr_handle->template value<double>());
        break;
      case GUNDAM::BasicDataType::kTypeString:
        dst_handle->SetAttribute(
            dst_attr_key, src_attr_handle->template value<std::string>());
        break;
      default:
        // unknown data type
        assert(false);
        break;
    }
    if (diff_vertex_attr_set != nullptr) {
      diff_vertex_attr_set->insert(dst_handle->id());
    }
    return true;
  }

 public:
  VariableLiteral(Pattern &pattern, PatternVertexIDType x_id,
                  DataGraphVertexAttributeKeyType x_attr_key,
                  PatternVertexIDType y_id,
                  DataGraphVertexAttributeKeyType y_attr_key)
      : x_handle_(pattern.FindVertex(x_id)),
        y_handle_(pattern.FindVertex(y_id)),
        x_attr_key_(x_attr_key),
        y_attr_key_(y_attr_key) {
    assert(this->x_handle_);
    assert(this->y_handle_);
    return;
  }

  VariableLiteral(Pattern &pattern, const LiteralInfoType &literal_info)
      : x_handle_(pattern.FindVertex(literal_info.x_id())),
        y_handle_(pattern.FindVertex(literal_info.y_id())),
        x_attr_key_(literal_info.x_attr_key()),
        y_attr_key_(literal_info.y_attr_key()) {
    assert(literal_info.literal_type() == LiteralType::kVariableLiteral);
    assert(this->x_handle_);
    assert(this->y_handle_);
    return;
  }

  ~VariableLiteral() { return; }

  virtual bool Satisfy(
      const std::map<PatternVertexHandle, 
                   DataGraphVertexHandle> &match_result) const override {
    assert(this->MappedBy(match_result));
    DataGraphVertexHandle match_x_handle =
        match_result.find(this->x_handle_)->second;
    DataGraphVertexHandle match_y_handle =
        match_result.find(this->y_handle_)->second;
    assert(match_x_handle);
    assert(match_y_handle);
    return this->Satisfy(match_x_handle,
                         match_y_handle);
  }

  virtual bool Satisfy(
      const GUNDAM::Match<Pattern, DataGraph> &match) const override {
    assert(this->MappedBy(match));
    DataGraphVertexHandle match_x_handle
            = match.MapTo(this->x_handle_);
    DataGraphVertexHandle match_y_handle
            = match.MapTo(this->y_handle_);
    assert(match_x_handle);
    assert(match_y_handle);
    return this->Satisfy(match_x_handle,
                         match_y_handle);
  }

  virtual bool MappedBy(
      const std::map<PatternVertexHandle, DataGraphVertexHandle> &match_result)
      const override {
    return (match_result.find(this->x_handle_) != match_result.end()) &&
           (match_result.find(this->y_handle_) != match_result.end());
  }

  virtual bool MappedBy(
      const GUNDAM::Match<Pattern, DataGraph> &match) const override {
    return match.HasMap(this->x_handle_) && match.HasMap(this->y_handle_);
  }

  virtual bool Update(
      const std::map<PatternVertexHandle, DataGraphVertexHandle> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<DataGraphEdgeIDType> &edge_id_gen,
      std::set<DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<DataGraphEdgeIDType>   *diff_edge_set        = nullptr) const override {
    assert(this->MappedBy(match_result));
    assert(!this->Satisfy(match_result));
    DataGraphVertexHandle x_handle = match_result.find(this->x_handle_)->second;
    DataGraphVertexHandle y_handle = match_result.find(this->y_handle_)->second;
    assert(x_handle);
    assert(y_handle);

    return this->InnerUpdate(x_handle,
                             y_handle,
                             diff_vertex_attr_set);
  }

  virtual bool Update(
      const GUNDAM::Match<Pattern, DataGraph> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<DataGraphEdgeIDType> &edge_id_gen,
      std::set<DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<DataGraphEdgeIDType>   *diff_edge_set        = nullptr) const override {
    assert(this->MappedBy(match_result));
    assert(!this->Satisfy(match_result));
    DataGraphVertexHandle x_handle = match_result.MapTo(this->x_handle_);
    DataGraphVertexHandle y_handle = match_result.MapTo(this->y_handle_);
    assert(x_handle);
    assert(y_handle);

    return this->InnerUpdate(x_handle,
                             y_handle,
                             diff_vertex_attr_set);
  }

  virtual Literal<Pattern, DataGraph> *Copy(Pattern &pattern) const override {
    assert(this->x_handle_);
    assert(this->y_handle_);
    PatternVertexIDType x_id = this->x_handle_->id();
    PatternVertexIDType y_id = this->y_handle_->id();
    return new VariableLiteral<Pattern, DataGraph>(
        pattern, x_id, this->x_attr_key_, y_id, this->y_attr_key_);
  }

  virtual void Print(std::ostream &out,
                     std::string gar_name = std::string()) const override {
    out << "Variable," << ProtectedSeparatorVal(this->x_handle_->id()) << ","
        << ProtectedSeparatorVal(this->x_attr_key_) << ","
        << ProtectedSeparatorVal(this->y_handle_->id()) << ","
        << ProtectedSeparatorVal(this->y_attr_key_) << ","
        << ",";
    if (gar_name != "") {
      out << "," << gar_name;
    }
    out << std::endl;
    return;
  }

  virtual void CalPivot(std::set<PatternVertexHandle> &pivot) const override {
    pivot.insert(this->x_handle_);
    pivot.insert(this->y_handle_);
    return;
  }

  virtual void CalPivot(
      std::vector<PatternVertexHandle> &pivot) const override {
    pivot.emplace_back(this->x_handle_);
    pivot.emplace_back(this->y_handle_);
    return;
  }

  inline LiteralInfoType info() const override {
    return LiteralInfoType(this->x_handle_->id(), this->x_attr_key_,
                           this->y_handle_->id(), this->y_attr_key_);
  }

  inline LiteralStandAloneInfoType stand_alone_info() const override {
    return LiteralStandAloneInfoType(
        this->x_handle_->label(), this->x_attr_key_, this->y_handle_->label(),
        this->y_attr_key_);
  }

  virtual enum LiteralType type() const override {
    return LiteralType::kVariableLiteral;
  }
};

// x.A = c
template <class Pattern, class DataGraph>
class ConstantLiteral : public Literal<Pattern, DataGraph> {
 private:
  friend class LiteralInfo<Pattern, DataGraph>;

  using LiteralInfoType = typename Literal<Pattern, DataGraph>::LiteralInfoType;

  using LiteralStandAloneInfoType = LiteralStandAloneInfo<Pattern, DataGraph>;

  using BaseLiteralType = Literal<Pattern, DataGraph>;
  using PatternVertexHandle = typename BaseLiteralType::PatternVertexHandle;
  using PatternVertexIDType = typename BaseLiteralType::PatternVertexIDType;

  using DataGraphVertexHandle = typename BaseLiteralType::DataGraphVertexHandle;
  using DataGraphVertexIDType = typename BaseLiteralType::DataGraphVertexIDType;
  using DataGraphVertexAttributeKeyType =
      typename BaseLiteralType::DataGraphVertexAttributeKeyType;

  using DataGraphEdgeHandle = typename BaseLiteralType::DataGraphEdgeHandle;
  using DataGraphEdgeIDType = typename BaseLiteralType::DataGraphEdgeIDType;

  inline static bool Satisfy(const DataGraphVertexHandle& match_x_handle,
                             const DataGraphVertexAttributeKeyType& x_attr_key,
                             const std::string& c_str,
                             const GUNDAM::BasicDataType &data_type) {
    assert(match_x_handle);
    auto x_attr_handle = match_x_handle->FindAttribute(x_attr_key);
    if (x_attr_handle.IsNull()) {
      // does not have this literal, does not satisfy this literal
      return false;
    }
    if (x_attr_handle->value_type() != data_type) {
      // type are not the same
      return false;
    }
    if (x_attr_handle->value_str() != c_str) {
      // type is the same, value is not the same
      return false;
    }
    return true;
  }

  inline bool Satisfy(const DataGraphVertexHandle match_x_handle) const {
    return ConstantLiteral::Satisfy(match_x_handle,
                                    this->x_attr_key_,
                                    this->c_str_,
                                    this->data_type_);
  }

  inline bool InnerUpdate(const DataGraphVertexHandle &x_handle,
                       std::set<DataGraphVertexIDType> *diff_vertex_attr_set) const {
    assert(x_handle);

    auto x_attr_handle = x_handle->FindAttribute(this->x_attr_key_);
    if (!x_attr_handle) {
      // does not have this attribute
      bool x_attr_ret = false;
      std::tie(x_attr_handle, 
               x_attr_ret) = x_handle->AddAttribute(this->x_attr_key_, 
                                                    this->data_type_, 
                                                    this->c_str_);
      // should added successfully
      assert(x_attr_ret);
      return true;
    }
    assert(x_attr_handle);
    // have this value
    GUNDAM::BasicDataType x_value_type = x_attr_handle->value_type();
    bool x_empty_flag = (x_value_type == GUNDAM::BasicDataType::kTypeString &&
                         x_attr_handle->template value<std::string>() == "#");
    if (!x_empty_flag) {
      // whether there are x_attr_handle->value_str() != this->c_str_
      // or not
      return false;
    }
    auto [x_add_attr_handle, 
          x_add_attr_ret]
        = x_handle->SetAttribute(this->x_attr_key_, this->c_str_);
    assert(x_add_attr_ret);
    if (diff_vertex_attr_set != nullptr) {
      diff_vertex_attr_set->insert(x_handle->id());
    }
    return true;
  }

 private:
  PatternVertexHandle x_handle_;

  DataGraphVertexAttributeKeyType x_attr_key_;

  std::string c_str_;

  GUNDAM::BasicDataType data_type_;

 public:
  template <typename ConstantType>
  ConstantLiteral(Pattern &pattern, PatternVertexIDType x_id,
                  DataGraphVertexAttributeKeyType attr_key, ConstantType c)
      : x_handle_(pattern.FindVertex(x_id)),
        x_attr_key_(attr_key),
        c_str_(GUNDAM::ToString(c)),
        data_type_(GUNDAM::TypeToEnum(c)) {
    assert(this->x_handle_);
    return;
  }

  ConstantLiteral(Pattern &pattern, PatternVertexIDType x_id,
                  DataGraphVertexAttributeKeyType attr_key, std::string c_str,
                  GUNDAM::BasicDataType data_type)
      : x_handle_(pattern.FindVertex(x_id)),
        x_attr_key_(attr_key),
        c_str_(c_str),
        data_type_(data_type) {
    assert(this->x_handle_);
    return;
  }

  ConstantLiteral(Pattern &pattern, const LiteralInfoType &literal_info)
      : x_handle_(pattern.FindVertex(literal_info.x_id())),
        x_attr_key_(literal_info.x_attr_key()),
        c_str_(literal_info.c_str()),
        data_type_(literal_info.data_type()) {
    assert(literal_info.literal_type() == LiteralType::kConstantLiteral);
    assert(this->x_handle_);
    return;
  }

  ~ConstantLiteral() { return; }

  virtual bool Satisfy(
      const std::map<PatternVertexHandle, 
                   DataGraphVertexHandle> &match_result) const override {
    assert(this->MappedBy(match_result));
    DataGraphVertexHandle match_x_handle =
        match_result.find(this->x_handle_)->second;
    assert(match_x_handle);
    return this->Satisfy(match_x_handle);
  }

  virtual bool Satisfy(
      const GUNDAM::Match<Pattern, DataGraph> &match) const override {
    assert(this->MappedBy(match));
    DataGraphVertexHandle match_x_handle = match.MapTo(this->x_handle_);
    assert(match_x_handle);
    return this->Satisfy(match_x_handle);
  }

  virtual bool MappedBy(
      const std::map<PatternVertexHandle, DataGraphVertexHandle> &match_result)
      const override {
    return match_result.find(this->x_handle_) != match_result.end();
  }

  virtual bool MappedBy(
      const GUNDAM::Match<Pattern, DataGraph> &match) const override {
    return match.HasMap(this->x_handle_->id());
  }

  virtual bool Update(
      const std::map<PatternVertexHandle, DataGraphVertexHandle> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<DataGraphEdgeIDType> &edge_id_gen,
      std::set<DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<DataGraphEdgeIDType>   *diff_edge_set        = nullptr) const override {
    assert(this->MappedBy(match_result));
    assert(!this->Satisfy(match_result));
    DataGraphVertexHandle x_handle = match_result.find(this->x_handle_)->second;
    assert(x_handle);

    return this->InnerUpdate(x_handle, diff_vertex_attr_set);
  }

  virtual bool Update(
      const GUNDAM::Match<Pattern, DataGraph> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<DataGraphEdgeIDType> &edge_id_gen,
      std::set<DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<DataGraphEdgeIDType>   *diff_edge_set        = nullptr) const override {
    assert(this->MappedBy(match_result));
    assert(!this->Satisfy(match_result));
    DataGraphVertexHandle x_handle = match_result.MapTo(this->x_handle_);
    assert(x_handle);

    return this->InnerUpdate(x_handle, diff_vertex_attr_set);
  }

  virtual Literal<Pattern, DataGraph> *Copy(Pattern &pattern) const override {
    assert(this->x_handle_);
    PatternVertexIDType x_id = this->x_handle_->id();
    return new ConstantLiteral<Pattern, DataGraph>(
        pattern, x_id, this->x_attr_key_, this->c_str_, this->data_type_);
  }

  virtual void Print(std::ostream &out,
                     std::string gar_name = std::string()) const override {
    std::string value_type = GUNDAM::EnumToString(this->data_type_);
    out << "Constant," << ProtectedSeparatorVal(this->x_handle_->id()) << ","
        << ProtectedSeparatorVal(this->x_attr_key_) << ","
        << ",,," << ProtectedSeparatorVal(this->c_str_) << ";" << value_type;
    if (gar_name != "") {
      out << "," << gar_name;
    }
    out << std::endl;
  }

  virtual void CalPivot(std::set<PatternVertexHandle> &pivot) const override {
    pivot.insert(this->x_handle_);
    return;
  }

  virtual void CalPivot(
      std::vector<PatternVertexHandle> &pivot) const override {
    pivot.emplace_back(this->x_handle_);
    return;
  }

  inline LiteralInfoType info() const override {
    return LiteralInfoType(this->x_handle_->id(), this->x_attr_key_,
                           GUNDAM::ToString(this->c_str_), this->data_type_);
  }

  inline LiteralStandAloneInfoType stand_alone_info() const override {
    return LiteralStandAloneInfoType(
        this->x_handle_->label(), this->x_attr_key_,
        GUNDAM::ToString(this->c_str_), this->data_type_);
  }

  virtual enum LiteralType type() const override {
    return LiteralType::kConstantLiteral;
  }
};

// x.id = y.id
template <class Pattern, class DataGraph>
class KeyLiteral : public Literal<Pattern, DataGraph> {
 private:
  friend class LiteralInfo<Pattern, DataGraph>;
  
  using LiteralInfoType = typename Literal<Pattern, DataGraph>::LiteralInfoType;

  using LiteralStandAloneInfoType 
      = LiteralStandAloneInfo<Pattern, DataGraph>;

  using BaseLiteralType = Literal<Pattern, DataGraph>;
  using PatternVertexHandle = typename BaseLiteralType::PatternVertexHandle;
  using PatternVertexIDType = typename BaseLiteralType::PatternVertexIDType;

  using DataGraphVertexHandle = typename BaseLiteralType::DataGraphVertexHandle;
  using DataGraphVertexIDType = typename BaseLiteralType::DataGraphVertexIDType;
  using DataGraphVertexAttributeKeyType =
      typename BaseLiteralType::DataGraphVertexAttributeKeyType;

  using DataGraphEdgeHandle = typename BaseLiteralType::DataGraphEdgeHandle;
  using DataGraphEdgeIDType = typename BaseLiteralType::DataGraphEdgeIDType;

  inline static bool Satisfy(const DataGraphVertexHandle &match_x_handle,
                             const DataGraphVertexHandle &match_y_handle) {
    return match_x_handle == match_y_handle;
  }

 private:
  PatternVertexHandle x_handle_,
                      y_handle_;

 public:
  KeyLiteral(Pattern &pattern, 
             PatternVertexIDType x_id,
             PatternVertexIDType y_id)
      : x_handle_(pattern.FindVertex(x_id)),
        y_handle_(pattern.FindVertex(y_id)) {
    assert(this->x_handle_);
    assert(this->y_handle_);
    return;
  }

  KeyLiteral(Pattern &pattern, const LiteralInfoType &literal_info)
          : x_handle_(pattern.FindVertex(literal_info.x_id())),
            y_handle_(pattern.FindVertex(literal_info.x_id())) {
    assert(literal_info.literal_type() == LiteralType::kKeyLiteral);
    assert(this->x_handle_);
    assert(this->y_handle_);
    return;
  }

  ~KeyLiteral() { 
    return;
  }

  virtual bool Satisfy(
      const std::map<PatternVertexHandle, 
                   DataGraphVertexHandle> &match_result) const override {
    assert(this->MappedBy(match_result));
    return match_result.find(this->x_handle_)->second->id()
        == match_result.find(this->y_handle_)->second->id();
  }

  virtual bool Satisfy(
      const GUNDAM::Match<Pattern, DataGraph> &match) const override {
    assert(this->MappedBy(match));
    return match.MapTo(this->x_handle_)->id()
        == match.MapTo(this->y_handle_)->id();
  }

  virtual bool MappedBy(
      const std::map<PatternVertexHandle, 
                   DataGraphVertexHandle> &match_result) const override {
    return (match_result.find(this->x_handle_) != match_result.end())
        && (match_result.find(this->y_handle_) != match_result.end());
  }

  virtual bool MappedBy(
      const GUNDAM::Match<Pattern, DataGraph> &match) const override {
    return match.HasMap(this->x_handle_)
        && match.HasMap(this->y_handle_);
  }

  virtual bool Update(
      const std::map<PatternVertexHandle, 
                   DataGraphVertexHandle> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<DataGraphEdgeIDType> &edge_id_gen,
      std::set<DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<DataGraphEdgeIDType>   *diff_edge_set        = nullptr) const override {
    // merge two vertexes, does not supported now
    assert(false);
    return false;
  }

  virtual bool Update(
      const GUNDAM::Match<Pattern, DataGraph> &match_result,
      DataGraph &data_graph,
      GUNDAM::SimpleArithmeticIDGenerator<DataGraphEdgeIDType> &edge_id_gen,
      std::set<DataGraphVertexIDType> *diff_vertex_attr_set = nullptr,
      std::set<DataGraphEdgeIDType>   *diff_edge_set        = nullptr) const override {
      // merge two vertexes, does not supported now
    assert(false);
    return false;
  }

  virtual Literal<Pattern, DataGraph> *Copy(Pattern &pattern) const override {
    assert(this->x_handle_);
    assert(this->y_handle_);
    PatternVertexIDType x_id = this->x_handle_->id();
    PatternVertexIDType y_id = this->y_handle_->id();
    return new KeyLiteral<Pattern, DataGraph>(pattern, x_id, y_id);
  }

  virtual void Print(std::ostream &out,
                     std::string gar_name = std::string()) const override {
  }

  virtual void CalPivot(std::set<PatternVertexHandle> &pivot) const override {
    pivot.emplace(this->x_handle_);
    pivot.emplace(this->y_handle_);
    return;
  }

  virtual void CalPivot(
      std::vector<PatternVertexHandle> &pivot) const override {
    pivot.emplace_back(this->x_handle_);
    pivot.emplace_back(this->y_handle_);
    return;
  }

  inline LiteralInfoType info() const override {
    return LiteralInfoType(this->x_handle_->id(), 
                           this->y_handle_->id());
  }

  inline LiteralStandAloneInfoType stand_alone_info() const override {
    return LiteralStandAloneInfoType(
        this->x_handle_->label(),
        this->y_handle_->label());
  }

  virtual enum LiteralType type() const override {
    return LiteralType::kKeyLiteral;
  }
};

}  // namespace gar

#endif