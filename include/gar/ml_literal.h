#ifndef EXAMPLES_ANALYTICAL_ML_LITERAL_H_
#define EXAMPLES_ANALYTICAL_ML_LITERAL_H_

#include "literal.h"
#include "literal_info.h"

//#include <xmlrpc-c/base.hpp>
//#include <xmlrpc-c/client_simple.hpp>
//#include <xmlrpc-c/girerr.hpp>

namespace gar {

template <class Pattern, class DataGraph>
class MlLiteral : public Literal<Pattern, DataGraph> {
 private:
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
  MlLiteral(Pattern &pattern, PatternVertexIDType src_id,
            PatternVertexIDType dst_id, DataGraphEdgeLabelType edge_label,
            std::string module_url) {
    this->src_handle_ = pattern.FindVertex(src_id);
    this->dst_handle_ = pattern.FindVertex(dst_id);
    this->edge_label_ = edge_label;
    this->module_url_ = module_url;
  }

  MlLiteral(Pattern &pattern, const LiteralInfoType &literal_info)
      : src_handle_{pattern.FindVertex(literal_info.x_id())},
        dst_handle_{pattern.FindVertex(literal_info.y_id())},
        edge_label_{literal_info.edge_label()},
        module_url_{literal_info.module_url()} {
    assert(literal_info.literal_type() == LiteralType::kMlLiteral);
    assert(src_handle_);
    assert(dst_handle_);
    return;
  }

  ~MlLiteral() {
    return;
  }

  virtual bool Satisfy(const std::map<PatternVertexHandle, 
                                    DataGraphVertexHandle> &match_result) const override {
    //xmlrpc_c::clientSimple rpc_client;
    //xmlrpc_c::value result;
    //xmlrpc_c::paramList para;

    //rpc_client.call(module_url_, method_, para, &result);

    //string const str = xmlrpc_c::value_string(result);
    return false;
  }

  virtual bool Satisfy(const GUNDAM::Match<Pattern, DataGraph> &match) const override {
    //xmlrpc_c::clientSimple rpc_client;
    //xmlrpc_c::value result;
    //xmlrpc_c::paramList para;

    //rpc_client.call(module_url_, method_, para, &result);

    //string const str = xmlrpc_c::value_string(result);
    return false;
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
    PatternVertexIDType src_id = this->src_handle_->id();
    PatternVertexIDType dst_id = this->dst_handle_->id();
    Literal<Pattern, DataGraph> *literal_ptr =
        new MlLiteral<Pattern, DataGraph>(pattern, src_id, dst_id,
                                          this->edge_label_, 
                                          this->module_url_);
    return literal_ptr;
  }

  virtual void Print(std::ostream &out,
                     std::string gar_name = std::string()) const override {
    out << "Ml," << ProtectedSeparatorVal(this->src_handle_->id()) << ",,"
                 << ProtectedSeparatorVal(this->dst_handle_->id()) << ",,"
                 << ProtectedSeparatorVal(this->edge_label_) << ",";
    if (gar_name != "") {
      out << "," << gar_name;
    }
    out << std::endl;
    return;
  }

  virtual void CalPivot(std::set<PatternVertexHandle> &pivot) const override {
    pivot.emplace(this->src_handle_);
    pivot.emplace(this->dst_handle_);
    return;
  }

  virtual void CalPivot(std::vector<PatternVertexHandle> &pivot) const override{
    pivot.emplace_back(this->src_handle_);
    pivot.emplace_back(this->dst_handle_);
    return;
  }

  virtual LiteralInfoType info() const override {
    return LiteralInfoType(this->src_handle_->id(),
                           this->dst_handle_->id(),
                           this->edge_label_,
                           this->module_url_);
  }

  virtual LiteralStandAloneInfoType stand_alone_info() const override {
    return LiteralStandAloneInfoType(this->src_handle_->label(),
                                     this->dst_handle_->label(),
                                     this->edge_label_,
                                     this->module_url_);
  }

  virtual enum LiteralType type() const override{
    return LiteralType::kMlLiteral;
  }
        
  virtual bool MappedBy(
      const std::map<PatternVertexHandle, 
                   DataGraphVertexHandle> &match_result) const override {
    assert(this->src_handle_);
    assert(this->dst_handle_);
    return (match_result.find(this->src_handle_) != match_result.end())
        && (match_result.find(this->dst_handle_) != match_result.end());
  }

  virtual bool MappedBy(const GUNDAM::Match<Pattern, DataGraph> &match) const override{
    return match.HasMap(this->src_handle_) 
        && match.HasMap(this->dst_handle_);
  }

 private:
  std::string module_url_;
  std::string method_;
  PatternVertexHandle src_handle_, 
                      dst_handle_;
  DataGraphEdgeLabelType edge_label_;
};

}  // namespace gar

#endif