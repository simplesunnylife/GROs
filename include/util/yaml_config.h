#ifndef UTIL_YAML_CONFIG_H_
#define UTIL_YAML_CONFIG_H_

#include "yaml-cpp/yaml.h"
#include "util/log.h"

#include "gar/literal_stand_alone_info.h"

namespace util{

inline std::pair<std::tuple<std::string,
                            std::string,
                            std::string>,
                  bool> GetPatternPathInfo(const YAML::Node& pattern_path_config) {
  if (!pattern_path_config["VFile"]) {
    util::Error( "cannot get gar v file" );
    return std::pair(std::tuple("","",""), false);
  }
  const std::string kPatternPathVFile 
         = pattern_path_config["VFile"].as<std::string>();

  if (!pattern_path_config["EFile"]) {
    util::Error( "cannot get gar e file" );
    return std::pair(std::tuple("","",""), false);
  }
  const std::string kPatternPathEFile  
         = pattern_path_config["EFile"].as<std::string>();

  std::string gar_pivot_id_str = "";
  if (pattern_path_config["PivotId"]) {
    gar_pivot_id_str = pattern_path_config["PivotId"].as<std::string>();
  }

  return std::pair(std::tuple(kPatternPathVFile,
                              kPatternPathEFile,
                              gar_pivot_id_str), true);
}
  
inline std::pair<std::tuple<std::string,
                            std::string,
                            std::string>, bool> 
  GetDataGraphInfoFromYaml(const YAML::Node& data_graph_info_node) {
  if (!data_graph_info_node["VFile"]) {
    util::Error("cannot get data graph v file");
    return std::pair(std::tuple("","",""), false);
  }
  const std::string kDataGraphPathVFile 
          = data_graph_info_node["VFile"].as<std::string>();

  if (!data_graph_info_node["EFile"]) {
    util::Error("cannot get data graph e file");
    return std::pair(std::tuple("","",""), false);
  }
  const std::string kDataGraphPathEFile 
          = data_graph_info_node["EFile"].as<std::string>();

  if (!data_graph_info_node["MlLiteralEdgesFile"]){
    return std::pair(std::tuple(kDataGraphPathVFile,
                                kDataGraphPathEFile, ""),
                      true);
  }
  const std::string kDataGraphPathMlLiteralEdgesFile
          = data_graph_info_node["MlLiteralEdgesFile"].as<std::string>();
  return std::pair(std::tuple(kDataGraphPathVFile,
                              kDataGraphPathEFile, 
                              kDataGraphPathMlLiteralEdgesFile),
                    true);
}

template <typename DataGraphType>
inline std::pair<std::tuple<typename DataGraphType::VertexType::LabelType,
                            typename DataGraphType::  EdgeType::LabelType,
                            typename DataGraphType::VertexType::LabelType>, bool>
  EdgeTypeFromYaml(const YAML::Node& edge_type_config) {

  using VertexLabelType = typename DataGraphType::VertexType::LabelType;
  using   EdgeLabelType = typename DataGraphType::  EdgeType::LabelType;

  using EdgeType = std::tuple<VertexLabelType,
                                EdgeLabelType,
                              VertexLabelType>;

  EdgeType edge_type;

  if (!edge_type_config["SrcLabel"]){
    std::cout << "cannot get the src label of specified edge type" 
              << std::endl;
    return std::pair(EdgeType(), false);
  }
  std::get<0>(edge_type) = edge_type_config["SrcLabel"].as<VertexLabelType>();

  if (!edge_type_config["EdgeLabel"]){
    std::cout << "cannot get the edge label of specified edge type" 
              << std::endl;
    return std::pair(EdgeType(), false);
  }
  std::get<1>(edge_type) = edge_type_config["EdgeLabel"].as<EdgeLabelType>();

  if (!edge_type_config["DstLabel"]){
    std::cout << "cannot get the dst label of specified edge type" 
              << std::endl;
    return std::pair(EdgeType(), false);
  }
  std::get<2>(edge_type) = edge_type_config["DstLabel"].as<VertexLabelType>();
  
  return std::pair(edge_type, true);
}

template <typename GraphPatternType,
          typename    DataGraphType>
inline std::pair<gar::LiteralStandAloneInfo<GraphPatternType,
                                               DataGraphType>, bool>
  LiteralStandAloneInfoFromYaml(const YAML::Node& yaml_config){

  using LiteralStandAloneInfoType
      = gar::LiteralStandAloneInfo<GraphPatternType,
                                      DataGraphType>;

  using VertexLabelType = typename GraphPatternType
                                       ::VertexType
                                        ::LabelType;

  using EdgeLabelType = typename GraphPatternType
                                       ::EdgeType
                                      ::LabelType;
                                  
  using VertexAttributeKeyType = typename DataGraphType
                                           ::VertexType
                                     ::AttributeKeyType;

  if (!yaml_config["Type"]){
    std::cout << "literal does not denote type" << std::endl;
    return std::pair(LiteralStandAloneInfoType(), false);
  }
  std::string literal_type_str
       = yaml_config["Type"].as<std::string>();

  if (literal_type_str != "variable_literal"
  //  && literal_type_str != "constant_literal"
   && literal_type_str !=     "edge_literal"){
    std::cout << "unsupported rhs literal type: "
              << literal_type_str << std::endl;
    return std::pair(LiteralStandAloneInfoType(), false);
  }

  enum gar::LiteralType literal_type = gar::LiteralType::kNoneLiteral;
  if (literal_type_str == "constant_literal"){
    literal_type = gar::LiteralType::kConstantLiteral;
  }
  else if (literal_type_str == "variable_literal"){
    literal_type = gar::LiteralType::kVariableLiteral;
  }
  else if (literal_type_str == "edge_literal"){
    literal_type = gar::LiteralType::kEdgeLiteral;
  }

  switch (literal_type) {
  case gar::LiteralType::kAttrValueLiteral: {
    // x.A
    if (!yaml_config["XLabel"]) {
      std::cout << "x label for attribute value literal is not specified" << std::endl;
      return std::pair(LiteralStandAloneInfoType(), false);
    }
    VertexLabelType x_vertex_label = yaml_config["XLabel"].as<VertexLabelType>();
    if (!yaml_config["XAttrKey"]) {
      std::cout << "x attribute key for attribute value literal is not specified" << std::endl;
      return std::pair(LiteralStandAloneInfoType(), false);
    }
    VertexAttributeKeyType x_attr_key = yaml_config["XAttrKey"].as<VertexAttributeKeyType>();
    return std::pair(LiteralStandAloneInfoType(x_vertex_label, x_attr_key), true);
  }
  case gar::LiteralType::kVariableLiteral: {
    // x.A == y.B
    if (!yaml_config["XLabel"]) {
      std::cout << "x label for variable value literal is not specified" << std::endl;
      return std::pair(LiteralStandAloneInfoType(), false);
    }
    VertexLabelType x_vertex_label = yaml_config["XLabel"].as<VertexLabelType>();
    if (!yaml_config["XAttrKey"]) {
      std::cout << "x attribute key for variable value literal is not specified" << std::endl;
      return std::pair(LiteralStandAloneInfoType(), false);
    }
    VertexAttributeKeyType x_attr_key = yaml_config["XAttrKey"].as<VertexAttributeKeyType>();
    
    if (!yaml_config["YLabel"]) {
      std::cout << "y label for variable value literal is not specified" << std::endl;
      return std::pair(LiteralStandAloneInfoType(), false);
    }
    VertexLabelType y_vertex_label = yaml_config["YLabel"].as<VertexLabelType>();
    if (!yaml_config["YAttrKey"]) {
      std::cout << "y attribute key for variable value literal is not specified" << std::endl;
      return std::pair(LiteralStandAloneInfoType(), false);
    }
    VertexAttributeKeyType y_attr_key = yaml_config["YAttrKey"].as<VertexAttributeKeyType>();
    return std::pair(LiteralStandAloneInfoType(x_vertex_label, x_attr_key,
                                               y_vertex_label, y_attr_key), true);
  }
  // case gar::LiteralType::kConstantLiteral:
    // x.A == c
    // if (!yaml_config["XLabel"]) {
    //   std::cout << "x label for constant literal is not specified" << std::endl;
    //   return;
    // }
    // VertexLabelType x_vertex_label = yaml_config["XLabel"].as<VertexLabelType>();

    // if (!yaml_config["XAttrKey"]) {
    //   std::cout << "x attribute key for constant literal is not specified" << std::endl;
    //   return;
    // }
    // std::string x_attr_key = yaml_config["XAttrKey"].as<std::string>();

    // if (!yaml_config["CStr"]) {
    //   std::cout << "c str for constant literal is not specified" << std::endl;
    //   return;
    // }
    // std::string c_str = yaml_config["CStr"].as<std::string>();
    
    // auto [ yaml_config_set_it,
    //        yaml_config_set_ret]
    //      = ctx.yaml_config_set_.emplace(x_vertex_label, x_attr_key,
    //                                     y_vertex_label, y_attr_key);
    // if (!yaml_config_set_ret){
    //   std::cout << "duplicate literal" << std::endl;
    //   return;
    // }
    // break;
    // std::cout << "does not support specified constant literal as rhs literal yet"
    //           << std::endl;
    // return;
  case gar::LiteralType::kEdgeLiteral: {
    // x_id_ - edge_label_ -> y_id_
    if (!yaml_config["XLabel"]) {
      std::cout << "x label for edge literal is not specified" << std::endl;
      return std::pair(LiteralStandAloneInfoType(), false);
    }
    VertexLabelType x_vertex_label = yaml_config["XLabel"].as<VertexLabelType>();

    if (!yaml_config["EdgeLabel"]) {
      std::cout << "edge label for edge literal is not specified" << std::endl;
      return std::pair(LiteralStandAloneInfoType(), false);
    }
    EdgeLabelType edge_label = yaml_config["EdgeLabel"].as<EdgeLabelType>();
    
    if (!yaml_config["YLabel"]) {
      std::cout << "y label for edge literal is not specified" << std::endl;
      return std::pair(LiteralStandAloneInfoType(), false);
    }
    VertexLabelType y_vertex_label = yaml_config["YLabel"].as<VertexLabelType>();

    return std::pair(LiteralStandAloneInfoType(x_vertex_label,
                                               y_vertex_label,
                                                   edge_label), true);
  }
  default:
    std::cout << "unsupported rhs literal type" << std::endl;
    break;
  }
  return std::pair(LiteralStandAloneInfoType(), false);
}

};
      
#endif // UTIL_YAML_CONFIG_H_