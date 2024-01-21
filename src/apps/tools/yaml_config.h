#ifndef EXAMPLES_ANALYTICAL_APPS_TOOLS_YAML_CONFIG_H_
#define EXAMPLES_ANALYTICAL_APPS_TOOLS_YAML_CONFIG_H_

#include <yaml-cpp/yaml.h>

#include <memory>
#include <vector>

#include <assert.h>

#include "src/apps/tools/tools_info.h"
#include "include/util/log.h"
#include "include/util/yaml_config.h"

namespace tools {

bool ReadGraphPackageInfoYAML(const YAML::Node &graph_node,
                                      ToolInfo &graph_info);

bool ReadGraphPackageInfoYAML(const YAML::Node &graph_node,
                                      ToolInfo &graph_info) {
  using namespace YAML;

  graph_info.Reset();

  if (!graph_node["Action"]) {
    util::Error("cannot find action!");
    return false;
  }

  const std::string kType = graph_node["Action"].as<std::string>();

  if (kType == "GarSetToGar") {
    graph_info.action_ = tools::ActionType::kGarSetToGar; 

    if (!graph_node["GarSetDir"]) {
      util::Error("cannot find gar set dir in GarSetToGar");
      return false;
    }
    graph_info.gar_set_dir_ 
             = graph_node["GarSetDir"].as<std::string>();

    graph_info.gar_size_in_each_set_file_ = 0;
    if (graph_node["GarSizeInEachSetFile"]) {
      graph_info.gar_size_in_each_set_file_ = graph_node["GarSizeInEachSetFile"].as<int>();
    }

    if (!graph_node["GarSetPrefix"]) {
      assert(graph_info.gar_set_prefix_set_.empty());
      // does not specify prefix, convert all legal gar set 
      // under this dir
      graph_info.gar_set_prefix_set_.emplace_back("");
      return true;
    }

    YAML::Node gar_set_prefix_config = graph_node["GarSetPrefix"];
    if (!gar_set_prefix_config.IsSequence()) {
      assert(graph_info.gar_set_prefix_set_.empty());
      graph_info.gar_set_prefix_set_.emplace_back(gar_set_prefix_config.as<std::string>());
      return true;
    }

    for (int i = 0; i < gar_set_prefix_config.size(); i++) {
      graph_info.gar_set_prefix_set_.emplace_back(gar_set_prefix_config[i].as<std::string>());
    }
    return true;
  }
  // else (kType == "GarToGarSet") {
  // }
  else if (kType == "GarInfoCollection") {
    graph_info.action_ = tools::ActionType::kGarInfoCollection;

    if (!graph_node["GarSetDir"]) {
      util::Error("cannot find gar set dir in GarInfoCollection");
      return false;
    }
    graph_info.gar_set_dir_ 
             = graph_node["GarSetDir"].as<std::string>();

    if (!graph_node["OutputFile"]) {
      util::Error("does not specify output file in GarInfoCollection");
      return false;
    }
    graph_info.output_file_
             = graph_node["OutputFile"].as<std::string>();

    if (!graph_node["InfoToCollect"]) {
      util::Error("does not specify infomation to collect in GarInfoCollection");
      return false;
    }

    YAML::Node gar_info_to_collect_config = graph_node["InfoToCollect"];
    for (int info_idx = 0; 
             info_idx < gar_info_to_collect_config.size(); 
             info_idx++) {
      const std::string kInfoType = gar_info_to_collect_config[info_idx].as<std::string>();
      if (kInfoType == "vertex_num") {
        if (graph_info.collect_vertex_num_)  {
          util::Error("re-specify infomation to collect: vertex_num");
        }
        graph_info.collect_vertex_num_ = true;
      }
      else if (kInfoType == "edge_num") {
        if (graph_info.collect_edge_num_)  {
          util::Error("re-specify infomation to collect: edge_num");
        }
        graph_info.collect_edge_num_ = true;
      }
      else if (kInfoType == "is_link") {
        if (graph_info.collect_is_link_)  {
          util::Error("re-specify infomation to collect: is_link");
        }
        graph_info.collect_is_link_ = true;
      }
      else if (kInfoType == "is_star") {
        if (graph_info.collect_is_star_)  {
          util::Error("re-specify infomation to collect: is_star");
        }
        graph_info.collect_is_star_ = true;
      }
      else if (kInfoType == "x_literal_num") {
        if (graph_info.x_literal_num_)  {
          util::Error("re-specify infomation to collect: x_literal_num");
        }
        graph_info.x_literal_num_ = true;
      }
      else if (kInfoType == "y_literal_num") {
        if (graph_info.y_literal_num_)  {
          util::Error("re-specify infomation to collect: y_literal_num");
        }
        graph_info.y_literal_num_ = true;
      }
      else {
        util::Error("unknown infomation to collect: " + kInfoType);
        return false;
      }
    }

    if (!graph_node["GarSetPrefix"]) {
      assert(graph_info.gar_set_prefix_set_.empty());
      // does not specify prefix, convert all legal gar set 
      // under this dir
      graph_info.gar_set_prefix_set_.emplace_back("");
      return true;
    }

    YAML::Node gar_set_prefix_config = graph_node["GarSetPrefix"];
    if (!gar_set_prefix_config.IsSequence()) {
      assert(graph_info.gar_set_prefix_set_.empty());
      graph_info.gar_set_prefix_set_.emplace_back(gar_set_prefix_config.as<std::string>());
      return true;
    }

    for (int i = 0; i < gar_set_prefix_config.size(); i++) {
      graph_info.gar_set_prefix_set_.emplace_back(gar_set_prefix_config[i].as<std::string>());
    }
    return true;
  }
  else if (kType == "GraphStatisticCollection") {
    graph_info.action_ = tools::ActionType::kGraphStatisticCollection;

    // collect output file, if not specified, output to std::out
    graph_info.output_file_ = "";
    if (graph_node["OutputFile"]) {
      graph_info.output_file_ = graph_node["OutputFile"].as<std::string>();
    }

    // collect data path
    if (!graph_node["GraphPath"]) {
      util::Error("cannot get graph path in GraphStatisticCollection!");
      return false;
    }

    auto [ data_graph_string, 
           data_graph_ret ] = util::GetDataGraphInfoFromYaml(graph_node["GraphPath"]);

    if (!data_graph_ret) {
      util::Error("load data graph error in GraphStatisticCollection!");
      return false;
    }

    graph_info.v_file_ = std::get<0>(data_graph_string);
    graph_info.e_file_ = std::get<1>(data_graph_string);

    if (!graph_node["StatisticToCollect"]) {
      util::Error("does not specify statistic to collect in GraphStatisticCollection");
      return false;
    }

    YAML::Node graph_statistic_to_collect_config = graph_node["StatisticToCollect"];
    for (int statistic_idx = 0;
             statistic_idx < graph_statistic_to_collect_config.size();
             statistic_idx++) {
      const std::string kInfoType = graph_statistic_to_collect_config[statistic_idx].as<std::string>();

      if (kInfoType == "vertex_label") {
        if (graph_info.collect_vertex_label_)  {
          util::Error("re-specify statistic to collect: vertex_label");
        }
        graph_info.collect_vertex_label_ = true;
      }
      else if (kInfoType == "edge_label") {
        if (graph_info.collect_edge_label_)  {
          util::Error("re-specify statistic to collect: edge_label");
        }
        graph_info.collect_edge_label_ = true;
      }
      else if (kInfoType == "edge_type") {
        if (graph_info.collect_edge_type_)  {
          util::Error("re-specify statistic to collect: edge_type");
        }
        graph_info.collect_edge_type_ = true;
      }
      else if (kInfoType == "vertex_attribute") {
        if (graph_info.collect_vertex_attribute_)  {
          util::Error("re-specify statistic to collect: vertex_attribute");
        }
        graph_info.collect_vertex_attribute_ = true;
      }
      else {
        util::Error("unknown infomation to collect: " + kInfoType);
        return false;
      }
    }
  }
  else {
    util::Error("unknown action: " + kType);
    // throw InvalidNode("Action");
    return false;
  }
  return true;
}

}  // namespace tools

#endif // EXAMPLES_ANALYTICAL_APPS_TOOLS_YAML_CONFIG_H_