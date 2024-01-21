#ifndef UTIL_LOAD_GRAPH_H_
#define UTIL_LOAD_GRAPH_H_

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "include/gundam/io/csvgraph.h"

namespace util {

template <typename GraphType>
inline std::pair<GraphType, bool> // ok to return in C11
       LoadGraph(const YAML::Node& graph_path_config) {
  std::string v_file, e_file;
  if (graph_path_config["VFile"]) {
    if (!graph_path_config["EFile"]) {
      // configure does not complete
      return std::pair(GraphType(), false);
    }
    v_file = graph_path_config["VFile"].as<std::string>();
    e_file = graph_path_config["EFile"].as<std::string>();         
  }
  else if (graph_path_config["GraphDir"]) {
    if (!graph_path_config["GraphName"]) {
      // configure does not complete
      return std::pair(GraphType(), false);
    }
    v_file = graph_path_config["GraphDir"].as<std::string>()
     + "/" + graph_path_config["GraphName"].as<std::string>()
     + "_v.csv";
    e_file = graph_path_config["GraphDir"].as<std::string>()
     + "/" + graph_path_config["GraphName"].as<std::string>()
     + "_e.csv";
  }
  else {
    return std::pair(GraphType(), false);
  }

  GraphType graph;
  if (GUNDAM::ReadCSVGraph(graph, v_file, 
                                  e_file) < 0) {
    return std::pair(GraphType(), false);
  }
  return std::pair(std::move(graph), true);
}

}; // util

#endif // UTIL_LOAD_GRAPH_H_