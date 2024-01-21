#include <string>
#include <ostream>
#include <iostream>

#include "include/util/log.h"

#include "include/gundam/io/csvgraph.h"
#include "include/gundam/tool/operator/duplicate_vertex.h"
#include "include/gundam/graph_type/large_graph.h"

#include <yaml-cpp/yaml.h>

using VertexIDType = int;
using VertexLabelType = int;
using VertexAttributeKeyType = std::string;
using   EdgeIDType = int;
using   EdgeLabelType = int;
using   EdgeAttributeKeyType = std::string;

using DataGraph =
    GUNDAM::LargeGraph<VertexIDType, VertexLabelType, VertexAttributeKeyType, 
                         EdgeIDType,   EdgeLabelType,   EdgeAttributeKeyType>;

int main(int argc, char **argv) {
  if (argc < 2) {
    util::Error("config yaml path must be given!");
    exit(-1);
  }
  if (argc > 2) {
    util::Error("parameter num is not correct! (0 or 1)");
    exit(-1);
  }

  std::string config_file_path = "duplicate_vertex.yaml";
  if (argc == 2) {
    config_file_path = argv[1];
  }

  YAML::Node config = YAML::LoadFile(config_file_path);

  if (!config["InputGraphPath"]) {
    util::Error("cannot get data graph path");
    return -1;
  }
  YAML::Node input_graph_path_config = config["InputGraphPath"];

  if (!input_graph_path_config["VFile"]) {
    util::Error("cannot get data graph v file");
    return -1;
  }
  const std::string kInputGraphPathVFile 
        = input_graph_path_config["VFile"].as<std::string>();

  if (!input_graph_path_config["EFile"]) {
    util::Error("cannot get data graph e file");
    return -1;
  }
  const std::string kInputGraphPathEFile 
        = input_graph_path_config["EFile"].as<std::string>();

  if (!config["OutputGraphPath"]) {
    util::Error("cannot get data graph path");
    return -1;
  }
  YAML::Node output_graph_path_config = config["OutputGraphPath"];

  if (!output_graph_path_config["VFile"]) {
    util::Error("cannot get data graph v file");
    return -1;
  }
  const std::string kOutputGraphPathVFile 
        = output_graph_path_config["VFile"].as<std::string>();

  if (!output_graph_path_config["EFile"]) {
    util::Error("cannot get data graph e file");
    return -1;
  }
  const std::string kOutputGraphPathEFile 
        = output_graph_path_config["EFile"].as<std::string>();

  if (!config["DuplicateVertexSize"]) {
    util::Error("does not specified duplicated vertex size");
    return -1;
  }
  const size_t kDuplicateVertexSize 
      = config["DuplicateVertexSize"].as<size_t>();

  DataGraph input_graph;

  GUNDAM::ReadCSVGraph(input_graph, kInputGraphPathVFile,
                                    kInputGraphPathEFile);

  using VertexHandleType = typename GUNDAM::VertexHandle<DataGraph>::type;
  std::vector<VertexHandleType> vertex_handle_set;
  vertex_handle_set.reserve(input_graph.CountVertex());
  if (!config["DuplicateVertexLabel"]) {
    for (auto vertex_it = input_graph.VertexBegin();
             !vertex_it.IsDone();
              vertex_it++) {
      vertex_handle_set.emplace_back(vertex_it);
    }
  }
  else {
    auto vertex_label = config["DuplicateVertexLabel"].as<VertexLabelType>();
    for (auto vertex_it = input_graph.VertexBegin();
             !vertex_it.IsDone();
              vertex_it++) {
      if (vertex_it->label() != vertex_label) {
        continue;
      }
      vertex_handle_set.emplace_back(vertex_it);
    }
  }
  if (vertex_handle_set.size() < kDuplicateVertexSize) {
    util::Error("does not have enough vertex to duplicate");
    return 0;
  }

  std::random_shuffle ( vertex_handle_set.begin(), 
                        vertex_handle_set.end() );

  vertex_handle_set.resize(kDuplicateVertexSize);

  auto vertex_id_map = GUNDAM::DuplicateVertex(input_graph, vertex_handle_set);

  GUNDAM::WriteCSVGraph(input_graph, kOutputGraphPathVFile,
                                     kOutputGraphPathEFile);

  if (!config["VertexIdMap"]) {
    for (const auto& [ori_vertex, new_vertex] : vertex_id_map) {
      std::cout << ori_vertex->id() << "," 
                << new_vertex->id() << std::endl;
    }
  }
  else {
    std::ofstream of(config["VertexIdMap"].as<std::string>());
    for (const auto& [ori_vertex, new_vertex] : vertex_id_map) {
      of << ori_vertex->id() << "," 
         << new_vertex->id() << std::endl;
    }
  }

  return 0;
}
