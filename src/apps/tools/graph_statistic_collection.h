#ifndef EXAMPLES_ANALYTICAL_APPS_TOOLS_GRAPH_STATISTIC_COLLECTION_H_
#define EXAMPLES_ANALYTICAL_APPS_TOOLS_GRAPH_STATISTIC_COLLECTION_H_

#include <iostream>

#include "gundam/io/csvgraph.h"
#include "gundam/graph_statistics/graph_basic_statistics.h"
#include "gundam/graph_type/large_graph.h"
#include "gundam/graph_type/graph.h"

#include "include/gar/gar.h"
#include "include/util/load_gar_path.h"

namespace tools {

bool GraphStatisticCollection(const ToolInfo tool_info) {

  using VertexIDType    = int;
  using VertexLabelType = int;
  using VertexAttributeKeyType = std::string;
  using   EdgeIDType    = int;
  using   EdgeLabelType = int;
  using   EdgeAttributeKeyType = std::string;

  using DataGraph = GUNDAM::Graph<
                    GUNDAM::SetVertexIDType<int>,
                    GUNDAM::SetVertexAttributeStoreType<GUNDAM::AttributeType::kGrouped>,
                    // GUNDAM::SetVertexAttributeStoreType<GUNDAM::AttributeType::kSeparated>,
                    GUNDAM::SetVertexLabelType<int>,
                    GUNDAM::SetVertexLabelContainerType<GUNDAM::ContainerType::Map>,
                    GUNDAM::SetVertexIDContainerType<GUNDAM::ContainerType::Map>,
                    GUNDAM::SetVertexPtrContainerType<GUNDAM::ContainerType::Map>,
                    GUNDAM::SetEdgeLabelContainerType<GUNDAM::ContainerType::Vector>,
                    GUNDAM::SetVertexAttributeKeyType<std::string>,
                    GUNDAM::SetEdgeIDType<int>,
                    GUNDAM::SetEdgeAttributeStoreType<GUNDAM::AttributeType::kGrouped>,
                    // GUNDAM::SetEdgeAttributeStoreType<GUNDAM::AttributeType::kSeparated>,
                    GUNDAM::SetEdgeLabelType<int>,
                    GUNDAM::SetEdgeAttributeKeyType<std::string>>;

  using Pattern =
      GUNDAM::LargeGraph<VertexIDType, VertexLabelType, VertexAttributeKeyType, 
                           EdgeIDType,   EdgeLabelType,   EdgeAttributeKeyType>;

  using GarType = gar::GraphAssociationRule<Pattern,
                                          DataGraph>;

  DataGraph data_graph;
                                          
  auto res = GUNDAM::ReadCSVGraph(data_graph, 
                           tool_info.v_file_,
                           tool_info.e_file_);

  if (res < 0) {
    util::Error("load data graph fail!");
    return false;
  }

  GUNDAM::GraphBasicStatistics<DataGraph> graph_statistic(data_graph);

  std::string output_info;

  if (tool_info.collect_vertex_label_) {
    output_info += "####################################################################\n" ;
    output_info += "# vertex label statistic collection, form: <vertex_label>, counter #\n" ;
    output_info += "####################################################################\n" ;
    for (const auto& [vertex_label, counter] : graph_statistic.vertex_label_counter()) {
      output_info += "<"  + std::to_string(vertex_label)
                   + "> " + std::to_string(counter)
                   + "\n";
    }
  }

  if (tool_info.collect_vertex_attribute_) {
    output_info += "##########################################################################################\n" ;
    output_info += "# vertex label statistic collection, form: <vertex_label, attr_key, attr_value:type>, counter #\n" ;
    output_info += "##########################################################################################\n" ;
    for (const auto& [vertex_label, attr_set] : graph_statistic.legal_attr_set()) {
      for (const auto& [attr_key, attr_value_counter] : attr_set) {
        for (const auto& [value, counter] : attr_value_counter) {
          output_info += "<"  + std::to_string(vertex_label)
                       + ", " + GUNDAM::ToString(attr_key)
                       + ", " + GUNDAM::ToString(value.first) + ":" + GUNDAM::EnumToString(value.second)
                       + "> " + std::to_string(counter)
                       + "\n";
        }
      }
    }
  }

  if (tool_info.collect_edge_label_) {
    output_info += "################################################################\n" ;
    output_info += "# edge label statistic collection, form: <edge_label>, counter #\n" ;
    output_info += "################################################################\n" ;
    for (const auto& [edge_label, counter] : graph_statistic.edge_label_counter()) {
      output_info += "<"  + std::to_string(edge_label)
                   + "> " + std::to_string(counter)
                   + "\n";
    }
  }

  if (tool_info.collect_edge_type_) {
    output_info += "######################################################################################\n" ;
    output_info += "# edge types statistic collection, form: <src_label, edge_label, dst_label>, counter #\n" ;
    output_info += "######################################################################################\n" ;
    for (auto edge_type_counter_it  = graph_statistic.edge_type_counter_begin();
              edge_type_counter_it != graph_statistic.edge_type_counter_end();
              edge_type_counter_it++) {
      output_info += "<"  + std::to_string(std::get<0>(edge_type_counter_it->first))
                   + ","  + std::to_string(std::get<1>(edge_type_counter_it->first))
                   + ","  + std::to_string(std::get<2>(edge_type_counter_it->first)) 
                   + "> " + std::to_string(edge_type_counter_it->second)
                   + "\n";
    }
  }

  if (tool_info.output_file_ == "")  {
    // export to standard output
    std::cout << output_info;
    return true;
  }

  std::ofstream output_file(tool_info.output_file_);
  output_file << output_info;

  return true;
}

};


#endif // EXAMPLES_ANALYTICAL_APPS_TOOLS_GRAPH_STATISTIC_COLLECTION_H_