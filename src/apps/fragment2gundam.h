#ifndef _GRAPE_FRAGMENT2GUNDAM_H
#define _GRAPE_FRAGMENT2GUNDAM_H
#include <iostream>
#include <string>
#include <vector>
namespace grape {
template <bool store_outervertex = true, class Fragment, class GUNDAMGraph>
void Fragment2GUNDAMGraph(Fragment& fragment, GUNDAMGraph& gundam_graph) {
  auto vertex_set = fragment.Vertices();
  using VertexIDType = typename GUNDAMGraph::VertexType::IDType;
  using VertexLabelType = typename GUNDAMGraph::VertexType::LabelType;
  using EdgeIDType = typename GUNDAMGraph::EdgeType::IDType;
  using EdgeLabelType = typename GUNDAMGraph::EdgeType::LabelType;
  using AttributeContainer = std::vector<std::string>;

  for (auto& vertex : vertex_set) {
    if (!store_outervertex && fragment.IsOuterVertex(vertex)) continue;
    VertexIDType vertex_id = fragment.GetId(vertex);
    auto& vertex_data = fragment.GetData(vertex);
    VertexLabelType vertex_label = vertex_data.label_;

    auto [vertex_ptr, add_flag] =
        gundam_graph.AddVertex(vertex_id, vertex_label);
    if (add_flag) {
      if (gundam_graph.vertex_has_attribute) {
        AttributeContainer vertex_attribute = vertex_data.attributes_;
        for (size_t i = 0; i < vertex_attribute.size(); i++) {
          vertex_ptr->AddAttribute(i, vertex_attribute[i]);
        }
      }
    }
  }
  EdgeIDType edge_id = 0;
  for (auto& vertex : vertex_set) {
    if (!fragment.IsInnerVertex(vertex)) continue;
    VertexIDType vertex_id = fragment.GetId(vertex);
    auto out_adj_list = fragment.GetOutgoingAdjList(vertex);
    for (const auto& single_adj : out_adj_list) {
      if (!store_outervertex && fragment.IsOuterVertex(single_adj.neighbor))
        continue;
      VertexIDType dst_id = fragment.GetId(single_adj.neighbor);
      EdgeLabelType edge_label = single_adj.data.label_;
      auto [edge_ptr, add_flag] =
          gundam_graph.AddEdge(vertex_id, dst_id, edge_label, ++edge_id);
      if (add_flag) {
      }
    }
    auto inner_adj_list = fragment.GetIncomingAdjList(vertex);
    if (store_outervertex) {
      for (const auto& single_adj : inner_adj_list) {
        if (!fragment.IsOuterVertex(single_adj.neighbor)) continue;
        VertexIDType dst_id = fragment.GetId(single_adj.neighbor);
        EdgeLabelType edge_label = single_adj.data.label_;
        auto [edge_ptr, add_flag] =
            gundam_graph.AddEdge(dst_id, vertex_id, edge_label, ++edge_id);
        if (add_flag) {
        }
      }
    }
  }
}
}  // namespace grape

#endif