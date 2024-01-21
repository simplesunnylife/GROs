#ifndef _GRAPE_FRAGMENTGRAPH_H
#define _GRAPE_FRAGMENTGRAPH_H
#include <type_traits>

#include "gundam/component/iterator2.h"
#include "gundam/component/pointer.h"
namespace grape {
template <class Fragment>
class FragmentGraph {
 public:
  FragmentGraph(Fragment &frag) : fragment(frag) {}
  using GraphType = FragmentGraph;
  using FragmentVertexType = typename Fragment::vertex_t;
  using FragmentEdgeType = typename Fragment::nbr_t;
  using VertexIDType = typename Fragment::oid_t;
  using EdgeIDType = typename Fragment::oid_t;
  using EdgeLabelType = typename Fragment::edata_t::LabelType;
  template <bool is_const>
  class _Vertex;
  template <bool is_const>
  class _Edge;
  using Vertex = _Vertex<false>;

  using ConstVertex = _Vertex<true>;

  using Edge = _Edge<false>;

  using ConstEdge = _Edge<true>;
  using VertexPtr = GUNDAM::GPointer<false, Vertex, ConstVertex>;

  using VertexConstPtr = GUNDAM::GPointer<true, Vertex, ConstVertex>;

  using EdgePtr = GUNDAM::GPointer<false, Edge, ConstEdge>;

  using EdgeConstPtr = GUNDAM::GPointer<true, Edge, ConstEdge>;
  using VertexIterator =
      GUNDAM::GIterator2<false, FragmentGraph, FragmentVertexType, Vertex,
                         VertexPtr>;

  using VertexConstIterator =
      GUNDAM::GIterator2<true, const FragmentGraph, FragmentVertexType,
                         ConstVertex, VertexConstPtr>;

  using EdgeIterator =
      GUNDAM::GIterator2<false, FragmentGraph,
                         typename Fragment::adj_list_t::iterator, Edge,
                         EdgePtr>;

  using EdgeConstIterator = GUNDAM::GIterator2<
      true, const FragmentGraph,
      typename std::remove_const_t<
          typename Fragment::const_adj_list_t::const_iterator>,
      ConstEdge, EdgeConstPtr>;

  template <bool is_const>
  class _Vertex {
   public:
    using GraphType = FragmentGraph;
    using IDType = VertexIDType;
    using LabelType = typename Fragment::vdata_t::LabelType;
    friend typename GraphType::EdgeIterator;
    friend typename GraphType::EdgeConstIterator;
    _Vertex() : graph_(nullptr){};
    _Vertex(typename std::conditional<is_const, const GraphType *,
                                      GraphType *>::type graph,
            const FragmentVertexType vertex)
        : graph_(graph), vertex_(vertex) {}

    const IDType id() const { return graph_->fragment.GetId(vertex_); }
    const LabelType &label() const {
      return graph_->fragment.GetData(vertex_).label_;
    }
    EdgeIterator OutEdgeBegin() {
      assert(HasValue());
      auto &&data = graph_->fragment.GetOutgoingAdjList(vertex_);
      auto it = EdgeIterator(graph_, data.begin(), data.end());
      return it;
    }
    EdgeConstIterator OutEdgeCBegin() {
      assert(HasValue());
      auto &&data = graph_->fragment.GetOutgoingAdjList(vertex_);
      auto it = EdgeConstIterator(graph_, data.begin(), data.end());
      return it;
    }
    EdgeIterator InEdgeBegin() {
      assert(HasValue());
      auto &&data = graph_->fragment.GetIncomingAdjList(vertex_);
      auto it = EdgeIterator(graph_, data.begin(), data.end());
      return it;
    }
    EdgeConstIterator InEdgeCBegin() {
      assert(HasValue());
      auto &&data = graph_->fragment.GetIncomingAdjList(vertex_);
      auto it = EdgeConstIterator(graph_, data.begin(), data.end());
      return it;
    }
    bool operator==(const _Vertex &b) const {
      if (!graph_) {
        return !b.graph_;
      } else if (!b.graph_) {
        return false;
      } else {
        return this->id() == b.id();
      }
    }
    bool operator<(const _Vertex &b) const {
      if (!graph_ || !b.graph_) {
        return false;
      } else {
        return this->id() < b.id();
      }
    }

   public:
    void Reset() { graph_ = nullptr; }

    bool HasValue() const { return (graph_ != nullptr); }
    typename std::conditional<is_const, const GraphType *, GraphType *>::type
        graph_;
    FragmentVertexType vertex_;
  };

  template <bool is_const>
  class _Edge {
   public:
    using GraphType = FragmentGraph;
    using IDType = EdgeIDType;
    _Edge() : graph_(nullptr){};
    _Edge(typename std::conditional<is_const, const GraphType *,
                                    GraphType *>::type graph,
          const FragmentEdgeType &edge)
        : graph_(graph), edge_(&edge) {}

   public:
    const EdgeIDType id() const { return edge_->data.id_; }
    const EdgeLabelType &label() const { return edge_->data.label_; }
    const VertexIDType src_id() const {
      return graph_->fragment.GetId(edge_->neighbor);
    }
    const VertexIDType dst_id() const {
      return graph_->fragment.GetId(edge_->neighbor);
    }
    VertexPtr src_handle() { return VertexPtr(Vertex(graph_, edge_->neighbor)); }
    VertexConstPtr const_src_handle() const {
      return VertexConstPtr(ConstVertex(graph_, edge_->neighbor));
    }
    VertexPtr dst_handle() { return VertexPtr(Vertex(graph_, edge_->neighbor)); }
    VertexConstPtr const_dst_handle() const {
      return VertexConstPtr(ConstVertex(graph_, edge_->neighbor));
    }
    bool operator==(const _Edge &b) const {
      if (!graph_) {
        return !b.graph_;
      } else if (!b.graph_) {
        return false;
      } else {
        return this->id() == b.id();
      }
    }

    bool operator<(const _Edge &b) const {
      if (!graph_ || !b.graph_) {
        return false;
      } else {
        return this->id() == b.id();
      }
    }

   private:
    void Reset() { graph_ = nullptr; }

    bool HasValue() const { return (graph_ != nullptr); }

    typename std::conditional<is_const, const GraphType *, GraphType *>::type
        graph_;
    const FragmentEdgeType *edge_;
  };

 public:
  VertexPtr FindVertex(const VertexIDType &id) {
    FragmentVertexType vertex;
    fragment.GetVertex(id, vertex);
    return VertexPtr(Vertex(this, vertex));
  }

  VertexConstPtr FindConstVertex(const VertexIDType &id) const {
    FragmentVertexType vertex;
    fragment.GetVertex(id, vertex);
    return VertexConstPtr(ConstVertex(this, vertex));
  }

  VertexIterator VertexBegin() {
    return VertexIterator(this, this->fragment.Vertices().begin(),
                          this->fragment.Vertices().end());
  }

  VertexConstIterator VertexCBegin() const {
    return VertexConstIterator(this, this->fragment.Vertices().begin(),
                               this->fragment.Vertices().end());
  }

 private:
  Fragment &fragment;
};
}  // namespace grape

#endif