#ifndef _GRAPE_FRAGMENTGRAPHWITHINDEX_H
#define _GRAPE_FRAGMENTGRAPHWITHINDEX_H
#include <type_traits>

#include "gundam/component/iterator2.h"
#include "gundam/component/pointer.h"

#include "gundam/type_getter/edge_handle.h"
#include "gundam/type_getter/vertex_handle.h"

namespace grape {
template <class Fragment>
class FragmentGraphWithIndex {
 public:
  FragmentGraphWithIndex() : fragment(nullptr) {}
  FragmentGraphWithIndex(Fragment &frag) : fragment(&frag) {}
  using GraphType = FragmentGraphWithIndex;
  using FragmentVertexType = typename Fragment::vertex_t;
  using FragmentEdgeType = typename Fragment::nbr_t;
  using VertexIDType = typename Fragment::vid_t;
  using VertexLabelType = typename Fragment::vdata_t::LabelType;
  using EdgeIDType = typename Fragment::oid_t;
  using EdgeLabelType = typename Fragment::edata_t::LabelType;
  using IndexEdgeContainer = std::vector<const FragmentEdgeType *>;
  using LabelIndexEdgeContainer = std::map<EdgeLabelType, IndexEdgeContainer>;
  using VertexLabelIndexEdgeContainer =
      std::map<FragmentVertexType, LabelIndexEdgeContainer>;
  using IndexVertexContainer = std::vector<FragmentVertexType>;
  using LabelIndexVertexContainer =
      std::map<EdgeLabelType, IndexVertexContainer>;

  using VertexLabelEdgeLabelIndexEdgeContainer =
      std::map<VertexLabelType, LabelIndexEdgeContainer>;

  using VertexLabelEdgeLabelIndexVertexContainer =
      std::map<VertexLabelType, LabelIndexVertexContainer>;

  template <bool is_const>
  class _VertexWithIndex;

  template <bool is_const>
  class _Edge;

  using Vertex = _VertexWithIndex<false>;

  using ConstVertex = _VertexWithIndex<true>;

  using Edge = _Edge<false>;

  using ConstEdge = _Edge<true>;

 private:
  friend class GUNDAM::VertexHandle<      FragmentGraphWithIndex>;
  friend class GUNDAM::VertexHandle<const FragmentGraphWithIndex>;

  friend class GUNDAM::EdgeHandle<      FragmentGraphWithIndex>;
  friend class GUNDAM::EdgeHandle<const FragmentGraphWithIndex>;

  using VertexPtr = GUNDAM::GPointer<false, Vertex, ConstVertex>;

  using VertexConstPtr = GUNDAM::GPointer<true, Vertex, ConstVertex>;

  using EdgePtr = GUNDAM::GPointer<false, Edge, ConstEdge>;

  using EdgeConstPtr = GUNDAM::GPointer<true, Edge, ConstEdge>;

 public:
  using VertexType = Vertex;

  using EdgeType = Edge;

  using VertexIterator =
      GUNDAM::GIterator2<false, FragmentGraphWithIndex, FragmentVertexType,
                         Vertex, VertexPtr>;

  using VertexConstIterator =
      GUNDAM::GIterator2<true, const FragmentGraphWithIndex, FragmentVertexType,
                         ConstVertex, VertexConstPtr>;

  using EdgeIterator =
      GUNDAM::GIterator2<false, FragmentGraphWithIndex,
                         typename Fragment::adj_list_t::iterator, Edge,
                         EdgePtr>;

  using EdgeConstIterator = GUNDAM::GIterator2<
      true, const FragmentGraphWithIndex,
      typename std::remove_const_t<
          typename Fragment::const_adj_list_t::const_iterator>,
      ConstEdge, EdgeConstPtr>;
  using EdgeIndexIterator =
      GUNDAM::GIterator3<false, FragmentGraphWithIndex,
                         typename IndexEdgeContainer::iterator, Edge, EdgePtr>;
  using EdgeIndexConstIterator =
      GUNDAM::GIterator3<true, const FragmentGraphWithIndex,
                         typename IndexEdgeContainer::const_iterator, ConstEdge,
                         EdgeConstPtr>;
  using VertexIndexIterator =
      GUNDAM::GIterator2<false, FragmentGraphWithIndex,
                         typename IndexVertexContainer::iterator, Vertex,
                         VertexPtr>;
  using VertexIndexConstIterator =
      GUNDAM::GIterator2<true, const FragmentGraphWithIndex,
                         typename IndexVertexContainer::const_iterator,
                         ConstVertex, VertexConstPtr>;
  template <bool is_const>
  class _VertexWithIndex {
   public:
    using GraphType = FragmentGraphWithIndex;
    using IDType = VertexIDType;
    using LabelType = typename Fragment::vdata_t::LabelType;
    friend typename GraphType::EdgeIterator;
    friend typename GraphType::EdgeConstIterator;
    friend typename GraphType::EdgeIndexIterator;
    friend typename GraphType::EdgeIndexConstIterator;
    friend class FragmentGraphWithIndex;
    friend typename GraphType::VertexPtr;
    _VertexWithIndex() : graph_(nullptr){};
    _VertexWithIndex(typename std::conditional<is_const, const GraphType *,
                                               GraphType *>::type graph,
                     const FragmentVertexType vertex)
        : graph_(graph), vertex_(vertex) {}

    const IDType id() const { return graph_->fragment->GetId(vertex_); }
    const LabelType &label() const {
      return graph_->fragment->GetData(vertex_).label_;
    }
    EdgeIterator OutEdgeBegin() {
      assert(HasValue());
      auto &&data = graph_->fragment->GetOutgoingAdjList(vertex_);
      auto it = EdgeIterator(graph_, data.begin(), data.end());
      return it;
    }
    EdgeIndexIterator OutEdgeBegin(const EdgeLabelType edge_label) {
      assert(HasValue());
      auto it = graph_->out_edge_label_index_.find(this->vertex_);
      if (it == graph_->out_edge_label_index_.end()) {
        return EdgeIndexIterator();
      }
      auto &&data1 = it->second;
      auto it1 = data1.find(edge_label);
      if (it1 == data1.end()) {
        return EdgeIndexIterator();
      }
      auto &&data = it1->second;
      auto it2 = EdgeIndexIterator(graph_, data.begin(), data.end());
      return it2;
    }
    EdgeIndexIterator OutEdgeBegin(const EdgeLabelType edge_label,
                                   const VertexPtr vertex_ptr) {
      assert(HasValue());
      auto it = graph_->out_vertex_edge_label_index_.find(this->vertex_);
      if (it == graph_->out_vertex_edge_label_index_.end()) {
        return EdgeIndexIterator();
      }
      auto it1 = it->second.find(vertex_ptr->vertex());
      if (it1 == it->second.end()) {
        return EdgeIndexIterator();
      }
      auto it2 = it1->second.find(edge_label);
      if (it2 == it1->second.end()) {
        return EdgeIndexIterator();
      }
      auto &&data = it2->second;
      return EdgeIndexIterator(graph_, data.begin(), data.end());
    }
    EdgeConstIterator OutEdgeCBegin() const {
      assert(HasValue());
      auto &&data = graph_->fragment->GetOutgoingAdjList(vertex_);
      auto it = EdgeConstIterator(graph_, data.begin(), data.end());
      return it;
    }
    EdgeIndexConstIterator OutEdgeCBegin(const EdgeLabelType edge_label) const {
      assert(HasValue());
      auto it = graph_->out_edge_label_index_.find(this->vertex_);
      if (it == graph_->out_edge_label_index_.end()) {
        return EdgeIndexConstIterator();
      }
      auto &&data1 = it->second;
      auto it1 = data1.find(edge_label);
      if (it1 == data1.end()) {
        return EdgeIndexConstIterator();
      }
      auto &&data = it1->second;
      auto it2 = EdgeIndexConstIterator(graph_, data.cbegin(), data.cend());
      return it2;
    }
    EdgeIndexConstIterator OutEdgeCBegin(
        const EdgeLabelType edge_label, const VertexConstPtr vertex_ptr) const {
      assert(HasValue());
      auto it = graph_->out_vertex_edge_label_index_.find(this->vertex_);
      if (it == graph_->out_vertex_edge_label_index_.end()) {
        return EdgeIndexConstIterator();
      }
      auto it1 = it->second.find(vertex_ptr->vertex());
      if (it1 == it->second.end()) {
        return EdgeIndexConstIterator();
      }
      auto it2 = it1->second.find(edge_label);
      if (it2 == it1->second.end()) {
        return EdgeIndexConstIterator();
      }
      auto &&data = it2->second;
      return EdgeIndexConstIterator(graph_, data.cbegin(), data.cend());
    }
    EdgeIterator InEdgeBegin() {
      assert(HasValue());
      auto &&data = graph_->fragment->GetIncomingAdjList(vertex_);
      auto it = EdgeIterator(graph_, data.begin(), data.end());
      return it;
    }
    EdgeIndexIterator InEdgeBegin(const EdgeLabelType edge_label) {
      assert(HasValue());
      auto it = graph_->in_edge_label_index_.find(this->vertex_);
      if (it == graph_->in_edge_label_index_.end()) {
        return EdgeIndexIterator();
      }
      auto &&data1 = it->second;
      auto it1 = data1.find(edge_label);
      if (it1 == data1.end()) {
        return EdgeIndexIterator();
      }
      auto &&data = it1->second;
      auto it2 = EdgeIndexIterator(graph_, data.begin(), data.end());
      return it2;
    }
    EdgeIndexIterator InEdgeBegin(const EdgeLabelType edge_label,
                                  const VertexPtr vertex_ptr) {
      assert(HasValue());
      auto it = graph_->in_vertex_edge_label_index_.find(this->vertex_);
      if (it == graph_->in_vertex_edge_label_index_.end()) {
        return EdgeIndexIterator();
      }
      auto it1 = it->second.find(vertex_ptr->vertex());
      if (it1 == it->second.end()) {
        return EdgeIndexIterator();
      }
      auto it2 = it1->second.find(edge_label);
      if (it2 == it1->second.end()) {
        return EdgeIndexIterator();
      }
      auto &&data = it2->second;
      return EdgeIndexIterator(graph_, data.begin(), data.end());
    }

    EdgeConstIterator InEdgeCBegin() const {
      assert(HasValue());
      auto &&data = graph_->fragment->GetIncomingAdjList(vertex_);
      auto it = EdgeConstIterator(graph_, data.begin(), data.end());
      return it;
    }
    EdgeIndexConstIterator InEdgeCBegin(const EdgeLabelType edge_label) {
      assert(HasValue());
      auto it = graph_->in_edge_label_index_.find(this->vertex_);
      if (it == graph_->in_edge_label_index_.end()) {
        return EdgeIndexConstIterator();
      }
      auto &&data1 = it->second;
      auto it1 = data1.find(edge_label);
      if (it1 == data1.end()) {
        return EdgeIndexConstIterator();
      }
      auto &&data = it1->second;
      auto it2 = EdgeIndexIterator(graph_, data.begin(), data.end());
      return it2;
    }
    EdgeIndexConstIterator InEdgeCBegin(const EdgeLabelType edge_label,
                                        const VertexConstPtr vertex_ptr) const {
      assert(HasValue());
      auto it = graph_->in_vertex_edge_label_index_.find(this->vertex_);
      if (it == graph_->in_vertex_edge_label_index_.end()) {
        return EdgeIndexConstIterator();
      }
      auto it1 = it->second.find(vertex_ptr->vertex());
      if (it1 == it->second.end()) {
        return EdgeIndexConstIterator();
      }
      auto it2 = it1->second.find(edge_label);
      if (it2 == it1->second.end()) {
        return EdgeIndexConstIterator();
      }
      auto &&data = it2->second;
      return EdgeIndexConstIterator(graph_, data.cbegin(), data.cend());
    }

    VertexIndexIterator OutVertexBegin() {
      auto it = graph_->out_index_vertex_.find(this->vertex_);
      if (it == graph_->out_index_vertex_.end()) {
        return VertexIndexIterator();
      }
      auto &&data = it->second;
      return VertexIndexIterator(graph_, data.begin(), data.end());
    }
    VertexIndexIterator OutVertexBegin(const EdgeLabelType edge_label) {
      auto it = graph_->out_edge_label_index_vertex_.find(this->vertex_);
      if (it == graph_->out_edge_label_index_vertex_.end()) {
        return VertexIndexIterator();
      }
      auto it1 = it->second.find(edge_label);
      if (it1 == it->second.end()) {
        return VertexIndexIterator();
      }
      auto &&data = it1->second;
      return VertexIndexIterator(graph_, data.begin(), data.end());
    }

    VertexIndexConstIterator OutVertexCBegin() const {
      auto it = graph_->out_index_vertex_.find(this->vertex_);
      if (it == graph_->out_index_vertex_.end()) {
        return VertexIndexConstIterator();
      }
      auto &&data = it->second;
      return VertexIndexConstIterator(graph_, data.cbegin(), data.cend());
    }
    VertexIndexConstIterator OutVertexCBegin(
        const EdgeLabelType edge_label) const {
      auto it = graph_->out_edge_label_index_vertex_.find(this->vertex_);
      if (it == graph_->out_edge_label_index_vertex_.end()) {
        return VertexIndexConstIterator();
      }
      auto it1 = it->second.find(edge_label);
      if (it1 == it->second.end()) {
        return VertexIndexConstIterator();
      }
      auto &&data = it1->second;
      return VertexIndexConstIterator(graph_, data.cbegin(), data.cend());
    }
    VertexIndexIterator InVertexBegin() {
      auto it = graph_->in_index_vertex_.find(this->vertex_);
      if (it == graph_->in_index_vertex_.end()) {
        return VertexIndexIterator();
      }
      auto &&data = it->second;
      return VertexIndexIterator(graph_, data.begin(), data.end());
    }
    VertexIndexIterator InVertexBegin(const EdgeLabelType edge_label) {
      auto it = graph_->in_edge_label_index_vertex_.find(this->vertex_);
      if (it == graph_->in_edge_label_index_vertex_.end()) {
        return VertexIndexIterator();
      }
      auto it1 = it->second.find(edge_label);
      if (it1 == it->second.end()) {
        return VertexIndexIterator();
      }
      auto &&data = it1->second;
      return VertexIndexIterator(graph_, data.begin(), data.end());
    }
    VertexIndexConstIterator InVertexCBegin() const {
      auto it = graph_->in_index_vertex_.find(this->vertex_);
      if (it == graph_->in_index_vertex_.end()) {
        return VertexIndexConstIterator();
      }
      auto &&data = it->second;
      return VertexIndexConstIterator(graph_, data.cbegin(), data.cend());
    }
    VertexIndexConstIterator InVertexCBegin(
        const EdgeLabelType edge_label) const {
      auto it = graph_->in_edge_label_index_vertex_.find(this->vertex_);
      if (it == graph_->in_edge_label_index_vertex_.end()) {
        return VertexIndexConstIterator();
      }
      auto it1 = it->second.find(edge_label);
      if (it1 == it->second.end()) {
        return VertexIndexConstIterator();
      }
      auto &&data = it1->second;
      return VertexIndexConstIterator(graph_, data.cbegin(), data.cend());
    }

    size_t CountOutEdge() const {
      return this->graph_->fragment->GetLocalOutDegree(this->vertex_);
    }
    size_t CountOutEdge(const EdgeLabelType edge_label) const {
      auto it = this->graph_->out_edge_label_index_.find(this->vertex_);
      if (it == this->graph_->out_edge_label_index_.end()) return 0;
      auto it1 = it->second.find(edge_label);
      if (it1 == it->second.end()) return 0;
      return it1->second.size();
    }
    size_t CountOutEdge(const EdgeLabelType edge_label,
                        VertexConstPtr vertex_ptr) const {
      auto it = this->graph_->out_vertex_edge_label_index_.find(this->vertex_);
      if (it == this->graph_->out_vertex_edge_label_index_.end()) {
        return 0;
      }
      auto it1 = it->second.find(vertex_ptr->vertex());
      if (it1 == it->second.end()) {
        return 0;
      }
      auto it2 = it1->second.find(edge_label);
      if (it2 == it1->second.end()) {
        return 0;
      }
      return it2->second.size();
    }
    size_t CountOutEdge(const EdgeLabelType edge_label,
                        const VertexLabelType vertex_label) const {
      auto it =
          this->graph_->out_vertex_label_edge_label_index_.find(this->vertex_);
      if (it == this->graph_->out_vertex_label_edge_label_index_.end()) {
        return 0;
      }
      auto it1 = it->second.find(vertex_label);
      if (it1 == it->second.end()) {
        return 0;
      }
      auto it2 = it1->second.find(edge_label);
      if (it2 == it1->second.end()) {
        return 0;
      }
      return it2->second.size();
    }
    size_t CountInEdge() const {
      return this->graph_->fragment->GetLocalInDegree(this->vertex_);
    }
    size_t CountInEdge(const EdgeLabelType edge_label) const {
      auto it = this->graph_->in_edge_label_index_.find(this->vertex_);
      if (it == this->graph_->in_edge_label_index_.end()) return 0;
      auto it1 = it->second.find(edge_label);
      if (it1 == it->second.end()) return 0;
      return it1->second.size();
    }
    size_t CountInEdge(const EdgeLabelType edge_label,
                       VertexConstPtr vertex_ptr) const {
      auto it = this->graph_->in_vertex_edge_label_index_.find(this->vertex_);
      if (it == this->graph_->in_vertex_edge_label_index_.end()) {
        return 0;
      }
      auto it1 = it->second.find(vertex_ptr->vertex());
      if (it1 == it->second.end()) {
        return 0;
      }
      auto it2 = it1->second.find(edge_label);
      if (it2 == it1->second.end()) {
        return 0;
      }
      return it2->second.size();
    }
    size_t CountInEdge(const EdgeLabelType edge_label,
                       const VertexLabelType vertex_label) const {
      auto it =
          this->graph_->in_vertex_label_edge_label_index_.find(this->vertex_);
      if (it == this->graph_->in_vertex_label_edge_label_index_.end()) {
        return 0;
      }
      auto it1 = it->second.find(vertex_label);
      if (it1 == it->second.end()) {
        return 0;
      }
      auto it2 = it1->second.find(edge_label);
      if (it2 == it1->second.end()) {
        return 0;
      }
      return it2->second.size();
    }
    size_t CountOutVertex() const {
      auto it = this->graph_->out_index_vertex_.find(this->vertex_);
      if (it == this->graph_->out_index_vertex_.end()) {
        return 0;
      }
      return it->second.size();
    }
    size_t CountOutVertex(const EdgeLabelType edge_label) const {
      auto it = this->graph_->out_edge_label_index_vertex_.find(this->vertex_);
      if (it == this->graph_->out_edge_label_index_vertex_.end()) {
        return 0;
      }
      auto it1 = it->second.find(edge_label);
      if (it1 == it->second.end()) {
        return 0;
      }
      return it1->second.size();
    }
    size_t CountOutVertex(const EdgeLabelType edge_label,
                          const VertexLabelType vertex_label) const {
      auto it = this->graph_->out_vertex_label_edge_label_index_vertex_.find(
          this->vertex_);
      if (it == this->graph_->out_vertex_label_edge_label_index_vertex_.end()) {
        return 0;
      }
      auto it1 = it->second.find(vertex_label);
      if (it1 == it->second.end()) {
        return 0;
      }
      auto it2 = it1->second.find(edge_label);
      if (it2 == it1->second.end()) {
        return 0;
      }
      return it2->second.size();
    }

    size_t CountInVertex() const {
      auto it = this->graph_->in_index_vertex_.find(this->vertex_);
      if (it == this->graph_->in_index_vertex_.end()) {
        return 0;
      }
      return it->second.size();
    }
    size_t CountInVertex(const EdgeLabelType edge_label) const {
      auto it = this->graph_->in_edge_label_index_vertex_.find(this->vertex_);
      if (it == this->graph_->in_edge_label_index_vertex_.end()) {
        return 0;
      }
      auto it1 = it->second.find(edge_label);
      if (it1 == it->second.end()) {
        return 0;
      }
      return it1->second.size();
    }
    size_t CountInVertex(const EdgeLabelType edge_label,
                         const VertexLabelType vertex_label) const {
      auto it = this->graph_->in_vertex_label_edge_label_index_vertex_.find(
          this->vertex_);
      if (it == this->graph_->in_vertex_label_edge_label_index_vertex_.end()) {
        return 0;
      }
      auto it1 = it->second.find(vertex_label);
      if (it1 == it->second.end()) {
        return 0;
      }
      auto it2 = it1->second.find(edge_label);
      if (it2 == it1->second.end()) {
        return 0;
      }
      return it2->second.size();
    }
    bool operator==(const _VertexWithIndex &b) const {
      if (!graph_) {
        return !b.graph_;
      } else if (!b.graph_) {
        return false;
      } else {
        return this->vertex_ == b.vertex_;
      }
    }
    bool operator<(const _VertexWithIndex &b) const {
      if (!graph_ || !b.graph_) {
        return false;
      } else {
        return this->vertex_ < b.vertex_;
      }
    }

   public:
    void Reset() { graph_ = nullptr; }

    bool HasValue() const { return (graph_ != nullptr); }
    FragmentVertexType vertex() const { return this->vertex_; }
    typename std::conditional<is_const, const GraphType *, GraphType *>::type
        graph_;
    FragmentVertexType vertex_;
  };

  template <bool is_const>
  class _Edge {
   public:
    using GraphType = FragmentGraphWithIndex;
    using IDType = EdgeIDType;
    using LabelType = EdgeLabelType;
    _Edge() : graph_(nullptr){};
    _Edge(typename std::conditional<is_const, const GraphType *,
                                    GraphType *>::type graph,
          const FragmentEdgeType &edge)
        : graph_(graph), edge_(&edge) {}

   public:
    const EdgeIDType id() const { return edge_->data.id_; }
    const EdgeLabelType &label() const { return edge_->data.label_; }
    const VertexIDType src_id() const {
      return graph_->fragment->GetId(edge_->neighbor);
    }
    const VertexIDType dst_id() const {
      return graph_->fragment->GetId(edge_->neighbor);
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
        return this->edge_ == b.edge_;
      }
    }

    bool operator<(const _Edge &b) const {
      if (!graph_ || !b.graph_) {
        return false;
      } else {
        return this->edge_ < b.edge_;
      }
    }

    void Reset() { graph_ = nullptr; }

    bool HasValue() const { return (graph_ != nullptr) && (edge_ != nullptr); }

   private:
    typename std::conditional<is_const, const GraphType *, GraphType *>::type
        graph_;
    const FragmentEdgeType *edge_;
  };

 public:
  VertexPtr FindVertex(const VertexIDType &id) {
    FragmentVertexType vertex;
    fragment->GetVertex(id, vertex);
    return VertexPtr(Vertex(this, vertex));
  }

  VertexConstPtr FindConstVertex(const VertexIDType &id) const {
    FragmentVertexType vertex;
    fragment->GetVertex(id, vertex);
    return VertexConstPtr(ConstVertex(this, vertex));
  }

  VertexIterator VertexBegin() {
    return VertexIterator(this, this->fragment->Vertices().begin(),
                                this->fragment->Vertices().end());
  }
  VertexConstIterator VertexBegin() const {
    return VertexConstIterator(this, this->fragment->Vertices().begin(),
                                     this->fragment->Vertices().end());
  }

  VertexIndexIterator VertexBegin(const VertexLabelType vertex_label) {
    auto it = this->label_index_vertex_.find(vertex_label);
    if (it == this->label_index_vertex_.end()) {
      return VertexIndexIterator();
    }
    auto &&data = it->second;
    return VertexIndexIterator(this, data.begin(), data.end());
  }
  VertexIndexConstIterator VertexBegin(
      const VertexLabelType vertex_label) const {
    auto it = this->label_index_vertex_.find(vertex_label);
    if (it == this->label_index_vertex_.end()) {
      return VertexIndexConstIterator();
    }
    auto &&data = it->second;
    return VertexIndexConstIterator(this, data.cbegin(), data.cend());
  }

  void BuildIndex() {
    for (auto vertex_it = this->VertexBegin(); 
             !vertex_it.IsDone();
              vertex_it++) {
      VertexIDType vertex_id = vertex_it->id();
      VertexLabelType vertex_label = vertex_it->label();
      auto vertex1 = vertex_it->vertex();
      this->label_index_vertex_[vertex_label].push_back(vertex1);
      FragmentVertexType vertex;
      fragment->GetVertex(vertex_id, vertex);
      for (auto nbt_ptr =
               this->fragment->GetOutgoingAdjList(vertex).begin_pointer();
           nbt_ptr != this->fragment->GetOutgoingAdjList(vertex).end_pointer();
           nbt_ptr++) {
        auto edge_label = (*nbt_ptr).data.label_;
        auto adj_vertex = (*nbt_ptr).neighbor;
        auto adj_vertex_label = fragment->GetData(adj_vertex).label_;
        this->out_edge_label_index_[vertex][edge_label].push_back((nbt_ptr));
        this->out_vertex_edge_label_index_[vertex][adj_vertex][edge_label]
            .push_back(nbt_ptr);
        this->out_edge_label_index_vertex_[vertex][edge_label].push_back(
            adj_vertex);
        this->out_index_vertex_[vertex].push_back(adj_vertex);
        this->out_vertex_label_edge_label_index_[vertex][adj_vertex_label]
                                                [edge_label]
                                                    .push_back(nbt_ptr);
        this
            ->out_vertex_label_edge_label_index_vertex_[vertex]
                                                       [adj_vertex_label]
                                                       [edge_label]
            .push_back(adj_vertex);
      }
      for (auto nbt_ptr =
               this->fragment->GetIncomingAdjList(vertex).begin_pointer();
           nbt_ptr != this->fragment->GetIncomingAdjList(vertex).end_pointer();
           nbt_ptr++) {
        auto edge_label = (*nbt_ptr).data.label_;
        auto adj_vertex = (*nbt_ptr).neighbor;
        auto adj_vertex_label = fragment->GetData(adj_vertex).label_;
        this->in_edge_label_index_[vertex][edge_label].push_back((nbt_ptr));
        this->in_vertex_edge_label_index_[vertex][adj_vertex][edge_label]
            .push_back(nbt_ptr);
        this->in_edge_label_index_vertex_[vertex][edge_label].push_back(
            adj_vertex);
        this->in_index_vertex_[vertex].push_back(adj_vertex);
        this->in_vertex_label_edge_label_index_[vertex][adj_vertex_label]
                                               [edge_label]
                                                   .push_back(nbt_ptr);
        this
            ->in_vertex_label_edge_label_index_vertex_[vertex][adj_vertex_label]
                                                      [edge_label]
            .push_back(adj_vertex);
      }
      for (auto &vertex_list : this->out_edge_label_index_vertex_[vertex]) {
        std::sort(vertex_list.second.begin(), vertex_list.second.end());
        auto unique_it =
            std::unique(vertex_list.second.begin(), vertex_list.second.end());
        vertex_list.second.erase(unique_it, vertex_list.second.end());
      }

      for (auto &vertex_list1 :
           this->out_vertex_label_edge_label_index_vertex_[vertex]) {
        for (auto &vertex_list : vertex_list1.second) {
          std::sort(vertex_list.second.begin(), vertex_list.second.end());
          auto unique_it =
              std::unique(vertex_list.second.begin(), vertex_list.second.end());
          vertex_list.second.erase(unique_it, vertex_list.second.end());
        }
      }

      auto &out_vertex_list = this->out_index_vertex_[vertex];
      std::sort(out_vertex_list.begin(), out_vertex_list.end());
      auto unique_it =
          std::unique(out_vertex_list.begin(), out_vertex_list.end());
      out_vertex_list.erase(unique_it, out_vertex_list.end());

      for (auto &vertex_list : this->in_edge_label_index_vertex_[vertex]) {
        std::sort(vertex_list.second.begin(), vertex_list.second.end());
        auto unique_it =
            std::unique(vertex_list.second.begin(), vertex_list.second.end());
        vertex_list.second.erase(unique_it, vertex_list.second.end());
      }

      for (auto &vertex_list1 :
           this->in_vertex_label_edge_label_index_vertex_[vertex]) {
        for (auto &vertex_list : vertex_list1.second) {
          std::sort(vertex_list.second.begin(), vertex_list.second.end());
          auto unique_it =
              std::unique(vertex_list.second.begin(), vertex_list.second.end());
          vertex_list.second.erase(unique_it, vertex_list.second.end());
        }
      }

      auto &in_vertex_list = this->in_index_vertex_[vertex];
      std::sort(in_vertex_list.begin(), in_vertex_list.end());
      unique_it = std::unique(in_vertex_list.begin(), in_vertex_list.end());
      in_vertex_list.erase(unique_it, in_vertex_list.end());
    }
  }

 public:
  Fragment *fragment;
  std::map<FragmentVertexType, LabelIndexEdgeContainer> out_edge_label_index_;
  std::map<FragmentVertexType, LabelIndexEdgeContainer>  in_edge_label_index_;
  std::map<FragmentVertexType, VertexLabelIndexEdgeContainer>
      out_vertex_edge_label_index_;
  std::map<FragmentVertexType, VertexLabelIndexEdgeContainer>
      in_vertex_edge_label_index_;
  std::map<FragmentVertexType, LabelIndexVertexContainer>
      out_edge_label_index_vertex_;
  std::map<FragmentVertexType, LabelIndexVertexContainer>
      in_edge_label_index_vertex_;
  std::map<FragmentVertexType, VertexLabelEdgeLabelIndexEdgeContainer>
      out_vertex_label_edge_label_index_;
  std::map<FragmentVertexType, VertexLabelEdgeLabelIndexEdgeContainer>
      in_vertex_label_edge_label_index_;
  std::map<FragmentVertexType, VertexLabelEdgeLabelIndexVertexContainer>
      out_vertex_label_edge_label_index_vertex_;
  std::map<FragmentVertexType, VertexLabelEdgeLabelIndexVertexContainer>
      in_vertex_label_edge_label_index_vertex_;

  std::map<FragmentVertexType, IndexVertexContainer> out_index_vertex_;
  std::map<FragmentVertexType, IndexVertexContainer> in_index_vertex_;
  std::map<VertexLabelType, IndexVertexContainer> label_index_vertex_;
};
}  // namespace grape

#endif