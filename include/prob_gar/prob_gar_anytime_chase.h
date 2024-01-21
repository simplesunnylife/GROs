#ifndef _PROB_GAR_PROB_GAR_CHASE_ANYTIME_H
#define _PROB_GAR_PROB_GAR_CHASE_ANYTIME_H

#include <vector>
#include <chrono>

#include "include/gundam/tool/find_edge.h"
#include "include/gundam/tool/map_edge_to.h"
#include "include/gundam/tool/is_match_to.h"

#include "include/gundam/type_getter/edge_id.h"

#include "include/gundam/algorithm/match_using_match.h"

#include "include/logic/disjunctive_paradigm.h"

#include "include/util/is_number.h"
#include "include/util/log.h"

namespace prob_gar {

namespace _prob_gar_anytime_chase {

template <typename GraphPatternType,
          typename    DataGraphType>
class Atom {
 private:
  using EdgeIdType = typename GUNDAM::EdgeID<DataGraphType>::type;

  inline static std::map<std::string, uint32_t> match_hash_str_to_int_;

  inline static std::map<uint32_t, std::string> match_hash_int_to_str_;

  static uint32_t get_match_hash_str_to_int(const std::string& match_hash_str) {
    auto [ match_hash_str_to_int_it,
           match_hash_str_to_int_ret ] = Atom::match_hash_str_to_int_.emplace(match_hash_str, 
                                                                        Atom::match_hash_str_to_int_.size());
    if (match_hash_str_to_int_ret) {
      // does not have this match before
      Atom::match_hash_int_to_str_.emplace(match_hash_str_to_int_it->second, 
                                           match_hash_str);
    }
    return match_hash_str_to_int_it->second;
  }

 public:
  Atom() : is_prob_gar_match_(true) {
    return;
  }

  ~Atom() = default;

  Atom(const Atom&) = default;
  Atom(Atom&&) = default;

  Atom& operator=(const Atom&) = default;	  
  Atom& operator=(Atom&&) = default;

  // match of prob gar in data graph
  Atom(size_t   gar_idx, 
       const std::string& match_hash) 
    : is_prob_gar_match_(true),
               gar_idx_ (  gar_idx),
         match_hash_int_(Atom::get_match_hash_str_to_int(match_hash)) {
    assert( this->is_prob_gar_match());
    assert(!this->is_edge());
    return;
  }

  // edge 
  Atom(EdgeIdType edge_id) 
  : is_prob_gar_match_(false),
              edge_id_(edge_id) {
    assert(!this->is_prob_gar_match());
    assert( this->is_edge());
    return;
  }

  inline bool is_prob_gar_match() const {
    return this->is_prob_gar_match_;
  }

  inline size_t prob_gar_idx() const {
    assert(this->is_prob_gar_match());
    return this->gar_idx_;
  }

  inline const auto& match_hash() const {
    assert(this->is_prob_gar_match());
    return this->match_hash_int_;
  }

  inline bool is_edge() const {
    return !(this->is_prob_gar_match());
  }

  inline EdgeIdType edge_id() const {
    // assert(this->is_edge());
    return this->edge_id_;
  }

  inline bool operator<(const Atom& atom) const {
    // first compare atom type
    if (this->is_prob_gar_match_ < atom.is_prob_gar_match_) {
      return true;
    }
    if (this->is_prob_gar_match_ > atom.is_prob_gar_match_) {
      return false;
    }
    if (this->is_edge()) {
      assert(atom.is_edge());
      // is an edge in data graph, compare only edge id
      // do not care gar_idx_ and match_hash_int_
      return this->edge_id_ < atom.edge_id_;
    }
    assert(this->is_prob_gar_match()
         && atom.is_prob_gar_match());
    // is a match of gar in data graph
    // first compare gar_idx_
    if (this->gar_idx_ < atom.gar_idx_) {
      return true;
    }
    if (this->gar_idx_ > atom.gar_idx_) {
      return false;
    }
    // then compare match_hash_int_, does not care
    // edge id
    return this->match_hash_int_ < atom.match_hash_int_;
  }

  inline bool operator==(const Atom& atom) const {
    // first compare atom type
    if (this->is_prob_gar_match_ != atom.is_prob_gar_match_) {
      return false;
    }
    if (this->is_edge()) {
      assert(atom.is_edge());
      // is an edge in data graph, compare only edge id
      // do not care gar_idx_ and match_hash_int_
      return (this->edge_id_ == atom.edge_id_);
    }
    return (this->gar_idx_ == atom.gar_idx_)
        && (this->match_hash_int_ == atom.match_hash_int_);
  }

  inline std::string ToString() const {
    if (this->is_edge()) {
      return "# edge: " + GUNDAM::ToString(this->edge_id_) + " #";
    }
    return "# prob gar match: <" + GUNDAM::ToString(this->gar_idx_) 
                          + ", " + Atom::match_hash_int_to_str_[this->match_hash_int_] + "> #";
    return "# prob gar match: <" + GUNDAM::ToString(this->gar_idx_) 
                            + ", " + std::to_string(this->match_hash_int_) + "> #";
  }

  static void ResetDict() {
    Atom::match_hash_int_to_str_.clear();
    Atom::match_hash_str_to_int_.clear();
    return;
  }

 private:
  bool is_prob_gar_match_;

  // match of gar
  size_t   gar_idx_;
  uint32_t match_hash_int_;

  // edge
  EdgeIdType edge_id_;
};

template <typename ConfType,
          typename GraphPatternType,
          typename    DataGraphType>
ConfType CalcConfidence(
    const 
     std::vector<ProbGraphAssociationRule<GraphPatternType, 
                                             DataGraphType>>& prob_gar_list, 
                        const logic::DisjunctiveParadigm<
                                     Atom<GraphPatternType,
                                             DataGraphType>>& disjunctive_paradigm) {

  util::Debug("disjunctive_paradigm.ToString(): "
             + disjunctive_paradigm.ToString());

  if (disjunctive_paradigm.Empty()) {
    return 1.0;
  }

  ConfType conf = 1.0;
  for (const auto& intersection_set : disjunctive_paradigm) {
    ConfType conjunctive_conf = 1.0;
    util::Debug("conf before: " + std::to_string(conf));
    util::Debug("intersection_set.ToString(): "
                + intersection_set.ToString());
    for (const auto& atom : intersection_set) {
      util::Debug("\tatom: " + atom.ToString());
      assert(atom.is_prob_gar_match());
      assert(atom.prob_gar_idx() >= 0 
          && atom.prob_gar_idx() <  prob_gar_list.size());
      conjunctive_conf *= prob_gar_list[atom.prob_gar_idx()].prob();
      util::Debug("\tconjunctive_conf: " + std::to_string(conjunctive_conf));
    }
    conf *= 1 - conjunctive_conf;
    util::Debug("conf after: " + std::to_string(conf));
  }
  return (1.0 - conf);
}

}; // _prob_gar_anytime_chase

// forward inference
template <typename GraphPatternType, 
          typename    DataGraphType>
int ProbGARAnytimeChase(
          std::vector<ProbGraphAssociationRule<GraphPatternType, 
                                                  DataGraphType>> &prob_gar_list, 
    DataGraphType &data_graph,
    std::map<typename GUNDAM::EdgeID<DataGraphType>::type, float>& diff_edge_set,
    float conf_bound =  1.0,
    float time_limit = -1.0) {

  using EdgeIDType = typename GUNDAM::EdgeID<DataGraphType>::type;
  using EdgeLabelType = typename GUNDAM::EdgeLabel<DataGraphType>::type;

  using VertexIDType = typename GUNDAM::VertexID<DataGraphType>::type;

  using VertexHandleType = typename GUNDAM::VertexHandle<DataGraphType>::type;
  
  using PatternEdgeHandleType = typename GUNDAM::EdgeHandle<GraphPatternType>::type;

  EdgeIDType max_data_graph_edge_id = 0;
  for (auto vertex_it = data_graph.VertexBegin();
           !vertex_it.IsDone();
            vertex_it++) {
    for (auto out_edge_it = vertex_it->OutEdgeBegin();
             !out_edge_it.IsDone();
              out_edge_it++) {
      max_data_graph_edge_id = max_data_graph_edge_id > out_edge_it->id() ?
                               max_data_graph_edge_id : out_edge_it->id();
    } 
  }
  EdgeIDType edge_id_counter = max_data_graph_edge_id;
  edge_id_counter++;

  using namespace _prob_gar_anytime_chase;

  using namespace logic;

  // atom store in the disjunctive paradigm
  // has two types: 
  //    1, edge
  //    2, match of a gar in data graph
  using AtomType = Atom<GraphPatternType,
                           DataGraphType>;

  using IntersectionSetType 
      = IntersectionSet<AtomType>;

  using DisjunctiveParadigmType 
      = DisjunctiveParadigm<AtomType>;

  using MatchType = GUNDAM::Match<GraphPatternType,
                                     DataGraphType>;

  // use hash code to distinguish the different between the different
  // 
  using MatchHashType = typename MatchType::HashType;

  //  tcp for each atom
  //  the atom can only contain two kinds of item:
  //  a match of gar in data graph, or a chased edge
  std::map<AtomType, DisjunctiveParadigmType> tcp_set;
  
  // set of atoms updated or added in last round
  std::vector<AtomType> atom_set_modified_in_last_round;

  // to mark whether is the first round,
  // some optimization can be added in first round
  bool is_first_round = true;

  auto begin_time = std::chrono::high_resolution_clock::now();

  size_t round = 0;

  using ConfType = float;

  while (time_limit < 0 // no time limit
     || (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - begin_time).count()) < time_limit // still have time
        ) {

    util::Info("==========================================");
    util::Info("    round: " + std::to_string(round) );
    round++;

    auto match_begin_time = std::chrono::high_resolution_clock::now();

    // if is not the first round, then only needs to collect
    // the vertexs set that connects to the updated atoms(edges)
    std::vector<VertexHandleType> delta_vertex_set;
    if (!atom_set_modified_in_last_round.empty()) {
      util::Debug("wenzhi: atom_set_modified_in_last_round.size(): "
          + std::to_string(atom_set_modified_in_last_round.size()));
      // is not first round 
      assert(!is_first_round);
      for (const auto& atom : atom_set_modified_in_last_round) {
        assert(atom.is_edge());
        if (!atom.is_edge()) {
          assert(false);
          continue;
        }
        auto edge_handle = data_graph.FindEdge(atom.edge_id());
        assert(edge_handle);
        delta_vertex_set.emplace_back(data_graph.FindVertex(edge_handle->src_handle()->id()));
        assert(delta_vertex_set.back());
        delta_vertex_set.emplace_back(data_graph.FindVertex(edge_handle->dst_handle()->id()));
        assert(delta_vertex_set.back());
      }
      std::sort(delta_vertex_set.begin(),
                delta_vertex_set.end());
      delta_vertex_set.erase(std::unique(delta_vertex_set.begin(), 
                                         delta_vertex_set.end()),
                                         delta_vertex_set.end());
      // since there are atom modified, there should be vertex updated
      assert(!delta_vertex_set.empty());
      // util::Info("delta_vertex_set: ");
      // for (const auto& delta_vertex_handle : delta_vertex_set) {
      //   util::Info("\t" + std::to_string(delta_vertex_handle->id()));
      // }
      // util::Info("atom_set_modified_in_last_round: ");
      // for (const auto& atom : atom_set_modified_in_last_round) {
      //   assert(atom.is_edge());
      //   if (!atom.is_edge()) {
      //     assert(false);
      //     continue;
      //   }
      //   util::Info("\t" + std::to_string(atom.edge_id()));
      // }
    }

    // since the data graph cannot be modified in the callback
    // during matching process, then store all udpated edges and
    // add them into data graph at once
    std::map<std::tuple<VertexIDType,  //  src_id
                        VertexIDType,  //  dst_id
                       EdgeLabelType>, // edge_label
              size_t> // edge_id 
              new_edges_to_add;

    // collects all atoms that are updated in this round
    std::map<AtomType, DisjunctiveParadigmType> current_round_tcp_set;

    // first collects all atom that needs to update
    // e.g. to find all matchs from gar to data graph that contains 
    //      all modified edges (all matches if is first round)
    for (size_t prob_gar_idx = 0;
                prob_gar_idx < prob_gar_list.size();
                prob_gar_idx++) {
      auto& prob_gar = prob_gar_list[prob_gar_idx];
      const auto rhs_literal_info = (*prob_gar.y_literal_set().begin())->info();
      // only support gar set with one rhs literal and it needs to be edge literal
      assert(prob_gar.y_literal_set().Count() == 1);
      assert(rhs_literal_info.literal_type() == gar::LiteralType::kEdgeLiteral);

      size_t match_counter = 0;

      // // #ifndef NDEBUG 
      // if (data_graph.FindVertex(75502)
      //  && data_graph.FindVertex(61290)
      //  && data_graph.FindVertex(55974)
      //  && data_graph.FindVertex(8826)
      //  && data_graph.FindVertex(81527)) {

      //   MatchType temp_match;
      //   temp_match.AddMap(prob_gar.pattern().FindVertex(0),
      //                             data_graph.FindVertex(75502));
      //   temp_match.AddMap(prob_gar.pattern().FindVertex(2),
      //                             data_graph.FindVertex(61290));
      //   temp_match.AddMap(prob_gar.pattern().FindVertex(3),
      //                             data_graph.FindVertex(55974));
      //   temp_match.AddMap(prob_gar.pattern().FindVertex(4),
      //                             data_graph.FindVertex(8826));
      //   temp_match.AddMap(prob_gar.pattern().FindVertex(5),
      //                             data_graph.FindVertex(81527));

      //   if (GUNDAM::IsMatchTo(prob_gar.pattern(),
      //                                 data_graph,
      //                                 temp_match)) {
      //     util::Info("found map:");
      //     util::Info("\t<0,1006,2,2> to <75502,1006,61290,302224>");
      //     util::Info("\t<0,1013,3,3> to <75502,1013,55974,131462>");
      //     util::Info("\t<4,1013,3,4> to <8826,1013,55974,176447>");
      //     util::Info("\t<4,1001,5,5> to <8826,1001,81527,30078>");
      //     util::Info("\t<3,1001,5,6> to <55974,1001,81527,116250>");

      //     if (!std::binary_search(delta_vertex_set.begin(),
      //                             delta_vertex_set.end(),
      //                             data_graph.FindVertex(75502))
      //      && !std::binary_search(delta_vertex_set.begin(),
      //                             delta_vertex_set.end(),
      //                             data_graph.FindVertex(61290))
      //      && !std::binary_search(delta_vertex_set.begin(),
      //                             delta_vertex_set.end(),
      //                             data_graph.FindVertex(55974))
      //      && !std::binary_search(delta_vertex_set.begin(),
      //                             delta_vertex_set.end(),
      //                             data_graph.FindVertex(8826))
      //      && !std::binary_search(delta_vertex_set.begin(),
      //                             delta_vertex_set.end(),
      //                             data_graph.FindVertex(81527))) {
      //       util::Info("\tnot contained in delta_vertex_set");
      //     }
      //   }
      //   else {
      //     util::Info("does not find map:");
      //     util::Info("\t<0,1006,2,2> to <75502,1006,61290,302224>");
      //     util::Info("\t<0,1013,3,3> to <75502,1013,55974,131462>");
      //     util::Info("\t<4,1013,3,4> to <8826,1013,55974,176447>");
      //     util::Info("\t<4,1001,5,5> to <8826,1001,81527,30078>");
      //     util::Info("\t<3,1001,5,6> to <55974,1001,81527,116250>");
      //   }
      // }
      // else {
      //   util::Info("does not find map:");
      //   util::Info("\t<0,1006,2,2> to <75502,1006,61290,302224>");
      //   util::Info("\t<0,1013,3,3> to <75502,1013,55974,131462>");
      //   util::Info("\t<4,1013,3,4> to <8826,1013,55974,176447>");
      //   util::Info("\t<4,1001,5,5> to <8826,1001,81527,30078>");
      //   util::Info("\t<3,1001,5,6> to <55974,1001,81527,116250>");
      // }
      // // #endif // NDEBUG

      // for each match from gar to data graph, consider:
      //    1, whether it can generate new edges except the ones in original data graph
      //    2, whether it maps the pattern in gar to edges in data graph modified in last round
      // further process only if it satisfy both conditions
      std::function<bool(const MatchType&)> 
           match_callback = [&is_first_round,
                             &prob_gar,
                             &prob_gar_idx,
                             &data_graph,
                             &current_round_tcp_set,
                             &new_edges_to_add,
                             &max_data_graph_edge_id,
                             &edge_id_counter,
                             &  match_counter,
                             &atom_set_modified_in_last_round] (const MatchType& match) -> bool {

        // first find whether such a match can chase new edges that 
        // are not contained in the original data graph
        
        assert(prob_gar.y_literal_set().Count() == 1);
        // to find the chased edge 
        const auto  y_literal_info
        =(*prob_gar.y_literal_set().begin())->info();
        // support only edge literal
        assert(y_literal_info.literal_type() == gar::LiteralType::kEdgeLiteral);

        AtomType atom_to_update;

        // if (match.MapTo(prob_gar.pattern().FindVertex(0))->id() == 75502 
        //  && match.MapTo(prob_gar.pattern().FindVertex(2))->id() == 61290
        //  && match.MapTo(prob_gar.pattern().FindVertex(3))->id() == 55974 
        //  && match.MapTo(prob_gar.pattern().FindVertex(4))->id() == 8826
        //  && match.MapTo(prob_gar.pattern().FindVertex(5))->id() == 81527) {
        //   util::Info("found map in callback:");
        //   util::Info("\t<0,1006,2,2> to <75502,1006,61290,302224>");
        //   util::Info("\t<0,1013,3,3> to <75502,1013,55974,131462>");
        //   util::Info("\t<4,1013,3,4> to <8826,1013,55974,176447>");
        //   util::Info("\t<4,1001,5,5> to <8826,1001,81527,30078>");
        //   util::Info("\t<3,1001,5,6> to <55974,1001,81527,116250>");
        // }

        util::Debug("wenzhi here0");
               
        assert(match.MapTo(y_literal_info.x_id()));
        assert(match.MapTo(y_literal_info.y_id()));

        assert(match.MapTo(y_literal_info.x_id())
            == data_graph.FindVertex(match.MapTo(y_literal_info.x_id())->id()));
        assert(match.MapTo(y_literal_info.y_id())
            == data_graph.FindVertex(match.MapTo(y_literal_info.y_id())->id()));
            
        // to find whether it can generate new edge in data graph
        auto edge_handle = GUNDAM::FindEdge<DataGraphType>(
                                    match.MapTo(y_literal_info.x_id()),
                                                y_literal_info.edge_label(),
                                    match.MapTo(y_literal_info.y_id()));

        if (edge_handle) {
          assert( data_graph.FindEdge(edge_handle->id()) );
          assert(edge_handle->label()      ==             y_literal_info.edge_label());
          assert(edge_handle->src_handle() == match.MapTo(y_literal_info.x_id()));
          assert(edge_handle->dst_handle() == match.MapTo(y_literal_info.y_id()));
          // has satisfied edges
          if (is_first_round) {
            // if is first round, then it means that the chased
            // edge is contained in the original graph
            return true;
          }
          assert(!is_first_round);
          // if is not first round, then it needs to find whether
          // such a match containes new edges modified in last round
          bool contains_edge_modified_in_last_round = false;
          for (auto edge_it = prob_gar.pattern().EdgeBegin();
                   !edge_it.IsDone();
                    edge_it++) {

            PatternEdgeHandleType pattern_edge_handle = edge_it;

             auto  edge_handle_in_data_graph = GUNDAM::MapEdgeTo(match, pattern_edge_handle);
            assert(edge_handle_in_data_graph);
            assert(edge_handle_in_data_graph->label()      ==             pattern_edge_handle->label());
            assert(edge_handle_in_data_graph->src_handle() == match.MapTo(pattern_edge_handle->src_handle()));
            assert(edge_handle_in_data_graph->dst_handle() == match.MapTo(pattern_edge_handle->dst_handle()));

            AtomType edge_atom(edge_handle_in_data_graph->id());
            
            assert(std::is_sorted(atom_set_modified_in_last_round.begin(),
                                  atom_set_modified_in_last_round.end()));
            if (!std::binary_search(atom_set_modified_in_last_round.begin(),
                                    atom_set_modified_in_last_round.end(),
                                        edge_atom)) {
              // is not an edge updated or added in last round
              continue;
            }
            contains_edge_modified_in_last_round = true;
            break;
          }
          if (!contains_edge_modified_in_last_round) {
            return true;
          }
          // to find whether it is contained in the original data graph
          util::Debug("#### wenzhi: has satisfy edge ####");

          // has found this edge in data graph
          if (edge_handle->id() <= max_data_graph_edge_id) {
            // this edge is contained in original data graph
            // this match does not need to be considered
            util::Debug("#### wenzhi: is original edge ####");
            return true;
          }
          util::Debug("#### wenzhi: is not original edge ####");
          atom_to_update = AtomType(edge_handle->id());
          assert( data_graph.FindEdge(atom_to_update.edge_id()) );
          assert((data_graph.FindEdge(atom_to_update.edge_id())->src_handle() == match.MapTo(y_literal_info.x_id()))
              && (data_graph.FindEdge(atom_to_update.edge_id())->dst_handle() == match.MapTo(y_literal_info.y_id())));
        }
        else {
          // #ifndef NDEBUG
          if (!is_first_round) {
            bool contains_edge_modified_in_last_round = false;
            for (auto edge_it = prob_gar.pattern().EdgeBegin();
                     !edge_it.IsDone();
                      edge_it++) {

              PatternEdgeHandleType pattern_edge_handle = edge_it;

               auto  edge_handle_in_data_graph = GUNDAM::MapEdgeTo(match, pattern_edge_handle);
              assert(edge_handle_in_data_graph);
              assert(edge_handle_in_data_graph->label()      ==             pattern_edge_handle->label());
              assert(edge_handle_in_data_graph->src_handle() == match.MapTo(pattern_edge_handle->src_handle()));
              assert(edge_handle_in_data_graph->dst_handle() == match.MapTo(pattern_edge_handle->dst_handle()));

              AtomType edge_atom(edge_handle_in_data_graph->id());
              
              assert(std::is_sorted(atom_set_modified_in_last_round.begin(),
                                    atom_set_modified_in_last_round.end()));
              if (!std::binary_search(atom_set_modified_in_last_round.begin(),
                                      atom_set_modified_in_last_round.end(),
                                          edge_atom)) {
                // is not an edge updated or added in last round
                continue;
              }
              contains_edge_modified_in_last_round = true;
              break;
            }
            if (!contains_edge_modified_in_last_round) {
              util::Info("not contains_edge_modified_in_last_round!"
                        + std::to_string(contains_edge_modified_in_last_round)
                        + "<" + std::to_string(match.MapTo(y_literal_info.x_id())->id())
                        + "," + std::to_string(            y_literal_info.edge_label())
                        + "," + std::to_string(match.MapTo(y_literal_info.y_id())->id())
                        + ">");
              for (auto edge_it = prob_gar.pattern().EdgeBegin();
                       !edge_it.IsDone();
                        edge_it++) {
                PatternEdgeHandleType pattern_edge_handle = edge_it;

                 auto  edge_handle_in_data_graph = GUNDAM::MapEdgeTo(match, pattern_edge_handle);
                assert(edge_handle_in_data_graph);
                assert(edge_handle_in_data_graph->label()      ==             pattern_edge_handle->label());
                assert(edge_handle_in_data_graph->src_handle() == match.MapTo(pattern_edge_handle->src_handle()));
                assert(edge_handle_in_data_graph->dst_handle() == match.MapTo(pattern_edge_handle->dst_handle()));

                util::Info("\t\t<" + std::to_string(edge_it->src_handle()->id())
                             + "," + std::to_string(edge_it->label())
                             + "," + std::to_string(edge_it->dst_handle()->id())
                             + "," + std::to_string(edge_it->id())
                        + "> to <" + std::to_string(edge_handle_in_data_graph->src_handle()->id())
                             + "," + std::to_string(edge_handle_in_data_graph->label())
                             + "," + std::to_string(edge_handle_in_data_graph->dst_handle()->id())
                             + "," + std::to_string(edge_handle_in_data_graph->id())
                             + ">");
              }
              for (auto out_edge_it = match.MapTo(y_literal_info.x_id())->OutEdgeBegin();
                       !out_edge_it.IsDone();
                        out_edge_it++) {
                util::Info("\t<" + std::to_string(out_edge_it->src_handle()->id())
                           + "," + std::to_string(out_edge_it->label())
                           + "," + std::to_string(out_edge_it->dst_handle()->id())
                           + ">");
              }
            }
            assert(contains_edge_modified_in_last_round);
          }
          // #endif // NDEBUG
          util::Debug("######################################");
          util::Debug("#### wenzhi: has not satisfy edge ####");

          util::Debug("#### wenzhi: match.MapTo(y_literal_info.x_id())->id(): "
                   + std::to_string(match.MapTo(y_literal_info.x_id())->id()));

          util::Debug("#### wenzhi: match.MapTo(y_literal_info.y_id())->id(): "
                   + std::to_string(match.MapTo(y_literal_info.y_id())->id()));

          util::Debug("#### wenzhi: y_literal_info.edge_label(): "
                   + std::to_string(y_literal_info.edge_label()));

          util::Debug("#### new_edges_to_add.size(): "
           + std::to_string(new_edges_to_add.size()));

          assert(match.MapTo(y_literal_info.x_id()));
          assert(match.MapTo(y_literal_info.y_id()));

          assert(match.MapTo(y_literal_info.x_id())
             == data_graph.FindVertex(match.MapTo(y_literal_info.x_id())->id()));
          assert(match.MapTo(y_literal_info.y_id())
             == data_graph.FindVertex(match.MapTo(y_literal_info.y_id())->id()));

          // this edge is not contained in the data graph
          // add it to new_edges_to_add
          auto edge_it =  new_edges_to_add.find(std::tuple(match.MapTo(y_literal_info.x_id())->id(),
                                                           match.MapTo(y_literal_info.y_id())->id(),
                                                                       y_literal_info.edge_label()));
          if ( edge_it == new_edges_to_add.end() ) {
            util::Debug("#### has not been added to new_edges_to_add ####");
            // has not been added in edge to add
            edge_it = new_edges_to_add.emplace_hint(edge_it,
                          std::tuple(match.MapTo(y_literal_info.x_id())->id(),
                                     match.MapTo(y_literal_info.y_id())->id(),
                                                 y_literal_info.edge_label()),
                                                    edge_id_counter++);
            assert(edge_it->second == (edge_id_counter - 1));
          }
          util::Debug("******************************************");
          util::Debug("*** wenzhi: atom_to_update.edge_id() before: "
                  + std::to_string(atom_to_update.edge_id()));
          util::Debug("******************************************");
          assert(edge_it != new_edges_to_add.end());
          assert(edge_it->first == std::tuple(match.MapTo(y_literal_info.x_id())->id(),
                                              match.MapTo(y_literal_info.y_id())->id(),
                                                          y_literal_info.edge_label()));
          assert(edge_it->second > max_data_graph_edge_id);
          util::Debug("#### edge_it->second: " + std::to_string(edge_it->second) + " ####");
          atom_to_update = AtomType(edge_it->second);
          assert(atom_to_update.is_edge());
          assert(atom_to_update.edge_id() == edge_it->second);
          util::Debug("******************************************");
          util::Debug("*** wenzhi: atom_to_update.edge_id() after: "
                  + std::to_string(atom_to_update.edge_id()));
          util::Debug("******************************************");
          assert(!data_graph.FindEdge(atom_to_update.edge_id()));
        }

        util::Debug("##########################################");
        util::Debug("### wenzhi: atom_to_update.edge_id(): "
                + std::to_string(atom_to_update.edge_id()));
        util::Debug("##########################################");

        assert(atom_to_update.is_edge());

        std::vector<AtomType> modified_edge_set_contained_in_match;
        
        // to find whether there contains edges modified,
        // updated or added, in the last round if is not first
        // round matching
        if (!is_first_round) {
          // would be unnecessary for the first round
          // since all edges are considered at the first time
          bool contain_new_edges = false;
          // iterate over all edges in pattern, to find the edge it
          // maps to in the data graph
          for (auto edge_it = prob_gar.pattern().EdgeBegin();
                   !edge_it.IsDone();
                    edge_it++) {

            PatternEdgeHandleType edge_handle = edge_it;

             auto  edge_handle_in_data_graph = GUNDAM::MapEdgeTo(match, edge_handle);
            assert(edge_handle_in_data_graph);

            AtomType edge_atom(edge_handle_in_data_graph->id());
            
            assert(std::is_sorted(atom_set_modified_in_last_round.begin(),
                                  atom_set_modified_in_last_round.end()));
            if (!std::binary_search(atom_set_modified_in_last_round.begin(),
                                    atom_set_modified_in_last_round.end(),
                                        edge_atom)) {
              // is not an edge updated or added in last round
              continue;
            }
            modified_edge_set_contained_in_match.emplace_back(edge_handle_in_data_graph->id());
          }
          #ifndef NDEBUG 
          if (modified_edge_set_contained_in_match.empty()) {
            // does not contain modified edge, does
            // not need to be considered
            util::Debug("##########################################");
            util::Debug("### wenzhi: atom_to_update.edge_id(): "
                    + std::to_string(atom_to_update.edge_id()));
            util::Debug("### does not contain modified edge");
            util::Debug("##########################################");
            assert(data_graph.FindEdge(atom_to_update.edge_id()));
            for (const auto& [edge_type, edge_id] : new_edges_to_add) {
              assert(atom_to_update.edge_id() != edge_id);
            }
            assert(data_graph.FindEdge(atom_to_update.edge_id()));
            assert(false);
          }
          #endif // NDEBUG
          assert(modified_edge_set_contained_in_match.size() > 0
             && (modified_edge_set_contained_in_match.size() <= prob_gar.pattern().CountEdge()));
        }

        assert(is_first_round || !modified_edge_set_contained_in_match.empty());
        
        if (is_first_round) {
          // as an optimization, the first round match does not need
          // to store the hash, store the index would be sufficient
          // since such a match would not be considered again
          modified_edge_set_contained_in_match.emplace_back(prob_gar_idx, std::to_string(match_counter));
          match_counter++;

          // modified_edge_set_contained_in_match.emplace_back(prob_gar_idx, match.hash());

          IntersectionSetType intersection_set(modified_edge_set_contained_in_match);

          current_round_tcp_set[atom_to_update].Union(intersection_set);
          return true;
        }

        MatchHashType prob_gar_match_hash = match.hash();
        //  = match.hash(prob_gar.pivot_set());

        assert(!util::is_number(prob_gar_match_hash));

        modified_edge_set_contained_in_match.emplace_back(prob_gar_idx,
                                                          prob_gar_match_hash);

        IntersectionSetType intersection_set(modified_edge_set_contained_in_match);

        current_round_tcp_set[atom_to_update].Union(intersection_set);

        return true;
      };

      util::Debug("# wenzhi prob_gar_idx: " 
           + std::to_string(prob_gar_idx));

      if (is_first_round) {
        std::function<bool(const MatchType&)> 
          prune_nothing_callback = [] (const MatchType& match) -> bool {
          // prune nothing
          return false;
        };

        auto ret = GUNDAM::MatchUsingMatch(prob_gar.pattern(),
                                           data_graph,
                                           prune_nothing_callback,
                                                   match_callback);
                            
        util::Debug("# wenzhi ret: " 
             + std::to_string(ret));
        continue;
      }
      // only need to consider the vertexes set that connects
      // to the updated edges 
      GUNDAM::IncrementalMatchUsingMatch(
                      prob_gar.pattern(),
                      data_graph,
                     delta_vertex_set,
                     match_callback);
    }
    atom_set_modified_in_last_round.clear();

    util::Info("    match time(ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - match_begin_time).count()) );

    auto logic_begin_time = std::chrono::high_resolution_clock::now();

    util::Debug("wenzhi: new_edges_to_add.size(): "
        + std::to_string(new_edges_to_add.size()));

    util::Debug("wenzhi: current_round_tcp_set.size(): "
        + std::to_string(current_round_tcp_set.size()));

    #ifndef NDEBUG
    if (new_edges_to_add.size() > current_round_tcp_set.size()) {
      util::Debug("###########################");
      util::Debug("##   new_edges_to_add    ##");
      util::Debug("###########################");
      for (const auto& [edge_type, edge_id] : new_edges_to_add) {
        util::Debug(std::to_string(edge_id));
      }
      util::Debug("###############################");
      util::Debug("##   current_round_tcp_set   ##");
      util::Debug("###############################");
      for (const auto& [atom, tcp_set] : current_round_tcp_set) {
        assert(atom.is_edge());
        util::Debug(std::to_string(atom.edge_id()));
      }
      std::vector<EdgeIDType> edge_id_not_in_current_round_tcp_set;
      for (const auto& [edge_type, edge_id] : new_edges_to_add) {
        AtomType edge_atom(edge_id);
        if (current_round_tcp_set.find(edge_atom)
         == current_round_tcp_set.end()) {
          edge_id_not_in_current_round_tcp_set.emplace_back(edge_id);
        }
      }
      util::Debug("##############################################");
      util::Debug("##   edge_id_not_in_current_round_tcp_set   ##");
      util::Debug("##############################################");
      for (const auto& edge_id : edge_id_not_in_current_round_tcp_set) {
        util::Debug(std::to_string(edge_id));
      }
    }
    #endif // NDEBUG

    assert(new_edges_to_add.size() <= current_round_tcp_set.size());
    #ifndef NDEBUG
    size_t edge_in_current_round_tcp_set = 0;
    for (const auto& [atom, tcp] : current_round_tcp_set) {
      if (atom.is_edge()) {
        edge_in_current_round_tcp_set++;
      }
    }
    assert(new_edges_to_add.size() <= edge_in_current_round_tcp_set);
    #endif // NDEBUG

    // add all new edges into the data graph
    for (const auto& [edge_info, edge_id] : new_edges_to_add) {
      auto [ edge_handle,
             edge_ret ] = data_graph.AddEdge(
                            std::get<0>(edge_info), 
                            std::get<1>(edge_info), 
                            std::get<2>(edge_info), edge_id);
      assert(edge_handle);
      assert(edge_handle->src_handle()->id() == std::get<0>(edge_info));
      assert(edge_handle->dst_handle()->id() == std::get<1>(edge_info));
      assert(edge_handle->label() == std::get<2>(edge_info));
      assert(edge_ret);

      // util::Info("\t<" + std::to_string(std::get<0>(edge_info))
      //            + "," + std::to_string(std::get<2>(edge_info))
      //            + "," + std::to_string(std::get<1>(edge_info))
      //            + ">");
    }
    new_edges_to_add.clear();

    std::map<AtomType, DisjunctiveParadigmType> last_round_tcp_set(tcp_set);

    for (const auto& [current_atom, tcp] : current_round_tcp_set) {
      assert(current_atom.is_edge());
      auto [ tcp_set_it, 
             tcp_set_ret ] = tcp_set.emplace(current_atom, DisjunctiveParadigmType());
      // if tcp_set_ret == true, means it is added successfully,
      // this is a new atom
      if (!tcp_set_ret && tcp_set_it->second.Empty()) {
        // this confidence has already over confidence bound 
        // and have been set to 1
        continue;
      }
      DisjunctiveParadigmType disjunctive_paradigm;
      // util::Debug("###########################");
      // util::Debug("# tcp.Size(): " + std::to_string(tcp.Size()));
      // util::Debug("# tcp.ToString(): " + tcp.ToString());
      // util::Debug("# disjunctive_paradigm.Size(): " + std::to_string(disjunctive_paradigm.Size()));
      // util::Debug("# disjunctive_paradigm.ToString(): " + disjunctive_paradigm.ToString());
      // util::Debug("###########################");
      for (const auto& intersection : tcp) {
        DisjunctiveParadigmType temp_disjunctive_paradigm;
        for (const auto& atom : intersection) {
          // should have been contained in tcp_set
          assert((atom.is_prob_gar_match() && last_round_tcp_set.find(atom) == last_round_tcp_set.end())
              || (atom.is_edge()           && last_round_tcp_set.find(atom) != last_round_tcp_set.end()));
          if (atom.is_prob_gar_match()) {
            assert(tcp_set.find(atom) == tcp_set.end());
            // util::Debug("\t@@@@@@@@@@@@@@@");
            // util::Debug("\t@ atom: " + atom.ToString());
            // util::Debug("\t@ temp_disjunctive_paradigm before Intersect: " + temp_disjunctive_paradigm.ToString());
            temp_disjunctive_paradigm.Intersect(atom);
            // util::Debug("\t@ temp_disjunctive_paradigm after Intersect: "  + temp_disjunctive_paradigm.ToString());
            // util::Debug("\t@@@@@@@@@@@@@@@");
            continue;
          }
          assert(last_round_tcp_set.find(atom) 
              != last_round_tcp_set.end());
          // util::Debug("\t@@@@@@@@@@@@@@@");
          // util::Debug("\t@ intersect tcp of atom: " + atom.ToString());
          // util::Debug("\t@ temp_disjunctive_paradigm before Intersect: " + temp_disjunctive_paradigm.ToString());
          // util::Debug("\t@ tcp_set.find(atom)->second to intersect: " + tcp_set.find(atom)->second.ToString());
          #ifndef NDEBUG
          for (const auto& _intersection : last_round_tcp_set.find(atom)->second) {
            for (const auto& _atom : _intersection) {
              // util::Debug("atom.ToString(): " + _atom.ToString());
              assert(_atom.is_prob_gar_match());
            }
          }
          #endif // NDEBUG
          temp_disjunctive_paradigm.Intersect(last_round_tcp_set.find(atom)->second);
          // util::Debug("\t@ temp_disjunctive_paradigm after Intersect: "  + temp_disjunctive_paradigm.ToString());
          // util::Debug("\t@@@@@@@@@@@@@@@");
        }
        // util::Debug("@@@@@@@@@@@@@@@");
        // util::Debug("@ before Union: " + disjunctive_paradigm.ToString());
        // util::Debug("@ temp_disjunctive_paradigm to Union: " + temp_disjunctive_paradigm.ToString());
        disjunctive_paradigm.Union(temp_disjunctive_paradigm);
        // util::Debug("@ after Union: " + disjunctive_paradigm.ToString());
        // util::Debug("@@@@@@@@@@@@@@@");
      }
      bool modified = tcp_set_it->second.Union(disjunctive_paradigm);
      if (modified || tcp_set_ret) {
        // util::Debug("*********************************************");
        // util::Debug("*  wenzhi: modified!  atom: " + current_atom.ToString());
        // util::Debug("*********************************************");
        atom_set_modified_in_last_round.emplace_back(current_atom);
      }
      #ifndef NDEBUG
      else {
        // util::Debug("*********************************************");
        // util::Debug("*  wenzhi: not modified!  atom: " + current_atom.ToString());
        // util::Debug("*********************************************");
        assert(_prob_gar_anytime_chase::CalcConfidence<ConfType>(prob_gar_list, 
                                                                  tcp_set_it->second) < conf_bound);
      }
      #endif // NDEBUG
      if (_prob_gar_anytime_chase::CalcConfidence<ConfType>(prob_gar_list, 
                                                             tcp_set_it->second) >= conf_bound) {
        // util::Debug("CalcConfidence: " + std::to_string(
        //   _prob_gar_anytime_chase::CalcConfidence<ConfType>(prob_gar_list, 
        //                                                      tcp_set_it->second)));
        tcp_set_it->second.Clear();
        assert(tcp_set_it->second.Empty());
      }
    }

    assert(std::is_sorted(atom_set_modified_in_last_round.begin(),
                          atom_set_modified_in_last_round.end()));

    util::Info("    logic time(ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - logic_begin_time).count()) );
    util::Info("==========================================");

    if (atom_set_modified_in_last_round.empty()) {
      break;
    }

    is_first_round = false;
  }

  util::Debug("====================");
  for (const auto& [atom_to_update, disjunctive_paradigm] : tcp_set) {
    util::Debug("\t=  atom_to_update: " + atom_to_update.ToString());
    util::Debug("\t=  disjunctive_paradigm: " + disjunctive_paradigm.ToString());
  }
  util::Debug("====================");

  // calculate the confidence based on the 
  assert(diff_edge_set.empty());
  diff_edge_set.clear();
  for (const auto& [atom_to_update, disjunctive_paradigm] : tcp_set) {
    assert(atom_to_update.is_edge());
    util::Debug("====================");
    util::Debug("atom_to_update.ToString(): " + atom_to_update.ToString());
    util::Debug("====================");

    ConfType conf = _prob_gar_anytime_chase::CalcConfidence<ConfType>(prob_gar_list, disjunctive_paradigm);

    diff_edge_set.emplace(atom_to_update.edge_id(), conf);
  }

  AtomType::ResetDict();
  
  return diff_edge_set.size();
}

// forward inference
template <typename GraphPatternType, 
          typename    DataGraphType>
int ProbGARAnytimeChase(
    const ProbGraphAssociationRule<GraphPatternType, 
                                      DataGraphType> &prob_gar, 
    DataGraphType &data_graph,
    std::map<typename GUNDAM::EdgeID<DataGraphType>::type, float>& diff_edge_set,
    float conf_bound =  1.0,
    float time_limit = -1.0) {
  std::vector<ProbGraphAssociationRule<GraphPatternType, 
                                          DataGraphType>> prob_gar_list;
  prob_gar_list.emplace_back(prob_gar);
  return ProbGARAnytimeChase(prob_gar_list, data_graph, diff_edge_set, conf_bound, time_limit);
}

}; // namespace prob_gar

#endif  // _PROB_GAR_PROB_GAR_CHASE_ANYTIME_H