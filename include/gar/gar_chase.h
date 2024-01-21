#ifndef GAR_CHASE_H_
#define GAR_CHASE_H_

#include <map>
#include <vector>

#include "gar.h"
#include "gar_match.h"
//#include "literalset.h"

namespace gar {

template <class GraphPatternType, 
          class    DataGraphType, 
          class EdgeIDGen>
int GARChase(
    std::vector<GraphAssociationRule<GraphPatternType, 
                                        DataGraphType>> &gar_list, 
    DataGraphType &data_graph,
    EdgeIDGen &edge_id_gen,
    std::set<typename DataGraphType::VertexType::IDType> *diff_vertex_set = nullptr,
    std::set<typename DataGraphType::  EdgeType::IDType> *diff_edge_set   = nullptr) {
  using GraphPatternVertexHandleType = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandleType = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using MatchMap = std::map<GraphPatternVertexHandleType, 
                               DataGraphVertexHandleType>;

  using MatchResultList = std::vector<MatchMap>;

  int count = 0;
  bool update_flag;
  do {
    update_flag = false;
    for (auto &gar : gar_list) {
      MatchResultList match_result;
      int res = GARMatch<false>(gar, data_graph, match_result);
      if (res < 0) return res;

      auto &y_literal_set = gar.y_literal_set();
      for (auto &single_literal_ptr : y_literal_set) {
        for (auto &single_match_result : match_result) {
          if (single_literal_ptr->Satisfy(single_match_result)) {
            continue;
          }
          bool r = single_literal_ptr->Update(single_match_result, data_graph,
                                              edge_id_gen, diff_vertex_set,
                                              diff_edge_set);

          update_flag |= r;
          if (r) ++count;
        }
      }
    }
  } while (update_flag);

  return count;
}

}      // namespace gar

#endif  // GAR_CHASE_H_