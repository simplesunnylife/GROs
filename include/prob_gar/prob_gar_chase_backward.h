#ifndef _PROB_GAR_PROB_GAR_CHASE_BACKWARD_H
#define _PROB_GAR_PROB_GAR_CHASE_BACKWARD_H

#include "prob_gar.h"

#include "include/gundam/algorithm/match_using_match.h"
#include "include/gundam/match/match.h"

#include "include/gundam/type_getter/vertex_id.h"

namespace prob_gar {

// backward chase
template <typename GraphPatternType, 
          typename    DataGraphType>
float ProbGARChase(
    const std::vector<ProbGraphAssociationRule<GraphPatternType, 
                                                  DataGraphType>> &prob_gar_list, 
    const DataGraphType &data_graph,
    const typename GUNDAM:: VertexID<DataGraphType>::type& src_vertex_id,
    const typename GUNDAM::EdgeLabel<DataGraphType>::type& edge_label,
    const typename GUNDAM:: VertexID<DataGraphType>::type& dst_vertex_id) {

  assert(data_graph.FindVertex(src_vertex_id));
  assert(data_graph.FindVertex(dst_vertex_id));
    
  using VertexLabelType = GUNDAM::VertexLabel<DataGraphType>::type;
  using   EdgeLabelType = GUNDAM::  EdgeLabel<DataGraphType>::type;

  using EdgeType = std::tuple<VertexLabelType,
                                EdgeLabelType,
                              VertexLabelType>;

  EdgeType edge_type_to_chase(data_graph.FindVertex(src_vertex_id)->label(),
                              edge_label,
                              data_graph.FindVertex(dst_vertex_id)->label());

  while () {
    /* code */
  }
  
  return ;
}

}; // prob_gar

#endif // _PROB_GAR_PROB_GAR_CHASE_BACKWARD_H