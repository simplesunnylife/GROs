#ifndef GAR_IMPLY_H_
#define GAR_IMPLY_H_

#include <map>

#include "gar/gar.h"

#include "gundam/component/generator.h"

#include "gundam/type_getter/vertex_handle.h"
#include "gundam/type_getter/edge_handle.h"

namespace gar{

// ###############################################################
// #  wenzhi: to do: modified this interface to:                 #
// #                                                             #
// #  template <typename   SrcPatternType,                       #
// #            typename SrcDataGraphType,                       #
// #            typename   DstPatternType,                       #
// #            typename DstDataGraphType>                       #
// #  bool GarImply(                                             #
// #      std::list<                                             #
// #        GraphAssociationRule<SrcPatternType,                 #
// #                           SrcDataGraphType>>& src_gar_list, #
// #        GraphAssociationRule<DstPatternType,                 #
// #                           DstDataGraphType>&  dst_gar);     #
// #                                                             #
// #  for better generalization                                  #
// ###############################################################
//
// return whether dst_gar can be implied from src_gar_list
template <typename   PatternType,
          typename DataGraphType>
bool GarImply(
    std::list<
      GraphAssociationRule<PatternType,
                         DataGraphType>>& src_gar_list, 
      GraphAssociationRule<PatternType,
                         DataGraphType>&  dst_gar) {

  using   PatternVertexConstPtr = typename GUNDAM::VertexHandle<  PatternType>::type;
  using DataGraphVertexConstPtr = typename GUNDAM::VertexHandle<DataGraphType>::type;

  using MatchMap = std::map<PatternVertexConstPtr, 
                          DataGraphVertexConstPtr>;

  using MatchResultList = std::vector<MatchMap>;

  using EdgeIDType = typename PatternType::EdgeType::IDType;

  PatternType dst_pattern (dst_gar.pattern());

  MatchMap dst_gar_to_dst_pattern;
  // add one-to-one match from dst_gar.pattern() to dst_pattern
  for (auto vertex_cit = dst_gar.pattern().VertexBegin();
           !vertex_cit.IsDone();
            vertex_cit++) {
    auto dst_ptr = dst_pattern.FindVertex(vertex_cit->id());
    assert(dst_ptr);
    dst_gar_to_dst_pattern.emplace(vertex_cit, dst_ptr);
  }

  GUNDAM::SimpleArithmeticIDGenerator<typename PatternType::EdgeType::IDType>
    edge_id_gen;

  EdgeIDType max_edge_id = 0;
  for (auto vertex_it = dst_pattern.VertexBegin(); 
           !vertex_it.IsDone();
            vertex_it++){
    for (auto edge_it = vertex_it->OutEdgeBegin();
             !edge_it.IsDone();
              edge_it++){
      if (max_edge_id < edge_it->id()){
        max_edge_id = edge_it->id();
      }
    }
  }

  for (int i = 0; i <= max_edge_id; i++){
    edge_id_gen.GetID();
  }
  
  for (const auto& x_literal_ptr 
         : dst_gar.x_literal_set()) {
    if (x_literal_ptr->Satisfy(dst_gar_to_dst_pattern)) {
      continue;
    }
    x_literal_ptr->Update(dst_gar_to_dst_pattern,
                                     dst_pattern,
                                     edge_id_gen);
  }
  
  while (true) {
    int update_counter = 0;
    for (auto& gar_cit : src_gar_list) {

      MatchResultList match_result;
      int res = gar::GARMatch<false>(gar_cit, dst_pattern, 
                                              match_result);

      assert(res >= 0);

      const auto &y_literal_set = gar_cit.y_literal_set();
      for (const auto &y_literal_ptr : y_literal_set) {
        for (const auto &single_match_result : match_result) {
          if (y_literal_ptr->Satisfy(single_match_result)) {
            continue;
          }
          if (!y_literal_ptr->Update(single_match_result,
                                            dst_pattern,
                                            edge_id_gen)){
            // violate
            std::cout << "violate" << std::endl;
            return true;
          }
          update_counter++;
        }
      }
    }
    if (update_counter == 0){
      break;
    }
  }

  for (const auto& y_literal_ptr 
         : dst_gar.y_literal_set()) {
    if (!y_literal_ptr->Satisfy(dst_gar_to_dst_pattern)) {
      // does not satisfy all literals in Y
      // cannot be implied from src_gar_list
      std::cout << "dst_gar violate" << std::endl;
      return false;
    }
  }
  // satisfy all literals in Y
  return true;
}

} // gar

#endif // GAR_IMPLY_H_