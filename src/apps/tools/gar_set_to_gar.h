#ifndef EXAMPLES_ANALYTICAL_APPS_TOOLS_GAR_SET_TO_GAR_H_
#define EXAMPLES_ANALYTICAL_APPS_TOOLS_GAR_SET_TO_GAR_H_

#include "gundam/graph_type/large_graph.h"

#include "include/gar/gar.h"
#include "include/util/load_gar_path.h"

namespace tools {

size_t GarSetToGar(const std::string&  gar_dir,
       const std::vector<std::string>& gar_prefix_set,
                             int gar_size_in_each_set_file = 0) {

  using VertexIDType    = int;
  using VertexLabelType = int;
  using VertexAttributeKeyType = std::string;
  using   EdgeIDType    = int;
  using   EdgeLabelType = int;
  using   EdgeAttributeKeyType = std::string;

  using DataGraph =
      GUNDAM::LargeGraph<VertexIDType, VertexLabelType, VertexAttributeKeyType, 
                           EdgeIDType,   EdgeLabelType,   EdgeAttributeKeyType>;

  using Pattern =
      GUNDAM::LargeGraph<VertexIDType, VertexLabelType, VertexAttributeKeyType, 
                           EdgeIDType,   EdgeLabelType,   EdgeAttributeKeyType>;

  using GarType = gar::GraphAssociationRule<Pattern,
                                          DataGraph>;

  assert(gar_size_in_each_set_file >= 0);

  size_t gar_counter = 0;

  for (const auto& gar_prefix : gar_prefix_set) {
    // if the prefix is empty, then there should only have this prefix
    assert(gar_prefix != ""
        || gar_prefix_set.size() == 1);
    std::vector<std::string> gar_file_set = util::LoadGarSetPath(gar_dir,
                                                                 gar_prefix);
    std::vector<GarType> gar_set;
    for (const auto& gar_file : gar_file_set) {
      std::vector<GarType> temp_gar_set;
      gar::ReadGARSet(temp_gar_set, gar_file + "_v_set.csv",
                                    gar_file + "_e_set.csv",
                                    gar_file + "_x_set.csv",
                                    gar_file + "_y_set.csv");
      gar_set.insert(gar_set.end(),
        std::make_move_iterator(temp_gar_set.begin()),
        std::make_move_iterator(temp_gar_set.end()));
    }
    gar_counter += gar_set.size();
    for (int i = 0; i < gar_set.size(); i++) {
      if (i >= gar_size_in_each_set_file
            && gar_size_in_each_set_file > 0 ) {
        break;
      }
      const auto& gar = gar_set[i];
      gar::WriteGAR(gar, 
                    gar_dir + gar_prefix + "_" + std::to_string(i) + "_v.csv", 
                    gar_dir + gar_prefix + "_" + std::to_string(i) + "_e.csv", 
                    gar_dir + gar_prefix + "_" + std::to_string(i) + "_x.csv", 
                    gar_dir + gar_prefix + "_" + std::to_string(i) + "_y.csv");
    }
  }
  return gar_counter;
}

};


#endif // EXAMPLES_ANALYTICAL_APPS_TOOLS_GAR_SET_TO_GAR_H_