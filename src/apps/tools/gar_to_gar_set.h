#ifndef EXAMPLES_ANALYTICAL_APPS_TOOLS_GAR_TO_GAR_SET_H_
#define EXAMPLES_ANALYTICAL_APPS_TOOLS_GAR_TO_GAR_SET_H_

#include "gundam/graph_type/large_graph.h"

#include "include/gar/gar.h"
#include "include/util/load_gar_path.h"

namespace tools {

size_t GarToGarSet(const std::string&  gar_dir,
       const std::vector<std::string>& gar_prefix_set) {

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

  size_t counter = 0;
                                        
  for (const auto& gar_prefix : gar_prefix_set){
    
    std::vector<GarType> gar_set;
    std::vector<std::string> gar_name_set;
    
    std::vector<std::string> gar_file_set 
       = util::LoadGarPath(gar_dir, gar_prefix);
    for (const auto& gar_file : gar_file_set){
      GarType gar;
      gar::ReadGAR(gar, gar_file + "_v.csv",
                        gar_file + "_e.csv",
                        gar_file + "_x.csv",
                        gar_file + "_y.csv");
      gar_set.emplace_back(std::move(gar));
      gar_name_set.emplace_back(gar_file);
    }

    gar::WriteGARSet(gar_set,  gar_name_set,
                    gar_dir + gar_prefix + "_v_set.csv", 
                    gar_dir + gar_prefix + "_e_set.csv", 
                    gar_dir + gar_prefix + "_x_set.csv", 
                    gar_dir + gar_prefix + "_y_set.csv");

    counter += gar_set.size();
  }

  return counter;
}

};


#endif // EXAMPLES_ANALYTICAL_APPS_TOOLS_GAR_TO_GAR_SET_H_