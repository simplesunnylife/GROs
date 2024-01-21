#ifndef EXAMPLES_ANALYTICAL_APPS_TOOLS_INFO_TOOLS_INFO_H_
#define EXAMPLES_ANALYTICAL_APPS_TOOLS_INFO_TOOLS_INFO_H_

#include <vector>
#include <string>

namespace tools {

enum class ActionType : uint8_t {
  kNull,
  kGarSetToGar,
  kGarToGarSet,
  kGarInfoCollection,
  kGraphStatisticCollection
};

struct ToolInfo {

  void Reset() {
    action_ = ActionType::kNull;
    // for GarSetToGar
    gar_set_dir_.clear();
    gar_set_prefix_set_.clear();
    gar_size_in_each_set_file_ = 0;
    collect_vertex_num_ = false;
    collect_edge_num_ = false;
    collect_is_link_ = false;
    collect_is_star_ = false;
    x_literal_num_ = false;
    y_literal_num_ = false;
    
    // for GraphStatisticCollection
    collect_vertex_label_ = false;  
    collect_edge_label_ = false;
    collect_edge_type_ = false;
    collect_vertex_attribute_ = false;
    return;
  }

  enum ActionType action_;

  // for GarSetToGar
  std::string gar_set_dir_;
  std::vector<std::string> gar_set_prefix_set_;

  int gar_size_in_each_set_file_;

  // for GarInfoCollection
  std::string output_file_;
  bool collect_vertex_num_,
         collect_edge_num_,
          collect_is_link_,
          collect_is_star_,
            x_literal_num_,
            y_literal_num_;

  // for GraphStatisticCollection
  std::string v_file_, e_file_;
  bool collect_vertex_label_,
       collect_edge_label_,
       collect_edge_type_,
       collect_vertex_attribute_;
};

}; // tools

#endif // EXAMPLES_ANALYTICAL_APPS_TOOLS_INFO_TOOLS_INFO_H_