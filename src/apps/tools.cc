
// #include "include/util/log.h"

#include "src/apps/tools/yaml_config.h"
#include "src/apps/tools/tools_info.h"
#include "src/apps/tools/gar_set_to_gar.h"
#include "src/apps/tools/gar_to_gar_set.h"
#include "src/apps/tools/gar_info_collection.h"
#include "src/apps/tools/graph_statistic_collection.h"

#include <string>

int main(int argc, char* argv[]) {
  if (argc > 2) {
    util::Error("parameter num is not correct! (0 or 1)");
    exit(-1);
  }

  std::string config_file_path = "tools_config.yaml";
  if (argc == 2) {
    config_file_path = argv[1];
  }

  YAML::Node config_node = YAML::LoadFile(config_file_path);

  tools::ToolInfo tool_info;

  bool res = tools::ReadGraphPackageInfoYAML(config_node, tool_info);
  if (!res) {
    return -1;
  }

  switch (tool_info.action_) {
   case tools::ActionType::kGarSetToGar :

    assert(!tool_info.gar_set_prefix_set_.empty());

    std::sort(tool_info.gar_set_prefix_set_.begin(),
              tool_info.gar_set_prefix_set_.end());
              
    tool_info.gar_set_prefix_set_.erase(
      std::unique(tool_info.gar_set_prefix_set_.begin(), 
                  tool_info.gar_set_prefix_set_.end()),
      tool_info.gar_set_prefix_set_.end()
    );

    tools::GarSetToGar(tool_info.gar_set_dir_,
                       tool_info.gar_set_prefix_set_,
                       tool_info.gar_size_in_each_set_file_);
    break;
    
   case tools::ActionType::kGarToGarSet :

    assert(!tool_info.gar_set_prefix_set_.empty());

    std::sort(tool_info.gar_set_prefix_set_.begin(),
              tool_info.gar_set_prefix_set_.end());
              
    tool_info.gar_set_prefix_set_.erase(
      std::unique(tool_info.gar_set_prefix_set_.begin(), 
                  tool_info.gar_set_prefix_set_.end()),
      tool_info.gar_set_prefix_set_.end()
    );

    tools::GarToGarSet(tool_info.gar_set_dir_,
                       tool_info.gar_set_prefix_set_);
    break;
    
   case tools::ActionType::kGarInfoCollection :
    tools::GarInfoCollection(tool_info);
    break;
    
   case tools::ActionType::kGraphStatisticCollection :
    tools::GraphStatisticCollection(tool_info);
    break;

   default:
    util::Error("illegal action");
    break;
  }

  return 0;
}
