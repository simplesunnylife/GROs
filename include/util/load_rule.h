#ifndef UTIL_LOAD_RULE_H_
#define UTIL_LOAD_RULE_H_

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "include/gundam/io/csvgraph.h"
#include "gar/gar.h"

namespace util {

template <typename GraphPatternType,
          typename    DataGraphType>
inline std::vector<
       std::pair<gar::GraphAssociationRule<GraphPatternType,
                                              DataGraphType> ,
                 std::string> > // ok to return in C11
       LoadRule(const YAML::Node& rule_path_config) {
  using RuleType = gar::GraphAssociationRule<GraphPatternType,
                                                DataGraphType>;
  
  std::vector<std::pair<RuleType, std::string>> rule_set;

  if (rule_path_config["VFile"]) {
    // is single rule
    if (!rule_path_config["EFile"]
     || !rule_path_config["XFile"]
     || !rule_path_config["YFile"]) {
      // configure does not complete
      return std::vector<std::pair<RuleType, std::string>>();
    }
    const std::string kVFile = rule_path_config["VFile"].as<std::string>(),
                      kEFile = rule_path_config["EFile"].as<std::string>(),
                      kXFile = rule_path_config["XFile"].as<std::string>(),
                      kYFile = rule_path_config["YFile"].as<std::string>();
            
    RuleType rule;
    if (gar::ReadGAR(rule, kVFile, kEFile,
                           kXFile, kYFile) < 0) {
      return std::vector<std::pair<RuleType, std::string>>();
    }
    std::string rule_name = "<" + kVFile + ","
                                + kEFile + ","
                                + kXFile + ","
                                + kYFile + ">";
    rule_set.reserve(1);
    rule_set.emplace_back(std::move(rule), std::move(rule_name));
    return rule_set;
  }
  if (rule_path_config["VSetFile"]) {
    // is a rule set
    if (!rule_path_config["ESetFile"]
     || !rule_path_config["XSetFile"]
     || !rule_path_config["YSetFile"]) {
      // configure does not complete
      return std::vector<std::pair<RuleType, std::string>>();
    }
    const std::string kVSetFile = rule_path_config["VSetFile"].as<std::string>(),
                      kESetFile = rule_path_config["ESetFile"].as<std::string>(),
                      kXSetFile = rule_path_config["XSetFile"].as<std::string>(),
                      kYSetFile = rule_path_config["YSetFile"].as<std::string>();
            
    std::vector<  RuleType > temp_rule_set;
    std::vector<std::string>      rule_name;
    if (gar::ReadGARSet(temp_rule_set, 
                             rule_name, kVSetFile, kESetFile,
                                        kXSetFile, kYSetFile) < 0) {
      return std::vector<std::pair<RuleType, std::string>>();
    }
    const std::string kRuleNamePrefix = "<" + kVSetFile + ","
                                            + kESetFile + ","
                                            + kXSetFile + ","
                                            + kYSetFile + ","; 
    rule_set.reserve(temp_rule_set.size());
    for (size_t i = 0; i < temp_rule_set.size(); i++) {
      rule_set.emplace_back(std::move(temp_rule_set [i]),
                         kRuleNamePrefix + rule_name[i] + ">");
    }
    return rule_set;
  }
  return std::vector<std::pair<RuleType, std::string>>();
}

}; // util

#endif // UTIL_LOAD_RULE_H_