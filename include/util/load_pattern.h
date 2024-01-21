#ifndef UTIL_LOAD_PATTERN_H_
#define UTIL_LOAD_PATTERN_H_

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "include/gundam/io/csvgraph.h"

namespace util {

template <typename PatternType>
inline std::vector<
       std::pair<PatternType,
                 std::string> > // ok to return in C11
       LoadPatternSet(const YAML::Node& pattern_path_config) {
  std::vector<std::pair<PatternType, std::string>> pattern_set;
  if (pattern_path_config["VFile"]) {
    // is a pattern 
    if (!pattern_path_config["EFile"]) {
      // configure does not complete
      return std::vector<std::pair<PatternType, std::string>>();
    }
    const std::string kVFile = pattern_path_config["VFile"].as<std::string>(),
                      kEFile = pattern_path_config["EFile"].as<std::string>();
            
    PatternType pattern;
    if (GUNDAM::ReadCSVGraph(pattern, kVFile, kEFile) < 0) {
      return std::vector<std::pair<PatternType, std::string>>();
    }
    std::string pattern_name = "<" + kVFile + ","
                                   + kEFile + ">";
    pattern_set.reserve(1);
    pattern_set.emplace_back(std::move(pattern), std::move(pattern_name));
    return pattern_set;
  }
  if (pattern_path_config["VSetFile"]) {
    // is a pattern set
    if (!pattern_path_config["ESetFile"]) {
      // configure does not complete
      return std::vector<std::pair<PatternType, std::string>>();
    }
    const std::string kVSetFile = pattern_path_config["VSetFile"].as<std::string>(),
                      kESetFile = pattern_path_config["ESetFile"].as<std::string>();
            
    std::vector<PatternType> temp_pattern_set;
    std::vector<std::string>      pattern_name;
    if (GUNDAM::ReadCSVGraphSet(temp_pattern_set, 
                                     pattern_name, kVSetFile, 
                                                   kESetFile) < 0) {
      return std::vector<std::pair<PatternType, std::string>>();
    }
    const std::string kPatternNamePrefix = "<" + kVSetFile + ","
                                               + kESetFile + ",";
    pattern_set.reserve(temp_pattern_set.size());
    for (size_t i = 0; i < temp_pattern_set.size(); i++) {
      pattern_set.emplace_back(std::move(temp_pattern_set [i]),
                         kPatternNamePrefix + pattern_name[i] + ">");
    }
    return pattern_set;
  }
  return std::vector<std::pair<PatternType, std::string>>();
}

}; // util

#endif // UTIL_LOAD_PATTERN_H_