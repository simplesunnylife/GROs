#ifndef UTIL_LOAD_GAR_PATH_H_
#define UTIL_LOAD_GAR_PATH_H_

#include "gundam/graph_type/large_graph.h"

#include "include/util/file.h"

#include "include/gar/gar.h"
#include "include/gar/csv_gar.h"

namespace util {

namespace _load_gar_path {

inline std::vector<std::string> // ok to return stl container after c11
  LoadGarPath(const std::string&  gar_dir,
              const std::string& rule_prefix,
              const std::string& rule_posfix) {

  assert(!rule_posfix.empty());

  std::set<std::string> dir_files;
  std::vector<std::string> qualified_gar_path_set_in_one_dir;

  util::GetFiles(gar_dir, dir_files);

  for (auto gar_file : dir_files){
    const bool isPrefix = (rule_prefix == "") // does not have require of the prefix of the rule
                        ||(rule_prefix.size() <= gar_file.size() 
         && std::mismatch(rule_prefix.begin(), 
                          rule_prefix.end(),
                             gar_file.begin(),    
                             gar_file.end()).first == rule_prefix.end());
    const bool isGarPosfix = rule_posfix.size() <= gar_file.size() 
            && std::mismatch(rule_posfix.begin(), 
                             rule_posfix. end (),
                                gar_file.end() - rule_posfix.size(), 
                                gar_file.end()).first == rule_posfix.end();
    if (!isPrefix || !isGarPosfix){
      continue;
    }

    // load seperated single gar
    gar_file = gar_file.substr(0, gar_file.length()
                              - rule_posfix.size());

    const std::string kGarName = gar_dir + "/" + gar_file;

    qualified_gar_path_set_in_one_dir.emplace_back(kGarName);
  }
  return qualified_gar_path_set_in_one_dir;
}

};

inline std::vector<std::string> // ok to return stl container after c11
  LoadGarPath(const std::string&  gar_dir,
              const std::string& rule_prefix) {
  return _load_gar_path::LoadGarPath(gar_dir,
                                    rule_prefix, "_v.csv");
}

inline std::vector<std::string> // ok to return stl container after c11
  LoadGarSetPath(const std::string&  gar_dir,
                 const std::string& rule_prefix) {
  return _load_gar_path::LoadGarPath(gar_dir,
                                    rule_prefix, "_v_set.csv");
}

};


#endif // UTIL_LOAD_GAR_PATH_H_