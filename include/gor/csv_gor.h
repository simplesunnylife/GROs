#ifndef GOR_CSV_GOR_H_
#define GOR_CSV_GOR_H_

#include "gor/gor.h"

#include "gar/gar.h"
#include "gar/csv_gar.h"

#include "gundam/io/csvgraph.h"

namespace gor {

template <typename   PatternType, 
          typename DataGraphType>
int ReadGORSet(
    std::vector<GraphOracleRule<PatternType, 
                              DataGraphType>> &gor_set,
    std::vector<std::string> &gor_name_set, 
    const std::string &v_set_file, const std::string &e_set_file, 
    const std::string &x_set_file, const std::string &y_set_file) {

  std::vector<GraphOracleRule<PatternType, 
                            DataGraphType>> gar_set;

  auto res = gar::ReadGARSet(gar_set, gor_name_set, 
                            v_set_file, e_set_file, 
                            x_set_file, y_set_file);  
  if (res < 0)
    return res;

  gor_set.reserve(gar_set.size());
  for (const auto& gar : gar_set) {
    gor_set.emplace_back(gar);
  }

  return res;
}

template <typename   PatternType, 
          typename DataGraphType>
int ReadGORSet(
    std::vector<GraphOracleRule<PatternType, 
                              DataGraphType>> &gor_set,
    const std::string &v_set_file, const std::string &e_set_file,
    const std::string &x_set_file, const std::string &y_set_file) {
  std::vector<std::string> gor_name_set;
  return ReadGORSet(gor_set, 
                    gor_name_set, 
                    v_set_file, 
                    e_set_file, 
                    x_set_file,
                    y_set_file);
}

template <typename   PatternType, 
          typename DataGraphType>
int ReadGOR(GraphOracleRule<PatternType, 
                          DataGraphType> &gor,
            const std::string &v_file, const std::string &e_file,
            const std::string &x_file, const std::string &y_file) {
  return gar::ReadGAR(gor, v_file, e_file,
                           x_file, y_file);
}

template <typename   PatternType, 
          typename DataGraphType>
int WriteGORSet(
    const std::vector<GraphOracleRule<PatternType, DataGraphType>> &gor_set,
    const std::string &v_file, const std::string &e_file,
    const std::string &x_file, const std::string &y_file) {
  std::vector<std::string> gor_name_set;
  gor_name_set.reserve(gor_set.size());
  for (size_t i = 0; i < gor_set.size(); i++) {
    gor_name_set.emplace_back(std::to_string(i));
  }

  return WriteGORSet(gor_set, gor_name_set, v_file, e_file, x_file, y_file);
}

template <typename   PatternType, 
          typename DataGraphType>
int WriteGORSet(
    const std::vector<GraphOracleRule<PatternType, 
                                    DataGraphType>> &gor_set,
    const std::vector<std::string> &gor_name_list, 
    const std::string &v_file, const std::string &e_file, 
    const std::string &x_file, const std::string &y_file) {

  std::vector<GraphOracleRule<PatternType, 
                            DataGraphType>> gar_set;

  gar_set.reserve(gor_set.size());
  for (const auto& gor : gor_set) {
    gar_set.emplace_back(gor);
  }

  return gar::WriteGARSet(gar_set, gor_name_list,
                                  v_file, e_file,
                                  x_file, y_file);
}

template <typename   PatternType,
          typename DataGraphType>
int WriteGOR(const GraphOracleRule<PatternType, 
                                  DataGraphType> &gor,
             const std::string &v_file, const std::string &e_file,
             const std::string &x_file, const std::string &y_file) {
  return gar::WriteGAR(gor, v_file, e_file,
                            x_file, y_file);
}

}  // namespace gor

#endif