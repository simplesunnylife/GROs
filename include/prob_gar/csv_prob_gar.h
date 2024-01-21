#ifndef _PROB_GAR_CSV_PROB_GAR_H
#define _PROB_GAR_CSV_PROB_GAR_H

#include "include/gundam/io/csvgraph.h"

#include "include/gar/gar.h"
#include "include/gar/csv_gar.h"

#include "include/prob_gar/prob_gar.h"

namespace prob_gar {

template <typename   PatternType, 
          typename DataGraphType>
int ReadProbGARSet(
    std::vector<ProbGraphAssociationRule<PatternType, 
                                       DataGraphType>> &prob_gar_set,
    std::vector<std::string> &prob_gar_name_set, 
    const std::string &v_set_file,
    const std::string &e_set_file, 
    const std::string &x_set_file,
    const std::string &y_set_file,
    const std::string &p_set_file) {

  int res;
  prob_gar_set.clear();
  prob_gar_name_set.clear();
  int count = 0;

  std::vector<PatternType> pattern_list;
  res = GUNDAM::ReadCSVGraphSet(pattern_list, 
                           prob_gar_name_set, 
                                  v_set_file,
                                  e_set_file);
  if (res < 0) 
    return res;
    
  count += res;

  std::vector<gar::GraphAssociationRule<PatternType, 
                                      DataGraphType>> gar_set;

  gar_set.reserve(pattern_list.size());
  const size_t kNumberOfGar = pattern_list.size();
  for (size_t i = 0; i < kNumberOfGar; i++) {
    gar_set.emplace_back(pattern_list[i]);
  }

  res = gar::ReadCSVLiteralSetFile<true>(gar_set, x_set_file);
  if (res < 0) 
    return res;
  count += res;

  res = gar::ReadCSVLiteralSetFile<false>(gar_set, y_set_file);
  if (res < 0) 
    return res;
  count += res;

  prob_gar_set.reserve(gar_set.size());

  for (const auto& gar : gar_set) {
    prob_gar_set.emplace_back(gar);
  }

  try {
    std::cout << p_set_file << std::endl;

    rapidcsv::Document p_csv(p_set_file, rapidcsv::LabelParams(0, -1));

    int count = 0;

    size_t row_count = p_csv.GetRowCount();
    for (size_t row = 0; row < row_count; row++) {
      float prob = p_csv.GetCell<float>(0, row);
      int gar_id = p_csv.GetCell< int >(1, row);
      if (gar_id < 0 || gar_id > prob_gar_set.size()) {
        return -1;
      }
      prob_gar_set[gar_id].set_prob(prob);
      ++count;
    }
  }
  catch (...) {
    return -1;
  }

  return count;
}

template <typename   PatternType, 
          typename DataGraphType>
int ReadProbGARSet(
    std::vector<ProbGraphAssociationRule<PatternType, 
                                       DataGraphType>> &prob_gar_set,
    const std::string &v_set_file, 
    const std::string &e_set_file,
    const std::string &x_set_file, 
    const std::string &y_set_file, 
    const std::string &p_set_file) {
  std::vector<std::string> prob_gar_name_set;
  return ReadProbGARSet(prob_gar_set, prob_gar_name_set, 
                                 v_set_file, e_set_file, 
                                 x_set_file, y_set_file,
                                 p_set_file);
}

template <typename   PatternType, 
          typename DataGraphType>
int ReadProbGAR(ProbGraphAssociationRule<PatternType, 
                                       DataGraphType> &prob_gar,
            const std::string &v_file, 
            const std::string &e_file,
            const std::string &x_file, 
            const std::string &y_file, 
            const std::string &p_file) {

  prob_gar.Reset();
  
  int res = gar::ReadGAR(prob_gar, v_file, e_file,
                                   x_file, y_file);

  if (res < 0) 
    return res;

  std::cout << p_file << std::endl;

  try {
    std::cout << p_file << std::endl;

    rapidcsv::Document p_csv(p_file, rapidcsv::LabelParams(0, -1));

    int count = 0;

    if (p_csv.GetRowCount() != 1) {
      return -1;
    }

    const float kProb = p_csv.GetCell<float>(0, 0);

    prob_gar.set_prob(kProb);
  }
  catch (...) {
    return -1;
  }
  return res + 1;
}

constexpr char kProbSetHead[] = "prob,gar_id";

template <typename   PatternType, 
          typename DataGraphType>
int WriteProbSetFile(
    const std::vector<ProbGraphAssociationRule<PatternType, 
                                             DataGraphType>> &prob_gar_set_list,
    const std::vector<std::string> &gar_name_set,
    const std::string &prob_file) {

  assert(gar_name_set.size() == prob_gar_set_list.size());

  std::ofstream f(prob_file);
  if (!f) 
    return -1;

  f << kProbSetHead << std::endl;
  int count = 0;

  for (size_t gar_idx = 0; gar_idx < prob_gar_set_list.size(); gar_idx++) {
    auto &gar = prob_gar_set_list[gar_idx];
    if (gar_name_set[gar_idx] == "") {
      std::cout << "gar name is empty" << std::endl;
      return -1;
    }

    f << gar.prob() << "," << gar_name_set[gar_idx] << std::endl;

    ++count;
  }

  return count;
}

template <typename   PatternType, 
          typename DataGraphType>
int WriteProbGARSet(
    const std::vector<ProbGraphAssociationRule<PatternType, 
                                             DataGraphType>> &prob_gar_set,
    const std::vector<std::string> &prob_gar_name_list, 
    const std::string &v_file, const std::string &e_file, 
    const std::string &x_file, const std::string &y_file, 
    const std::string &p_file) {

  int res;

  if (prob_gar_set.size() != prob_gar_name_list.size()) {
    std::cout << "prob_gar_set.size() and prob_gar_name_list.size() mismatch!"
              << std::endl;
    return -1;
  }

  int count = 0;
  std::vector<PatternType> pattern_list;
  for (auto &gar : prob_gar_set) {
    pattern_list.emplace_back(gar.pattern());
  }
  res = GUNDAM::WriteCSVGraphSet<false>(pattern_list, 
                                  prob_gar_name_list, 
                                      v_file, e_file);
  if (res < 0) 
    return res;
  count += res;

  std::vector<gar::GraphAssociationRule<PatternType, 
                                      DataGraphType>> gar_list;
  for (auto &gar : prob_gar_set) {
    gar_list.emplace_back(gar);
  }

  res = gar::WriteCSVLiteralSetFile<true>(gar_list, 
                                 prob_gar_name_list, x_file);
  if (res < 0) 
    return res;
  count += res;

  res = gar::WriteCSVLiteralSetFile<false>(gar_list, 
                                 prob_gar_name_list, y_file);
  if (res < 0) 
    return res;
  count += res;

  res = WriteProbSetFile(prob_gar_set, 
                         prob_gar_name_list, p_file);
  if (res < 0) 
    return res;
  count += res;

  return count;
}

template <typename   PatternType, 
          typename DataGraphType>
int WriteProbGARSet(
    const std::vector<ProbGraphAssociationRule<PatternType, 
                                             DataGraphType>> &prob_gar_set,
    const std::string &v_file, const std::string &e_file,
    const std::string &x_file, const std::string &y_file, 
    const std::string &p_file) {
  std::vector<std::string> prob_gar_name_set;
  prob_gar_name_set.reserve(prob_gar_set.size());

  for (size_t i = 0; i < prob_gar_set.size(); i++) {
    prob_gar_name_set.emplace_back(std::to_string(i));
  }

  return WriteProbGARSet(prob_gar_set, 
                         prob_gar_name_set, 
                            v_file, e_file, 
                            x_file, y_file, 
                            p_file);
}

template <typename   PatternType,
          typename DataGraphType>
int WriteProbGAR(const ProbGraphAssociationRule<PatternType, 
                                               DataGraphType> &prob_gar,
                 const std::string &v_file, 
                 const std::string &e_file,
                 const std::string &x_file, 
                 const std::string &y_file, 
                 const std::string &p_file) {

  int res = WriteGAR(prob_gar, v_file, e_file,
                               x_file, y_file);
  if (res < 0) 
    return res;

  std::ofstream p_fstream(p_file);

  if(!p_fstream) { 
    // file can not be opened
    return -1;
  }

  p_fstream << prob_gar.prob() << std::endl;

  return res + 1;
}

}  // namespace prob_gar

#endif // _PROB_GAR_CSV_PROB_GAR_H