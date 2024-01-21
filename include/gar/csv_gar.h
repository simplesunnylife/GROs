#ifndef EXAMPLES_ANALYTICAL_CSVGAR_H_
#define EXAMPLES_ANALYTICAL_CSVGAR_H_

#include "gar.h"
#include "gundam/io/csvgraph.h"

namespace gar {
template <bool is_x, typename PatternType, typename DataGraphType>
int ReadCSVEdgeLiteral(GraphAssociationRule<PatternType, DataGraphType> &gar,
                       rapidcsv::Document &literal_csv,
                       size_t row_pos) noexcept {
  using PatternVertexID = typename PatternType::VertexType::IDType;
  using EdgeLabel = typename DataGraphType::EdgeType::LabelType;

  try {
    PatternVertexID src_id = literal_csv.GetCell<PatternVertexID>(1, row_pos);
    PatternVertexID dst_id = literal_csv.GetCell<PatternVertexID>(3, row_pos);
    EdgeLabel edge_label = literal_csv.GetCell<EdgeLabel>(5, row_pos);

    if (is_x) {
      gar.template AddX<EdgeLiteral<PatternType, DataGraphType>>(src_id, dst_id,
                                                                 edge_label);
    } else {
      gar.template AddY<EdgeLiteral<PatternType, DataGraphType>>(src_id, dst_id,
                                                                 edge_label);
    }
    return 1;

  } catch (...) {
    return 0;
  }
}

template <bool is_x, 
          typename   PatternType, 
          typename DataGraphType>
int ReadCSVMlLiteral(GraphAssociationRule<PatternType, 
                                        DataGraphType> &gar,
                       rapidcsv::Document &literal_csv,
                       size_t row_pos) noexcept {
  using PatternVertexID = typename PatternType::VertexType::IDType;
  using EdgeLabel = typename DataGraphType::EdgeType::LabelType;

  try {
    PatternVertexID src_id = literal_csv.GetCell<PatternVertexID>(1, row_pos);
    PatternVertexID dst_id = literal_csv.GetCell<PatternVertexID>(3, row_pos);
    EdgeLabel edge_label = literal_csv.GetCell<EdgeLabel>(5, row_pos);

    std::string module_url = "#";

    if (is_x) {
      gar.template AddX<MlLiteral<PatternType, DataGraphType>>(src_id, dst_id,
                                                               edge_label,
                                                               module_url);
    } else {
      gar.template AddY<MlLiteral<PatternType, DataGraphType>>(src_id, dst_id,
                                                               edge_label,
                                                               module_url);
    }
    return 1;

  } catch (...) {
    return 0;
  }
}

template <bool is_x, typename PatternType, typename DataGraphType>
int ReadCSVAttributeLiteral(
    GraphAssociationRule<PatternType, DataGraphType> &gar,
    rapidcsv::Document &literal_csv, size_t row_pos) noexcept {
  using PatternVertexID = typename PatternType::VertexType::IDType;
  using VertexAttributeKey =
      typename DataGraphType::VertexType::AttributeKeyType;

  try {
    PatternVertexID x_id = literal_csv.GetCell<PatternVertexID>(1, row_pos);
    VertexAttributeKey x_attr_key =
        literal_csv.GetCell<VertexAttributeKey>(2, row_pos);

    if (is_x) {
      gar.template AddX<AttributeLiteral<PatternType, DataGraphType>>(
          x_id, x_attr_key);
    } else {
      gar.template AddY<AttributeLiteral<PatternType, DataGraphType>>(
          x_id, x_attr_key);
    }
    return 1;

  } catch (...) {
    return 0;
  }
}

template <bool is_x, typename PatternType, typename DataGraphType>
int ReadCSVVariableLiteral(
    GraphAssociationRule<PatternType, DataGraphType> &gar,
    rapidcsv::Document &literal_csv, size_t row_pos) noexcept {
  using PatternVertexID = typename PatternType::VertexType::IDType;
  using VertexAttributeKey =
      typename DataGraphType::VertexType::AttributeKeyType;

  try {
    PatternVertexID x_id = literal_csv.GetCell<PatternVertexID>(1, row_pos);
    VertexAttributeKey x_attr_key =
        literal_csv.GetCell<VertexAttributeKey>(2, row_pos);
    PatternVertexID y_id = literal_csv.GetCell<PatternVertexID>(3, row_pos);
    VertexAttributeKey y_attr_key =
        literal_csv.GetCell<VertexAttributeKey>(4, row_pos);

    if (is_x) {
      gar.template AddX<VariableLiteral<PatternType, DataGraphType>>(
          x_id, x_attr_key, y_id, y_attr_key);
    } else {
      gar.template AddY<VariableLiteral<PatternType, DataGraphType>>(
          x_id, x_attr_key, y_id, y_attr_key);
    }
    return 1;

  } catch (...) {
    return 0;
  }
}

template <bool is_x, typename PatternType, typename DataGraphType>
int ReadCSVConstantLiteral(
    GraphAssociationRule<PatternType, DataGraphType> &gar,
    rapidcsv::Document &literal_csv, size_t row_pos) noexcept {
  using PatternVertexID = typename PatternType::VertexType::IDType;
  using VertexAttributeKey =
      typename DataGraphType::VertexType::AttributeKeyType;

  try {
    PatternVertexID x_id = literal_csv.GetCell<PatternVertexID>(1, row_pos);
    VertexAttributeKey x_attr_key =
        literal_csv.GetCell<VertexAttributeKey>(2, row_pos);
    std::string value_with_type = literal_csv.GetCell<std::string>(6, row_pos);

    std::string value_str, type_str;
    bool partition_flag = false;
    for (const auto &c : value_with_type) {
      if (c == ';') {
        partition_flag = true;
        continue;
      }
      if (partition_flag)
        type_str.push_back(c);
      else
        value_str.push_back(c);
    }

    GUNDAM::BasicDataType value_type = GUNDAM::StringToEnum(type_str.c_str());

    switch (value_type) {
      case GUNDAM::BasicDataType::kTypeString: {
        if (is_x)
          gar.template AddX<ConstantLiteral<PatternType, DataGraphType>>(
              x_id, x_attr_key, value_str);
        else
          gar.template AddY<ConstantLiteral<PatternType, DataGraphType>>(
              x_id, x_attr_key, value_str);
        break;
      }
      case GUNDAM::BasicDataType::kTypeInt: {
        if (is_x)
          gar.template AddX<ConstantLiteral<PatternType, DataGraphType>>(
              x_id, x_attr_key, std::stoi(value_str));
        else
          gar.template AddY<ConstantLiteral<PatternType, DataGraphType>>(
              x_id, x_attr_key, std::stoi(value_str));
        break;
      }
      case GUNDAM::BasicDataType::kTypeInt64: {
        if (is_x)
          gar.template AddX<ConstantLiteral<PatternType, DataGraphType>>(
              x_id, x_attr_key, std::stoi(value_str));
        else
          gar.template AddY<ConstantLiteral<PatternType, DataGraphType>>(
              x_id, x_attr_key, std::stoi(value_str));
        break;
      }
      case GUNDAM::BasicDataType::kTypeFloat: {
        if (is_x)
          gar.template AddX<ConstantLiteral<PatternType, DataGraphType>>(
              x_id, x_attr_key, std::stof(value_str));
        else
          gar.template AddY<ConstantLiteral<PatternType, DataGraphType>>(
              x_id, x_attr_key, std::stof(value_str));
        break;
      }
      case GUNDAM::BasicDataType::kTypeDouble: {
        if (is_x)
          gar.template AddX<ConstantLiteral<PatternType, DataGraphType>>(
              x_id, x_attr_key, std::stod(value_str));
        else
          gar.template AddY<ConstantLiteral<PatternType, DataGraphType>>(
              x_id, x_attr_key, std::stod(value_str));
        break;
      }
      case GUNDAM::BasicDataType::kTypeDateTime: {
        if (is_x)
          gar.template AddX<ConstantLiteral<PatternType, DataGraphType>>(
              x_id, x_attr_key, GUNDAM::DateTime(value_str));
        else
          gar.template AddY<ConstantLiteral<PatternType, DataGraphType>>(
              x_id, x_attr_key, GUNDAM::DateTime(value_str));
        break;
      }
      case GUNDAM::BasicDataType::kTypeUnknown:
      default:
        return -1;
    }
    return 1;

  } catch (...) {
    return 0;
  }
}
template <bool is_x, 
          typename   PatternType, 
          typename DataGraphType>
int ReadCSVLiteralSetFile(
    std::vector<GraphAssociationRule<PatternType, 
                                   DataGraphType>> &gar_set,
    const std::string &literal_file) {
  try {
    std::cout << literal_file << std::endl;

    rapidcsv::Document literal_csv(literal_file, rapidcsv::LabelParams(0, -1));

    int count = 0;

    size_t row_count = literal_csv.GetRowCount();
    for (size_t row = 0; row < row_count; row++) {
      std::string literal_type = literal_csv.GetCell<std::string>(0, row);
      int gar_id = literal_csv.GetCell<int>(7, row);
      int res;
      if (literal_type == "Edge") {
        res = ReadCSVEdgeLiteral<is_x>(gar_set[gar_id], literal_csv, row);
      } else if (literal_type == "Attribute") {
        res = ReadCSVAttributeLiteral<is_x>(gar_set[gar_id], literal_csv, row);
      } else if (literal_type == "Variable") {
        res = ReadCSVVariableLiteral<is_x>(gar_set[gar_id], literal_csv, row);
      } else if (literal_type == "Constant") {
        res = ReadCSVConstantLiteral<is_x>(gar_set[gar_id], literal_csv, row);
      } else if (literal_type == "Ml") {
        res = ReadCSVMlLiteral<is_x>(gar_set[gar_id], literal_csv, row);
      } else {
        return -1;
      }
      if (res < 0) return res;

      ++count;
    }

    std::cout << "Literal: " << count << std::endl;

    return count;
  } catch (...) {
    return -1;
  }
}

template <bool is_x, typename PatternType, typename DataGraphType>
int ReadCSVLiteralFile(GraphAssociationRule<PatternType, DataGraphType> &gar,
                       const std::string &literal_file) {
  try {
    std::cout << literal_file << std::endl;

    rapidcsv::Document literal_csv(literal_file, rapidcsv::LabelParams(0, -1));

    int count = 0;

    size_t row_count = literal_csv.GetRowCount();
    for (size_t row = 0; row < row_count; row++) {
      std::string literal_type = literal_csv.GetCell<std::string>(0, row);

      int res;
      if (literal_type == "Edge") {
        res = ReadCSVEdgeLiteral<is_x>(gar, literal_csv, row);
      } else if (literal_type == "Attribute") {
        res = ReadCSVAttributeLiteral<is_x>(gar, literal_csv, row);
      } else if (literal_type == "Variable") {
        res = ReadCSVVariableLiteral<is_x>(gar, literal_csv, row);
      } else if (literal_type == "Constant") {
        res = ReadCSVConstantLiteral<is_x>(gar, literal_csv, row);
      } else if (literal_type == "Ml") {
        res = ReadCSVMlLiteral<is_x>(gar, literal_csv, row);
      } else {
        return -1;
      }
      if (res < 0) return res;

      ++count;
    }

    std::cout << "Literal: " << count << std::endl;

    return count;
  } catch (...) {
    return -1;
  }
}

template <typename   PatternType, 
          typename DataGraphType>
int ReadGARSet(
    std::vector<GraphAssociationRule<PatternType, 
                                   DataGraphType>> &gar_set,
    std::vector<std::string> &gar_name_set, 
    const std::string &v_set_file,
    const std::string &e_set_file, 
    const std::string &x_set_file,
    const std::string &y_set_file) {
  int res;
  gar_set.clear();
  gar_name_set.clear();
  int count = 0;

  std::vector<PatternType> pattern_list;
  res = GUNDAM::ReadCSVGraphSet(pattern_list, gar_name_set, v_set_file,
                                e_set_file);
  if (res < 0) return res;
  count += res;

  gar_set.reserve(pattern_list.size());
  const size_t kNumberOfGar = pattern_list.size();
  for (size_t i = 0; i < kNumberOfGar; i++) {
    gar_set.emplace_back(pattern_list[i]);
  }
  res = ReadCSVLiteralSetFile<true>(gar_set, x_set_file);
  if (res < 0) return res;
  count += res;

  res = ReadCSVLiteralSetFile<false>(gar_set, y_set_file);
  if (res < 0) return res;
  count += res;

  return count;
}

template <typename PatternType, typename DataGraphType>
int ReadGARSet(
    std::vector<GraphAssociationRule<PatternType, DataGraphType>> &gar_set,
    const std::string &v_set_file, const std::string &e_set_file,
    const std::string &x_set_file, const std::string &y_set_file) {
  std::vector<std::string> gar_name_set;
  return ReadGARSet(gar_set, gar_name_set, v_set_file, e_set_file, x_set_file,
                    y_set_file);
}

template <typename PatternType, typename DataGraphType>
int ReadGAR(GraphAssociationRule<PatternType, DataGraphType> &gar,
            const std::string &v_file, const std::string &e_file,
            const std::string &x_file, const std::string &y_file) {
  int res;

  gar.Reset();

  int count = 0;

  res = GUNDAM::ReadCSVGraph(gar.pattern(), v_file, e_file);
  if (res < 0) return res;
  count += res;

  res = ReadCSVLiteralFile<true>(gar, x_file);
  if (res < 0) return res;
  count += res;

  res = ReadCSVLiteralFile<false>(gar, y_file);
  if (res < 0) return res;
  count += res;

  return count;
}

constexpr char kCSVLiteralHead[] = "type,x_id,x_attr,y_id,y_attr,edge_label,c";
constexpr char kCSVLiteralSetHead[] =
    "type,x_id,x_attr,y_id,y_attr,edge_label,c,gar_id";

template <bool is_x, typename   PatternType, 
                     typename DataGraphType>
int WriteCSVLiteralSetFile(
    const std::vector<GraphAssociationRule<PatternType, DataGraphType>>
        &gar_set_list,
    const std::vector<std::string> &gar_name_set,
    const std::string &literal_file) {
  std::ofstream f(literal_file);
  if (!f) return -1;

  f << kCSVLiteralSetHead << std::endl;
  int count = 0;

  for (size_t gar_idx = 0; gar_idx < gar_set_list.size(); gar_idx++) {
    auto &gar = gar_set_list[gar_idx];
    if (gar_name_set[gar_idx] == "") {
      std::cout << "gar name is empty" << std::endl;
      return -1;
    }
    if (is_x) {
      for (const auto &l : gar.x_literal_set()) {
        l->Print(f, gar_name_set[gar_idx]);
      }
    } else {
      for (const auto &l : gar.y_literal_set()) {
        l->Print(f, gar_name_set[gar_idx]);
      }
    }
    ++count;
  }

  return count;
}

template <class LiteralSet>
int WriteCSVLiteralFile(const LiteralSet &literal_set,
                        const std::string &literal_file) {
  std::ofstream f(literal_file);
  if (!f) return -1;

  f << kCSVLiteralHead << std::endl;

  int count = 0;
  for (const auto &l : literal_set) {
    l->Print(f);
    ++count;
  }

  return count;
}

template <typename PatternType, typename DataGraphType>
int WriteGARSet(
    const std::vector<GraphAssociationRule<PatternType, DataGraphType>>
        &gar_set,
    const std::string &v_file, const std::string &e_file,
    const std::string &x_file, const std::string &y_file) {
  std::vector<std::string> gar_name_set;
  gar_name_set.reserve(gar_set.size());
  for (size_t i = 0; i < gar_set.size(); i++) {
    gar_name_set.emplace_back(std::to_string(i));
  }

  return WriteGARSet(gar_set, gar_name_set, v_file, e_file, x_file, y_file);
}

template <typename   PatternType, 
          typename DataGraphType>
int WriteGARSet(
    const std::vector<GraphAssociationRule<PatternType, 
                                         DataGraphType>> &gar_set,
    const std::vector<std::string> &gar_name_list, 
    const std::string &v_file, const std::string &e_file, 
    const std::string &x_file, const std::string &y_file) {

  int res;

  if (gar_set.size() != gar_name_list.size()) {
    std::cout << "gar_set.size() and gar_name_list.size() mismatch!"
              << std::endl;
    return -1;
  }

  int count = 0;
  std::vector<PatternType> pattern_list;
  for (auto &gar : gar_set) {
    pattern_list.emplace_back(gar.pattern());
  }
  res = GUNDAM::WriteCSVGraphSet<false>(pattern_list, gar_name_list, v_file,
                                        e_file);
  if (res < 0) return res;
  count += res;

  res = WriteCSVLiteralSetFile<true>(gar_set, gar_name_list, x_file);
  if (res < 0) return res;
  count += res;

  res = WriteCSVLiteralSetFile<false>(gar_set, gar_name_list, y_file);
  if (res < 0) return res;
  count += res;

  return count;
}

template <typename   PatternType,
          typename DataGraphType>
int WriteGAR(const GraphAssociationRule<PatternType, 
                                      DataGraphType> &gar,
             const std::string &v_file, const std::string &e_file,
             const std::string &x_file, const std::string &y_file) {
  int res;

  int count = 0;
  res = GUNDAM::WriteCSVGraph<false>(gar.pattern(), v_file, e_file);
  if (res < 0) return res;
  count += res;

  res = WriteCSVLiteralFile(gar.x_literal_set(), x_file);
  if (res < 0) return res;
  count += res;

  res = WriteCSVLiteralFile(gar.y_literal_set(), y_file);
  if (res < 0) return res;
  count += res;

  return count;
}

}  // namespace gar

#endif