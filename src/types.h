
#ifndef TYPES_H_
#define TYPES_H_

#include <istream>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

#include "examples/gnn_sampler/util.h"
#include "grape/serialization/in_archive.h"
#include "grape/serialization/out_archive.h"
#include "grape/types.h"

namespace grape {

struct Data {
 public:
  using LabelType = int32_t;
  using AttributeKeyType = size_t;
  using AttributeValueType = std::string;
  using AttributeContainer = std::vector<std::string>;

 private:
  static constexpr LabelType kDefaultVertexLabel = 0;

 public:
  Data() { 
    this->attributes_.resize(0); 
    return;
  }

  Data(const EmptyType& empty_data) {
    this->attributes_.resize(0);
    this->label_ = kDefaultVertexLabel;
    return;
  }

  ~Data() { return; }

  void Clear() {
    this->attributes_.clear();
    return;
  }

  Data(int32_t label) : label_(label) { this->attributes_.resize(0); }
  Data(int32_t label, std::string attrs) : label_(label), attributes_{attrs} {
    return;
  }
  Data(int32_t label, std::vector<std::string> attrs)
      : label_(label), attributes_(attrs) {
    return;
  }

  Data(const Data& data) : label_(data.label_), attributes_(data.attributes_) {
    return;
  }

  int32_t label_;
  std::vector<std::string> attributes_;

  bool operator==(const Data& vdata) const {
    if (this->attributes_.size() != vdata.attributes_.size() ||
        this->label_ != vdata.label_)
      return false;
    for (size_t i = 0; i < this->attributes_.size(); i++) {
      if (this->attributes_[i] != vdata.attributes_[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator<(const Data& vdata) const {
    return (this->label_ < vdata.label_);
  }
};
inline InArchive& operator<<(InArchive& in_archive, const Data& data) {
  in_archive << data.label_;
  in_archive << data.attributes_;
  return in_archive;
}
inline OutArchive& operator>>(OutArchive& out_archive, Data& data) {
  out_archive >> data.label_;
  out_archive >> data.attributes_;
  return out_archive;
}
namespace internal {
template <>
inline const char* match<grape::Data>(char const* str, grape::Data& r,
                                      char const* end) {
  r.Clear();
  const char* match_end = match(str, r.label_, end);
  while (match_end != end) {
    std::string attr;
    match_end = match(match_end, attr, end);
    if (attr == "" || attr == "\n") break;
    r.attributes_.emplace_back(attr);
  }
  return match_end;
}

}  // namespace internal

struct EdgeData {
  using LabelType = int32_t;
  using AttributeKeyType = size_t;
  using AttributeValueType = std::string;
  using AttributeContainer = std::vector<std::string>;
  EdgeData() { this->attributes_.resize(0); }
  EdgeData(int32_t id, int32_t label) : id_(id), label_(label) {
    this->attributes_.resize(0);
  }
  EdgeData(int32_t id, int32_t label, std::vector<std::string> attrs)
      : id_(id), label_(label), attributes_(attrs) {}
  int id_;
  int32_t label_;
  std::vector<std::string> attributes_;

  bool operator==(const EdgeData& edata) const {
    if (this->id_ != edata.id_ ||
        this->attributes_.size() != edata.attributes_.size() ||
        this->label_ != edata.label_)
      return false;
    for (size_t i = 0; i < this->attributes_.size(); i++) {
      if (this->attributes_[i] != edata.attributes_[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator<(const EdgeData& edata) const {
    return (this->label_ < edata.label_);
  }
};
inline InArchive& operator<<(InArchive& in_archive, const EdgeData& edata) {
  in_archive << edata.id_;
  in_archive << edata.label_;
  in_archive << edata.attributes_;
  return in_archive;
}
inline OutArchive& operator>>(OutArchive& out_archive, EdgeData& edata) {
  out_archive >> edata.id_;
  out_archive >> edata.label_;
  out_archive >> edata.attributes_;
  return out_archive;
}
namespace internal {
template <>
inline const char* match(char const* str, grape::EdgeData& r, char const*) {
  char* match_end;
  // std::cout << "str=" << str << std::endl;
  r.label_ = std::strtol(str, &match_end, 10);
  // std::cout << "match = " << match_end << std::endl;
  r.id_ = std::strtol(match_end, &match_end, 10);

  std::string tmp_str(match_end);
  tmp_str = tmp_str.substr(0, tmp_str.length() - 1);  // remove the last '\n'
  std::vector<std::string> attrs = split(tmp_str, " ");
  r.attributes_.clear();
  for (auto attr : attrs) {
    if (attr == "" || attr == "\n") continue;
    r.attributes_.emplace_back(attr);
  }
  return match_end;
}
}  // namespace internal

struct EdgeTimeData {
  using LabelType = int32_t;
  using AttributeKeyType = size_t;
  using AttributeValueType = std::string;
  using AttributeContainer = std::vector<std::string>;
  EdgeTimeData() {}
  EdgeTimeData(int32_t id, int32_t label) : id_(id), label_(label) {}
  EdgeTimeData(int32_t id, int32_t label, int32_t time1)
      : id_(id), label_(label), time_(time1) {}
  int id_;
  int32_t label_;
  int32_t time_;
  bool operator==(const EdgeTimeData& edata) const {
    if (this->id_ != edata.id_ || this->label_ != edata.label_ ||
        this->time_ != edata.time_)
      return false;
    return true;
  }

  bool operator<(const EdgeTimeData& edata) const {
    return (this->label_ < edata.label_);
  }
};
inline InArchive& operator<<(InArchive& in_archive, const EdgeTimeData& edata) {
  in_archive << edata.id_;
  in_archive << edata.label_;
  in_archive << edata.time_;
  return in_archive;
}
inline OutArchive& operator>>(OutArchive& out_archive, EdgeTimeData& edata) {
  out_archive >> edata.id_;
  out_archive >> edata.label_;
  out_archive >> edata.time_;
  return out_archive;
}
namespace internal {
template <>
inline const char* match(char const* str, grape::EdgeTimeData& r, char const*) {
  char* match_end;
  // std::cout << "str=" << str << std::endl;
  r.label_ = std::strtol(str, &match_end, 10);
  // std::cout << "match = " << match_end << std::endl;
  r.id_ = std::strtol(match_end, &match_end, 10);

  r.time_ = std::strtol(match_end, &match_end, 10);
  return nullptr;
}
}  // namespace internal
}  // namespace grape
#endif