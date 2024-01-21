#ifndef EXAMPLES_ANALYTICAL_LITERAL_INFO_H_
#define EXAMPLES_ANALYTICAL_LITERAL_INFO_H_

#include "gundam/data_type/datatype.h"
#include "gundam/component/generator.h"
#include "gundam/match/match.h"

#include "gundam/type_getter/edge_label.h"
#include "gundam/type_getter/edge_id.h"
#include "gundam/type_getter/vertex_label.h"
#include "gundam/type_getter/vertex_id.h"

#include "gar/literal.h"

#include <cstdint>

#include <map>
#include <set>

namespace gar {

template<typename   PatternType,
         typename DataGraphType>
class LiteralInfo;

template<typename   PatternType,
         typename DataGraphType>
std::string& operator<<(std::string& out_string, 
     const LiteralInfo<PatternType, 
                     DataGraphType>& literal_info) {
            
  using LiteralInfoType
      = LiteralInfo<PatternType, 
                  DataGraphType>;

  uint8_t literal_type = static_cast<uint8_t>(literal_info.literal_type());

  out_string = std::move(out_string) + " <l " 
             + std::to_string(literal_type);
  
  switch (literal_info.literal_type()){
   case LiteralType::kAttrValueLiteral:
    // x.A
    out_string = std::move(out_string) 
               + " " + GUNDAM::ToString(literal_info.x_id())
               + " " + GUNDAM::ToString(literal_info.x_attr_key());
    break;
   case LiteralType::kVariableLiteral:
    // x.A == y.B
    out_string = std::move(out_string) 
               + " " + GUNDAM::ToString(literal_info.x_id())
               + " " + GUNDAM::ToString(literal_info.x_attr_key())
               + " " + GUNDAM::ToString(literal_info.y_id())
               + " " + GUNDAM::ToString(literal_info.y_attr_key());
    break;
   case LiteralType::kConstantLiteral:
    // x.A == c
    out_string = std::move(out_string) 
               + " " + GUNDAM::ToString(literal_info.x_id())
               + " " + GUNDAM::ToString(literal_info.x_attr_key())
               + " " + GUNDAM::ToString(literal_info.c_str())
               + " " + GUNDAM::ToString(literal_info.data_type());
    break;
   case LiteralType::kEdgeLiteral:
    // x_id_ - edge_label_ -> y_id_
    out_string = std::move(out_string) 
               + " " + GUNDAM::ToString(literal_info.x_id())
               + " " + GUNDAM::ToString(literal_info.y_id())
               + " " + GUNDAM::ToString(literal_info.edge_label());
    break;
   case LiteralType::kMlLiteral:
    // x_id_ - (module_url_: edge_label_) -> y_id_
    out_string = std::move(out_string) 
               + " " + GUNDAM::ToString(literal_info.x_id())
               + " " + GUNDAM::ToString(literal_info.y_id())
               + " " + GUNDAM::ToString(literal_info.edge_label())
               + " " + GUNDAM::ToString(literal_info.module_url());
    break;
   case LiteralType::kKeyLiteral:
    // x_id_ == y_id_
    out_string = std::move(out_string) 
               + " " + GUNDAM::ToString(literal_info.x_id())
               + " " + GUNDAM::ToString(literal_info.y_id());
    break;
   default:
    // unknown literal type
    assert(false);
    break;
  }
  out_string = std::move(out_string) + " >";
  return out_string;
}

template<typename   PatternType,
         typename DataGraphType>
std::string& operator>>(std::string& in_string, 
             LiteralInfo<PatternType, 
                       DataGraphType>& literal) {  
  
  using namespace GUNDAM;
            
  using LiteralInfoType 
      = LiteralInfo<PatternType, 
                  DataGraphType>;

  using VertexIDType = typename GUNDAM::VertexID<PatternType>::type;

  using AttributeKeyType = typename DataGraphType
                                     ::VertexType
                               ::AttributeKeyType;

  using EdgeLabelType = typename GUNDAM::EdgeLabel<PatternType>::type;

  std::stringstream ss;
  ss << in_string;

  std::string str;
  
  ss>>str;
  assert(str == "<l");

  ss>>str;
  switch(static_cast<LiteralType>(std::stoi(str))){
   case LiteralType::kAttrValueLiteral:
    // x.A
    {
      ss>>str;
      VertexIDType x_id = GUNDAM::StringToDataType<VertexIDType>(str);
      ss>>str;
      AttributeKeyType x_attr_key = GUNDAM::StringToDataType<AttributeKeyType>(str);
      literal = LiteralInfoType(x_id, x_attr_key);
    }
    break;
   case LiteralType::kVariableLiteral:
    // x.A == y.B
    {
      ss>>str;
      VertexIDType x_id = GUNDAM::StringToDataType<VertexIDType>(str);
      ss>>str;
      AttributeKeyType x_attr_key = GUNDAM::StringToDataType<AttributeKeyType>(str);
      ss>>str;
      VertexIDType y_id = GUNDAM::StringToDataType<VertexIDType>(str);
      ss>>str;
      AttributeKeyType y_attr_key = GUNDAM::StringToDataType<AttributeKeyType>(str);
      literal = LiteralInfoType(x_id, x_attr_key,
                                y_id, y_attr_key);
    }
    break;
   case LiteralType::kConstantLiteral:
    // x.A == c
    {
      ss>>str;
      VertexIDType x_id = GUNDAM::StringToDataType<VertexIDType>(str);
      ss>>str;
      AttributeKeyType x_attr_key = GUNDAM::StringToDataType<AttributeKeyType>(str);
      ss>>str;
      std::string value_str = str;
      ss>>str;
      enum GUNDAM::BasicDataType value_type = GUNDAM::StringToEnum(str.c_str());
      literal = LiteralInfoType(x_id, x_attr_key, value_str, value_type);
    }
    break;
   case LiteralType::kEdgeLiteral:
    // x_id_ - edge_label_ -> y_id_
    {
      ss>>str;
      VertexIDType x_id = GUNDAM::StringToDataType<VertexIDType>(str);
      ss>>str;
      VertexIDType y_id = GUNDAM::StringToDataType<VertexIDType>(str);
      ss>>str;
      EdgeLabelType edge_label = GUNDAM::StringToDataType<EdgeLabelType>(str);
      literal = LiteralInfoType(x_id, y_id, edge_label);
    }
    break;
   case LiteralType::kMlLiteral:
    // x_id_ - (module_url_: edge_label_) -> y_id_
    {
      ss>>str;
      VertexIDType x_id = GUNDAM::StringToDataType<VertexIDType>(str);
      ss>>str;
      VertexIDType y_id = GUNDAM::StringToDataType<VertexIDType>(str);
      ss>>str;
      EdgeLabelType edge_label = GUNDAM::StringToDataType<EdgeLabelType>(str);
      ss>>str; // module_url 
      literal = LiteralInfoType(x_id, y_id, edge_label, str);
    }
    break;
   case LiteralType::kKeyLiteral:
    // x_id_ - (module_url_: edge_label_) -> y_id_
    {
      ss>>str;
      VertexIDType x_id = GUNDAM::StringToDataType<VertexIDType>(str);
      ss>>str;
      VertexIDType y_id = GUNDAM::StringToDataType<VertexIDType>(str);
      literal = LiteralInfoType(x_id, y_id);
    }
    break;
   default:
    // unknown literal type
    assert(false);
    break;
  }
  ss>>str;
  assert(str == ">");
  getline(ss, in_string);
  return in_string;
}

// assume that all value has the same type of value
template<typename   PatternType, 
         typename DataGraphType>
class LiteralInfo{
 private:
  enum LiteralType literal_type_;

  template <typename   PatternType2,
            typename DataGraphType2>
  friend class LiteralInfo;

  using VertexIDType = typename GUNDAM::VertexID<PatternType>::type;

  using AttributeKeyType = typename DataGraphType
                                     ::VertexType
                               ::AttributeKeyType;

  using EdgeLabelType = typename GUNDAM::EdgeLabel<DataGraphType>::type;

 public:
  LiteralInfo():literal_type_(LiteralType::kNoneLiteral){
    return;
  }

  template <typename   PatternType2,
            typename DataGraphType2>
  LiteralInfo(const LiteralInfo<PatternType2,
                              DataGraphType2>& literal) 
          :literal_type_(literal.literal_type_),
                   x_id_(literal.x_id_),
                   y_id_(literal.y_id_),
             x_attr_key_(literal.x_attr_key_),
             y_attr_key_(literal.y_attr_key_),
                  c_str_(literal.c_str_),
             edge_label_(literal.edge_label_),  
              data_type_(literal.data_type_){
    return;
  }

  // the construction method for variable literal
  LiteralInfo(const     VertexIDType& x_id,
              const AttributeKeyType& x_attr_key)
                       :literal_type_(LiteralType::kAttrValueLiteral),
                                x_id_(x_id),
                          x_attr_key_(x_attr_key){
    return;
  }

  // the construction method for variable literal
  LiteralInfo(const     VertexIDType& x_id,
              const AttributeKeyType& x_attr_key,
              const      std::string& c_str,
              const enum GUNDAM::BasicDataType& data_type)
                       :literal_type_(LiteralType::kConstantLiteral),
                                x_id_(x_id),
                          x_attr_key_(x_attr_key),
                               c_str_(c_str),
                           data_type_(data_type){
    return;
  }

  // the construction method for variable literal
  LiteralInfo(const     VertexIDType& x_id,
              const AttributeKeyType& x_attr_key,
              const     VertexIDType& y_id,
              const AttributeKeyType& y_attr_key)
                       :literal_type_(LiteralType::kVariableLiteral),
                                x_id_(x_id),
                          x_attr_key_(x_attr_key),
                                y_id_(y_id),
                          y_attr_key_(y_attr_key){
    return;
  }
  
  // the construction method for edge literal
  LiteralInfo(const  VertexIDType& x_id,
              const  VertexIDType& y_id,
              const EdgeLabelType& edge_label)
                    :literal_type_(LiteralType::kEdgeLiteral),
                             x_id_(x_id),
                             y_id_(y_id),
                       edge_label_(edge_label){
    return;
  }
  
  // the construction method for ml literal
  LiteralInfo(const  VertexIDType& x_id,
              const  VertexIDType& y_id,
              const EdgeLabelType& edge_label,
              const   std::string& module_url)
                    :literal_type_(LiteralType::kMlLiteral),
                             x_id_(x_id),
                             y_id_(y_id),
                       edge_label_(edge_label),
                       module_url_(module_url){
    return;
  }
  
  // the construction method for key literal
  LiteralInfo(const  VertexIDType& x_id,
              const  VertexIDType& y_id)
                    :literal_type_(LiteralType::kKeyLiteral),
                             x_id_(x_id),
                             y_id_(y_id){
    return;
  }

  LiteralInfo(const LiteralInfo&) = default;
  LiteralInfo(LiteralInfo&&) = default;

  LiteralInfo& operator=(const LiteralInfo&) = default;	
  LiteralInfo& operator=(LiteralInfo&&) = default;	

  inline bool operator< (const LiteralInfo& literal) const {
    // ##############
    if (this->literal_type_ != literal.literal_type_) {
      return this->literal_type_ < literal.literal_type_;
    }
    // ##############
    assert( this->literal_type_ == literal.literal_type_);
    switch (this->literal_type_){
     case LiteralType::kAttrValueLiteral:
      // x.A
      if (this->x_id_ != literal.x_id_) {
        return this->x_id_ < literal.x_id_;
      }
      assert(this->x_id_ == literal.x_id_);
      return this->x_attr_key_  < literal.x_attr_key_;
     case LiteralType::kVariableLiteral: {
        auto  this_id_0 = this->x_id_ < this->y_id_?
                          this->x_id_ : this->y_id_;
        auto  this_id_1 = this->x_id_ < this->y_id_?
                          this->y_id_ : this->x_id_;
        auto this_attr_key_0 = this->x_id_ < this->y_id_?
                               this->x_attr_key_ 
                             : this->y_attr_key_;
        auto this_attr_key_1 = this->x_id_ < this->y_id_?
                               this->y_attr_key_ 
                             : this->x_attr_key_;

        auto literal_id_0 = literal.x_id_ < literal.y_id_?
                            literal.x_id_ : literal.y_id_;
        auto literal_id_1 = literal.x_id_ < literal.y_id_?
                            literal.y_id_ : literal.x_id_;
        auto literal_attr_key_0 = literal.x_id_ < literal.y_id_?
                                  literal.x_attr_key_ 
                                : literal.y_attr_key_;
        auto literal_attr_key_1 = literal.x_id_ < literal.y_id_?
                                  literal.y_attr_key_ 
                                : literal.x_attr_key_;

        if (this_id_0 != literal_id_0) {
          return this_id_0 < literal_id_0;
        }
        if (this_id_1 != literal_id_1) {
          return this_id_1 < literal_id_1;
        }
        if (this_attr_key_0 != literal_attr_key_0) {
          return this_attr_key_0 < literal_attr_key_0;
        }
        return this_attr_key_1 < literal_attr_key_1;
      }
     case LiteralType::kConstantLiteral:
      // x.A == c
      if (this->x_id_ != literal.x_id_) {
        return this->x_id_ < literal.x_id_;
      }
      assert(this->x_id_ == literal.x_id_);
      
      if (this->x_attr_key_ != literal.x_attr_key_) {
        return this->x_attr_key_ < literal.x_attr_key_;
      }
      assert(this->x_attr_key_ == literal.x_attr_key_);
      
      return this->c_str_ < literal.c_str_;
     case LiteralType::kEdgeLiteral:
      if (this->x_id_ < literal.x_id_) {
        return true;
      }
      if (this->x_id_ > literal.x_id_) {
        return false;
      }
      assert(this->x_id_ == literal.x_id_);
      
      if (this->y_id_ < literal.y_id_) {
        return true;
      }
      if (this->y_id_ > literal.y_id_) {
        return false;
      }
      assert(this->y_id_ == literal.y_id_);
      
      return this->edge_label_ < literal.edge_label_;
     case LiteralType::kMlLiteral:
      // x_id_ - (module_url_: edge_label_) -> y_id_
      if (this->x_id_ < literal.x_id_) {
        return true;
      }
      if (this->x_id_ > literal.x_id_) {
        return false;
      }
      assert(this->x_id_ == literal.x_id_);
      
      if (this->y_id_ < literal.y_id_) {
        return true;
      }
      if (this->y_id_ > literal.y_id_) {
        return false;
      }
      assert(this->y_id_ == literal.y_id_);
      
      if (this->edge_label_ < literal.edge_label_) {
        return true;
      }
      if (this->edge_label_ > literal.edge_label_) {
        return false;
      }
      assert(this->edge_label_ == literal.edge_label_);
      
      return this->module_url_ < literal.module_url_;
     case LiteralType::kKeyLiteral:
      // x_id_ == y_id_
      if (this->x_id_ < literal.x_id_) {
        return true;
      }
      if (this->x_id_ > literal.x_id_) {
        return false;
      }
      assert(this->x_id_ == literal.x_id_);
      return this->y_id_ < literal.y_id_;
     default:
      return false;
    }
    return false;
  }

  inline bool operator==(const LiteralInfo& literal) const {
    if (this->literal_type_ != literal.literal_type_){
      // does not have same type
      return false;
    }
    switch (this->literal_type_){
     case LiteralType::kAttrValueLiteral:
      // x.A
      if (this->x_id_       == literal.x_id_
       && this->x_attr_key_ == literal.x_attr_key_){
        return true;
      }
      return false;
     case LiteralType::kVariableLiteral:
      // x.A == y.B

      //  x.A == y.B is equal to y.B == x.A
      if((this->x_id_       == literal.x_id_
       && this->y_id_       == literal.y_id_
       && this->x_attr_key_ == literal.x_attr_key_
       && this->y_attr_key_ == literal.y_attr_key_)
       ||(this->x_id_       == literal.y_id_
       && this->y_id_       == literal.x_id_
       && this->x_attr_key_ == literal.y_attr_key_
       && this->y_attr_key_ == literal.x_attr_key_)){
        return true;
      }
      return false;
     case LiteralType::kConstantLiteral:
      // x.A == c
      if (this->x_id_       == literal.x_id_
       && this->x_attr_key_ == literal.x_attr_key_
       && this->c_str_      == literal.c_str_){
        return true;
      }
      return false;
     case LiteralType::kEdgeLiteral:
      // x_id_ - edge_label_ -> y_id_
      if (this->x_id_       == literal.x_id_
       && this->y_id_       == literal.y_id_
       && this->edge_label_ == literal.edge_label_){
        return true;
      }
      return false;
     case LiteralType::kMlLiteral:
      // x_id_ - (module_url_: edge_label_) -> y_id_
      if (this->x_id_       == literal.x_id_
       && this->y_id_       == literal.y_id_
       && this->module_url_ == literal.module_url_
       && this->edge_label_ == literal.edge_label_){
        return true;
      }
      return false;
     case LiteralType::kKeyLiteral:
      // x_id_ == y_id_

      //  x_id_ == y_id_ is equal to y_id_ == x_id_
      if ((this->x_id_ == literal.x_id_
        && this->y_id_ == literal.y_id_)
       || (this->x_id_ == literal.y_id_
        && this->y_id_ == literal.x_id_)){
        return true;
      }
      return false;
     default:
      // unknown literal type
      assert(false);
      break;
    }
    // should not reach 
    assert(false);
    return false;
  }

  bool Satisfy(const GUNDAM::Match<PatternType, 
                                 DataGraphType>& match) const {

    assert(this->MapBy(match));

    switch (this->literal_type()){
    case LiteralType::kAttrValueLiteral: {
      // x.A 
      return AttributeLiteral<PatternType, DataGraphType>::Satisfy(match.MapTo(this->x_id_), 
                                                                               this->x_attr_key_);
    }
    case LiteralType::kVariableLiteral: {
      // x.A == y.B
      return VariableLiteral<PatternType, DataGraphType>::Satisfy(match.MapTo(this->x_id_), 
                                                                  match.MapTo(this->y_id_), 
                                                                              this->x_attr_key_, 
                                                                              this->y_attr_key_);
    }
    case LiteralType::kConstantLiteral: {
      // x.A == c
      return ConstantLiteral<PatternType, DataGraphType>::Satisfy(match.MapTo(this->x_id_),
                                                                              this->x_attr_key_,
                                                                              this->c_str_,
                                                                              this->data_type_);
    }
    case LiteralType::kEdgeLiteral: {
      // x_id_ - edge_label_ -> y_id_
      return EdgeLiteral<PatternType, DataGraphType>::Satisfy(match.MapTo(this->x_id_),
                                                              match.MapTo(this->y_id_), this->edge_label_);
    }
    case LiteralType::kMlLiteral: {
      // x_id_ - (module_url_: edge_label_) -> y_id_
      assert(false);
      return EdgeLiteral<PatternType, DataGraphType>::Satisfy(match.MapTo(this->x_id_),
                                                              match.MapTo(this->y_id_), this->edge_label_);
    }
    case LiteralType::kKeyLiteral: {
      // x_id_ == y_id_
      return KeyLiteral<PatternType, DataGraphType>::Satisfy(match.MapTo(this->x_id_),
                                                             match.MapTo(this->y_id_));
    }
    default:
      // unknown literal type
      assert(false);
      break;
    }
    return false;
  }

  bool MapBy(const std::map<typename GUNDAM::VertexHandle<  PatternType>::type,
                            typename GUNDAM::VertexHandle<DataGraphType>::type>& match) const {
    bool found_x = false,
         found_y = false;
    for (const auto& map : match) {
      if (map.first->id() == this->x_id()){
        found_x = true;
        if (found_y) {
          // both ids are found
          break;
        }
        if (this->literal_type() == LiteralType::  kConstantLiteral
         || this->literal_type() == LiteralType::kAttrValueLiteral) {
          // y id is not found but is not required
          break;
        }
      }
      if (map.first->id() == this->y_id()){
        found_y = true;
        if (found_x) {
          break;
        }
      }
    }
    switch (this->literal_type()){
     case LiteralType::kAttrValueLiteral:
      // x.A 
      return found_x;
     case LiteralType::kVariableLiteral:
      // x.A == y.B
      return found_x && found_y;
     case LiteralType::kConstantLiteral:
      // x.A == c
      return found_x;
     case LiteralType::kEdgeLiteral:
      // x_id_ - edge_label_ -> y_id_
      return found_x && found_y;
     case LiteralType::kMlLiteral:
      // x_id_ - (module_url_: edge_label_) -> y_id_
      return found_x && found_y;
     case LiteralType::kKeyLiteral:
      // x_id_ == y_id_
      return found_x && found_y;
     default:
      // unknown literal type
      assert(false);
      break;
    }
    return false;
  }

  // optimize me
  template<typename NewPatternType>
  LiteralInfo<NewPatternType,
               DataGraphType> MapTo(const std::map<
                                          typename GUNDAM::VertexHandle<   PatternType>::type,
                                          typename GUNDAM::VertexHandle<NewPatternType>::type>& match) const {
    using NewLiteralInfoType = LiteralInfo<NewPatternType,
                                            DataGraphType>;

    using NewPatternIDType = typename GUNDAM::VertexID<NewPatternType>::type;
    
    NewPatternIDType dst_x_id = 0,
                     dst_y_id = 0;

    bool found_x = false,
         found_y = false;

    for (const auto& [src_handle,
                      dst_handle] : match) {
      if (src_handle->id() == this->x_id()){
        found_x = true;
        dst_x_id = dst_handle->id();
        if (found_y) {
          // both ids are found
          break;
        }
        if (this->literal_type() == LiteralType:: kConstantLiteral
         || this->literal_type() == LiteralType::kAttrValueLiteral) {
          // y id is not found but is not required
          break;
        }
      }
      if (src_handle->id() == this->y_id()){
        found_y = true;
        dst_y_id = dst_handle->id();
        if (found_x) {
          break;
        }
      }
    }

    switch (this->literal_type_){
     case LiteralType::kAttrValueLiteral:
      // x.A
      if (!found_x) {
        return NewLiteralInfoType();
      }
      return NewLiteralInfoType(dst_x_id, this->x_attr_key());
     case LiteralType::kVariableLiteral:
      // x.A == y.B
      if (!found_x || !found_y) {
        return NewLiteralInfoType();
      }
      return NewLiteralInfoType(dst_x_id, this->x_attr_key(),
                                dst_y_id, this->y_attr_key());
     case LiteralType::kConstantLiteral:
      // x.A == c
      if (!found_x) {
        return NewLiteralInfoType();
      }
      return NewLiteralInfoType(dst_x_id, this->x_attr_key(),
                                          this->c_str(),
                                          this->data_type());
     case LiteralType::kEdgeLiteral:
      if (!found_x || !found_y) {
        return NewLiteralInfoType();
      }
      // x_id_ - edge_label_ -> y_id_
      return NewLiteralInfoType(dst_x_id, dst_y_id, this->edge_label());
     case LiteralType::kMlLiteral:
      if (!found_x || !found_y) {
        return NewLiteralInfoType();
      }
      // x_id_ - (module_url_: edge_label_) -> y_id_
      return NewLiteralInfoType(dst_x_id, dst_y_id, this->edge_label(),
                                                    this->module_url());
     case LiteralType::kKeyLiteral:
      // x_id_ == y_id_
      return NewLiteralInfoType(dst_x_id, dst_y_id);
     default:
      // unknown literal type
      assert(false);
      break;
    }
    return NewLiteralInfoType();
  }

  // optimize me
  template<typename NewPatternType,
           typename ThisPatternType,
           std::enable_if_t<std::is_same_v<      PatternType, ThisPatternType>
                         || std::is_same_v<const PatternType, ThisPatternType>,
                            bool> = false>
  LiteralInfo<NewPatternType,
               DataGraphType> MapTo(const GUNDAM::Match<ThisPatternType,
                                                         NewPatternType>& match) const {
    using NewLiteralInfoType = LiteralInfo<NewPatternType,
                                            DataGraphType>;

    using NewPatternIDType = typename GUNDAM::VertexID<NewPatternType>::type;
    
    NewPatternIDType dst_x_id = 0,
                     dst_y_id = 0;

    bool found_x = false,
         found_y = false;

    for (auto map_it = match.MapBegin();
             !map_it.IsDone();
              map_it++) {
      auto src_handle = map_it->src_handle();
      auto dst_handle = map_it->dst_handle();
      if (src_handle->id() == this->x_id()){
        found_x = true;
        dst_x_id = dst_handle->id();
        if (found_y) {
          // both ids are found
          break;
        }
        if (this->literal_type() == LiteralType::  kConstantLiteral
         || this->literal_type() == LiteralType::kAttrValueLiteral) {
          // y id is not found but is not required
          break;
        }
      }
      if (src_handle->id() == this->y_id()){
        found_y = true;
        dst_y_id = dst_handle->id();
        if (found_x) {
          break;
        }
      }
    }

    switch (this->literal_type_){
     case LiteralType::kAttrValueLiteral:
      // x.A
      if (!found_x) {
        return NewLiteralInfoType();
      }
      return NewLiteralInfoType(dst_x_id, this->x_attr_key());
     case LiteralType::kVariableLiteral:
      // x.A == y.B
      if (!found_x || !found_y) {
        return NewLiteralInfoType();
      }
      return NewLiteralInfoType(dst_x_id, this->x_attr_key(),
                                dst_y_id, this->y_attr_key());
     case LiteralType::kConstantLiteral:
      // x.A == c
      if (!found_x) {
        return NewLiteralInfoType();
      }
      return NewLiteralInfoType(dst_x_id, this->x_attr_key(),
                                          this->c_str(),
                                          this->data_type());
     case LiteralType::kEdgeLiteral:
      if (!found_x || !found_y) {
        return NewLiteralInfoType();
      }
      // x_id_ - edge_label_ -> y_id_
      return NewLiteralInfoType(dst_x_id, dst_y_id, this->edge_label());
     case LiteralType::kMlLiteral:
      if (!found_x || !found_y) {
        return NewLiteralInfoType();
      }
      // x_id_ - (module_url_: edge_label_) -> y_id_
      return NewLiteralInfoType(dst_x_id, dst_y_id, this->edge_label(),
                                                    this->module_url());
     case LiteralType::kKeyLiteral:
      // x_id_ == y_id_
      return NewLiteralInfoType(dst_x_id, dst_y_id);
     default:
      // unknown literal type
      assert(false);
      break;
    }
    return NewLiteralInfoType();
  }
  
  inline std::vector<VertexIDType> vertex_id_set() const {
    std::vector<VertexIDType> vertex_id_set;
    switch(this->literal_type_){
    case LiteralType::kAttrValueLiteral:
      // x.A
      vertex_id_set.emplace_back(this->x_id_);
      break;
    case LiteralType::kVariableLiteral:
      // x.A == y.B
      vertex_id_set.emplace_back(this->x_id_);
      vertex_id_set.emplace_back(this->y_id_);
      break;
    case LiteralType::kConstantLiteral:
      // x.A == c
      vertex_id_set.emplace_back(this->x_id_);
      break;
    case LiteralType::kEdgeLiteral:
      // x_id_ - edge_label_ -> y_id_
      vertex_id_set.emplace_back(this->x_id_);
      vertex_id_set.emplace_back(this->y_id_);
      break;
    case LiteralType::kMlLiteral:
      // x_id_ - (module_url_: edge_label_) -> y_id_
      vertex_id_set.emplace_back(this->x_id_);
      vertex_id_set.emplace_back(this->y_id_);
      break;
     case LiteralType::kKeyLiteral:
      // x_id_ == y_id_
      vertex_id_set.emplace_back(this->x_id_);
      vertex_id_set.emplace_back(this->y_id_);
      break;
    default:
      // unknown literal type
      assert(false);
      break;
    }
    return vertex_id_set;
  }

  inline bool operator!=(const LiteralInfo& literal) const {
    return !((*this) == literal);
  }

  inline const enum LiteralType literal_type() const{
    return this->literal_type_;
  }

  inline const VertexIDType& x_id() const{
    return this->x_id_;
  }

  inline const VertexIDType& y_id() const{
    return this->y_id_;
  }

  inline const AttributeKeyType& x_attr_key() const{
    return this->x_attr_key_;
  }

  inline const AttributeKeyType& y_attr_key() const{
    return this->y_attr_key_;
  }

  inline const std::string& c_str() const{
    return this->c_str_;
  }

  inline const EdgeLabelType& edge_label() const{
    return this->edge_label_;
  }

  inline enum GUNDAM::BasicDataType data_type() const{
    return this->data_type_;
  }

  inline const std::string& module_url() const{
    return this->module_url_;
  }

 private:
  VertexIDType x_id_, y_id_;

  AttributeKeyType x_attr_key_, 
                   y_attr_key_;

  std::string c_str_;

  EdgeLabelType edge_label_;

  enum GUNDAM::BasicDataType data_type_;

  std::string module_url_;
};

}  // namespace gar

#endif