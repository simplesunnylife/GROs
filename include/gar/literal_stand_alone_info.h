#ifndef EXAMPLES_ANALYTICAL_LITERAL_STAND_ALONE_INFO_H_
#define EXAMPLES_ANALYTICAL_LITERAL_STAND_ALONE_INFO_H_

#include "gundam/data_type/datatype.h"
#include "gundam/component/generator.h"

#include <cstdint>

#include <map>
#include <set>

namespace gar {
  
enum class LiteralType: char{
       kNoneLiteral = 'n',
   kVariableLiteral = 'v',
  kAttrValueLiteral = 'a',
   kConstantLiteral = 'c',
       kEdgeLiteral = 'e',
         kMlLiteral = 'm',
        kKeyLiteral = 'k'
};

template<typename   PatternType,
         typename DataGraphType>
class LiteralStandAloneInfo;

template<typename   PatternType,
         typename DataGraphType>
std::string& operator<<(std::string& out_string, 
     const LiteralStandAloneInfo<PatternType, 
                               DataGraphType>& literal_class) {

  uint8_t literal_type = static_cast<uint8_t>(literal_class.literal_type());

  out_string = std::move(out_string) + " <l " 
             + std::to_string(literal_type);
  
  switch (literal_class.literal_type()){
   case LiteralType::kAttrValueLiteral:
    // x.A
    out_string = std::move(out_string) 
               + " " + GUNDAM::ToString(literal_class.x_label())
               + " " + GUNDAM::ToString(literal_class.x_attr_key());
    break;
   case LiteralType::kVariableLiteral:
    // x.A == y.B
    out_string = std::move(out_string) 
               + " " + GUNDAM::ToString(literal_class.x_label())
               + " " + GUNDAM::ToString(literal_class.x_attr_key())
               + " " + GUNDAM::ToString(literal_class.y_label())
               + " " + GUNDAM::ToString(literal_class.y_attr_key());
    break;
   case LiteralType::kConstantLiteral:
    // x.A == c
    out_string = std::move(out_string) 
               + " " + GUNDAM::ToString(literal_class.x_label())
               + " " + GUNDAM::ToString(literal_class.x_attr_key())
               + " " + GUNDAM::ToString(literal_class.c_str());
    break;
   case LiteralType::kEdgeLiteral:
    // x_id_ - edge_label_ -> y_id_
    out_string = std::move(out_string) 
               + " " + GUNDAM::ToString(literal_class.x_label())
               + " " + GUNDAM::ToString(literal_class.y_label())
               + " " + GUNDAM::ToString(literal_class.edge_label());
    break;
   case LiteralType::kMlLiteral:
    // x_id_ - (module_url_: edge_label_) -> y_id_
    out_string = std::move(out_string) 
               + " " + GUNDAM::ToString(literal_class.x_label())
               + " " + GUNDAM::ToString(literal_class.y_label())
               + " " + GUNDAM::ToString(literal_class.edge_label())
               + " " + GUNDAM::ToString(literal_class.module_url());
    break;
   case LiteralType::kKeyLiteral:
    // x_id_ - (module_url_: edge_label_) -> y_id_
    out_string = std::move(out_string) 
               + " " + GUNDAM::ToString(literal_class.x_label())
               + " " + GUNDAM::ToString(literal_class.y_label());
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
             LiteralStandAloneInfo<PatternType, 
                                 DataGraphType>& literal) {  
  
  using namespace GUNDAM;
            
  using LiteralStandAloneInfoType 
      = LiteralStandAloneInfo<PatternType, 
                            DataGraphType>;

  using VertexLabelType = typename PatternType
                                  ::VertexType
                                   ::LabelType;

  using AttributeKeyType = typename DataGraphType
                                     ::VertexType
                               ::AttributeKeyType;

  using EdgeLabelType = typename DataGraphType
                                    ::EdgeType
                                   ::LabelType;

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
      VertexLabelType x_label = GUNDAM::StringToDataType<VertexLabelType>(str);
      ss>>str;
      AttributeKeyType x_attr_key = GUNDAM::StringToDataType<AttributeKeyType>(str);
      literal = LiteralStandAloneInfoType(x_label, x_attr_key);
    }
    break;
   case LiteralType::kVariableLiteral:
    // x.A == y.B
    {
      ss>>str;
      VertexLabelType x_label = GUNDAM::StringToDataType<VertexLabelType>(str);
      ss>>str;
      AttributeKeyType x_attr_key = GUNDAM::StringToDataType<AttributeKeyType>(str);
      ss>>str;
      VertexLabelType y_label = GUNDAM::StringToDataType<VertexLabelType>(str);
      ss>>str;
      AttributeKeyType y_attr_key = GUNDAM::StringToDataType<AttributeKeyType>(str);
      literal = LiteralStandAloneInfoType(x_label, x_attr_key,
                                y_label, y_attr_key);
    }
    break;
   case LiteralType::kConstantLiteral:
    // x.A == c
    {
      ss>>str;
      VertexLabelType x_label = GUNDAM::StringToDataType<VertexLabelType>(str);
      ss>>str;
      AttributeKeyType x_attr_key = GUNDAM::StringToDataType<AttributeKeyType>(str);
      ss>>str;
      std::string value_str = str;
      ss>>str;
      enum GUNDAM::BasicDataType value_type = GUNDAM::StringToEnum(str.c_str());
      literal = LiteralStandAloneInfoType(x_label, x_attr_key, value_str, value_type);
    }
    break;
   case LiteralType::kEdgeLiteral:
    // x_id_ - edge_label_ -> y_id_
    {
      ss>>str;
      VertexLabelType x_label = GUNDAM::StringToDataType<VertexLabelType>(str);
      ss>>str;
      VertexLabelType y_label = GUNDAM::StringToDataType<VertexLabelType>(str);
      ss>>str;
      EdgeLabelType edge_label = GUNDAM::StringToDataType<EdgeLabelType>(str);
      literal = LiteralStandAloneInfoType(x_label, y_label, edge_label);
    }
    break;
   case LiteralType::kMlLiteral:
    // x_id_ - (module_url_: edge_label_) -> y_id_
    {
      ss>>str;
      VertexLabelType x_label = GUNDAM::StringToDataType<VertexLabelType>(str);
      ss>>str;
      VertexLabelType y_label = GUNDAM::StringToDataType<VertexLabelType>(str);
      ss>>str;
      EdgeLabelType edge_label = GUNDAM::StringToDataType<EdgeLabelType>(str);
      ss>>str; // module_url 
      literal = LiteralStandAloneInfoType(x_label, y_label, edge_label, str);
    }
    break;
   case LiteralType::kKeyLiteral:
    // x_id_ == y_id_
    {
      ss>>str;
      VertexLabelType x_label = GUNDAM::StringToDataType<VertexLabelType>(str);
      ss>>str;
      VertexLabelType y_label = GUNDAM::StringToDataType<VertexLabelType>(str);
      literal = LiteralStandAloneInfoType(x_label, y_label);
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
class LiteralStandAloneInfo{
 private:
  enum LiteralType literal_type_;

  template <typename   PatternType2,
            typename DataGraphType2>
  friend class LiteralStandAloneInfo;

  using VertexLabelType = typename PatternType
                                  ::VertexType
                                      ::IDType;

  using AttributeKeyType = typename DataGraphType
                                     ::VertexType
                               ::AttributeKeyType;

  using EdgeLabelType = typename DataGraphType
                                    ::EdgeType
                                   ::LabelType;

 public:
  LiteralStandAloneInfo():literal_type_(LiteralType::kNoneLiteral){
    return;
  }

  template <typename   PatternType2,
            typename DataGraphType2>
  LiteralStandAloneInfo(const LiteralStandAloneInfo<PatternType2,
                                                  DataGraphType2>& literal) 
          :literal_type_(literal.literal_type_),
                x_label_(literal.x_label_),
                y_label_(literal.y_label_),
             x_attr_key_(literal.x_attr_key_),
             y_attr_key_(literal.y_attr_key_),
                  c_str_(literal.c_str_),
             edge_label_(literal.edge_label_),  
              data_type_(literal.data_type_){
    return;
  }

  // the construction method for variable literal
  LiteralStandAloneInfo(const  VertexLabelType& x_label,
                        const AttributeKeyType& x_attr_key)
                       :literal_type_(LiteralType::kAttrValueLiteral),
                             x_label_(x_label),
                          x_attr_key_(x_attr_key){
    return;
  }

  // the construction method for variable literal
  LiteralStandAloneInfo(const  VertexLabelType& x_label,
                        const AttributeKeyType& x_attr_key,
                        const      std::string& c_str,
                        const enum GUNDAM::BasicDataType& data_type)
                       :literal_type_(LiteralType::kConstantLiteral),
                             x_label_(x_label),
                          x_attr_key_(x_attr_key),
                               c_str_(c_str),
                           data_type_(data_type){
    return;
  }

  // the construction method for variable literal
  LiteralStandAloneInfo(const  VertexLabelType& x_label,
                        const AttributeKeyType& x_attr_key,
                        const  VertexLabelType& y_label,
                        const AttributeKeyType& y_attr_key)
                       :literal_type_(LiteralType::kVariableLiteral),
                             x_label_(x_label),
                          x_attr_key_(x_attr_key),
                             y_label_(y_label),
                          y_attr_key_(y_attr_key){
    return;
  }
  
  // the construction method for edge literal
  LiteralStandAloneInfo(const VertexLabelType&    x_label,
                        const VertexLabelType&    y_label,
                        const   EdgeLabelType& edge_label)
                    :literal_type_(LiteralType::kEdgeLiteral),
                          x_label_(x_label),
                          y_label_(y_label),
                       edge_label_(edge_label){
    return;
  }
  
  // the construction method for ml literal
  LiteralStandAloneInfo(const VertexLabelType&    x_label,
                        const VertexLabelType&    y_label,
                        const   EdgeLabelType& edge_label,
                        const     std::string& module_url)
                    :literal_type_(LiteralType::kMlLiteral),
                          x_label_(x_label),
                          y_label_(y_label),
                       edge_label_(edge_label),
                           module_url_(module_url){
    return;
  }
  
  // the construction method for ml literal
  LiteralStandAloneInfo(const VertexLabelType&    x_label,
                        const VertexLabelType&    y_label)
                    :literal_type_(LiteralType::kKeyLiteral),
                          x_label_(x_label),
                          y_label_(y_label){
    return;
  }

  inline bool operator< (const LiteralStandAloneInfo& literal) const {
    // ##############
    if (this->literal_type_ < literal.literal_type_) {
      return true;
    }
    else if (this->literal_type_ > literal.literal_type_) {
      return false;
    }
    // ##############
    assert( this->literal_type_ 
       == literal.literal_type_);
    switch (this->literal_type_){
     case LiteralType::kAttrValueLiteral:
      // x.A
      if (this->x_label_ < literal.x_label_) {
        return true;
      }
      if (this->x_label_ > literal.x_label_) {
        return false;
      }
      assert(this->x_label_ == literal.x_label_);
      if (this->x_attr_key_  < literal.x_attr_key_){
        return true;
      }
      return false;
     case LiteralType::kVariableLiteral:
      // x.A == y.B
      if (this->x_label_ < literal.x_label_) {
        return true;
      }
      if (this->x_label_ > literal.x_label_) {
        return false;
      }
      assert(this->x_label_ == literal.x_label_);
      
      if (this->y_label_ < literal.y_label_) {
        return true;
      }
      if (this->y_label_ > literal.y_label_) {
        return false;
      }
      assert(this->y_label_ == literal.y_label_);
      
      if (this->x_attr_key_ < literal.x_attr_key_) {
        return true;
      }
      if (this->x_attr_key_ > literal.x_attr_key_) {
        return false;
      }
      assert(this->x_attr_key_ == literal.x_attr_key_);
      
      if (this->y_attr_key_ < literal.y_attr_key_) {
        return true;
      }
      return false;
     case LiteralType::kConstantLiteral:
      // x.A == c
      if (this->x_label_ < literal.x_label_) {
        return true;
      }
      if (this->x_label_ > literal.x_label_) {
        return false;
      }
      assert(this->x_label_ == literal.x_label_);
      
      if (this->x_attr_key_ < literal.x_attr_key_) {
        return true;
      }
      if (this->x_attr_key_ > literal.x_attr_key_) {
        return false;
      }
      assert(this->x_attr_key_ == literal.x_attr_key_);
      
      if (this->c_str_ < literal.c_str_) {
        return true;
      }
      return false;
     case LiteralType::kEdgeLiteral:
      if (this->x_label_ < literal.x_label_) {
        return true;
      }
      if (this->x_label_ > literal.x_label_) {
        return false;
      }
      assert(this->x_label_ == literal.x_label_);
      
      if (this->y_label_ < literal.y_label_) {
        return true;
      }
      if (this->y_label_ > literal.y_label_) {
        return false;
      }
      assert(this->y_label_ == literal.y_label_);
      
      if (this->edge_label_ < literal.edge_label_) {
        return true;
      }
      return false;
     case LiteralType::kMlLiteral:
      // x_id_ - (module_url_: edge_label_) -> y_id_
      if (this->x_label_ < literal.x_label_) {
        return true;
      }
      if (this->x_label_ > literal.x_label_) {
        return false;
      }
      assert(this->x_label_ == literal.x_label_);
      
      if (this->y_label_ < literal.y_label_) {
        return true;
      }
      if (this->y_label_ > literal.y_label_) {
        return false;
      }
      assert(this->y_label_ == literal.y_label_);
      
      if (this->edge_label_ < literal.edge_label_) {
        return true;
      }
      if (this->edge_label_ > literal.edge_label_) {
        return false;
      }
      assert(this->edge_label_ == literal.edge_label_);
      
      if (this->module_url_ < literal.module_url_) {
        return true;
      }
      return false;
     case LiteralType::kKeyLiteral:
      // x_id_ == y_id_
      if (this->x_label_ < literal.x_label_) {
        return true;
      }
      if (this->x_label_ > literal.x_label_) {
        return false;
      }
      assert(this->x_label_ == literal.x_label_);
      return this->y_label_ < literal.y_label_;
     default:
      return false;
    }
    return false;
  }

  inline bool operator==(const LiteralStandAloneInfo& literal) const {
    if (this->literal_type_ != literal.literal_type_){
      // does not have same type
      return false;
    }
    switch (this->literal_type_){
     case LiteralType::kAttrValueLiteral:
      // x.A
      if (this->x_label_    == literal.x_label_
       && this->x_attr_key_ == literal.x_attr_key_){
        return true;
      }
      return false;
     case LiteralType::kVariableLiteral:
      // x.A == y.B
      if (this->x_label_    == literal.x_label_
       && this->y_label_    == literal.y_label_
       && this->x_attr_key_ == literal.x_attr_key_
       && this->y_attr_key_ == literal.y_attr_key_){
        return true;
      }
      return false;
     case LiteralType::kConstantLiteral:
      // x.A == c
      if (this->x_label_    == literal.x_label_
       && this->x_attr_key_ == literal.x_attr_key_
       && this->c_str_      == literal.c_str_){
        return true;
      }
      return false;
     case LiteralType::kEdgeLiteral:
      // x_label_ - edge_label_ -> y_label_
      if (this->x_label_    == literal.x_label_
       && this->y_label_    == literal.y_label_
       && this->edge_label_ == literal.edge_label_){
        return true;
      }
      return false;
     case LiteralType::kMlLiteral:
      // x_id_ - (module_url_: edge_label_) -> y_id_
      if (this->x_label_    == literal.x_label_
       && this->y_label_    == literal.y_label_
       && this->module_url_ == literal.module_url_
       && this->edge_label_ == literal.edge_label_){
        return true;
      }
      return false;
     case LiteralType::kKeyLiteral:
      // x_id_ == y_id_
      if (this->x_label_ == literal.x_label_
       && this->y_label_ == literal.y_label_){
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

  inline size_t VertexNum() const {
    switch (this->literal_type_){
     case LiteralType::kAttrValueLiteral:
      // x.A
      return 1;
     case LiteralType::kVariableLiteral:
      // x.A == y.B
      return 2;
     case LiteralType::kConstantLiteral:
      // x.A == c
      return 1;
     case LiteralType::kEdgeLiteral:
      // x_label_ - edge_label_ -> y_label_
      return 2;
     case LiteralType::kMlLiteral:
      // x_id_ - (module_url_: edge_label_) -> y_id_
      return 2;
     case LiteralType::kKeyLiteral:
      // x_id_ == y_id_
      return 2;
     default:
      // unknown literal type
      assert(false);
      break;
    }
    // should not reach 
    assert(false);
    return 0;
  }

  inline bool operator!=(const LiteralStandAloneInfo& literal) const {
    return !((*this) == literal);
  }

  inline const enum LiteralType literal_type() const{
    return this->literal_type_;
  }

  inline const VertexLabelType& x_label() const{
    return this->x_label_;
  }

  inline const VertexLabelType& y_label() const{
    return this->y_label_;
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
  VertexLabelType x_label_, 
                  y_label_;

  AttributeKeyType x_attr_key_, 
                   y_attr_key_;

  std::string c_str_;

  EdgeLabelType edge_label_;

  enum GUNDAM::BasicDataType data_type_;

  std::string module_url_;
};

}  // namespace gar

#endif