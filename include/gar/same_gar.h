#ifndef SAME_GAR_H_
#define SAME_GAR_H_

#include "gar/gar.h"

#include "gundam/tool/same_pattern.h"
#include "gundam/match/match.h"

#include "gundam/type_getter/vertex_handle.h"

namespace gar {

namespace _gar_same{

template <typename Pattern0, typename LiteralSetType0,
          typename Pattern1, typename LiteralSetType1>
bool LiteralSetIsSame(
  Pattern0& src_pattern, const LiteralSetType0& src_literal_set,
  Pattern0& dst_pattern, const LiteralSetType1& dst_literal_set,
  const GUNDAM::Match<Pattern0, Pattern1>& match){
  if (src_literal_set.Count() != dst_literal_set.Count()){
    // the size of two set is not equal
    return false;
  }
  bool literal_set_same = true;
  for (const auto& kSrcXLiteralPtr : src_literal_set){
    bool has_same = false;
    const auto kSrcXLiteralInfo = kSrcXLiteralPtr->info();
    for (const auto& kDstXLiteralPtr : dst_literal_set){
      if (kSrcXLiteralPtr->type()
       != kDstXLiteralPtr->type()) {
        // the two literal does not have the same type
        continue;
      }
      const auto kDstXLiteralInfo = kDstXLiteralPtr->info();
      switch (kSrcXLiteralInfo.literal_type()){
      case gar::LiteralType::kVariableLiteral:
        // x.A == y.B
        if (kSrcXLiteralInfo.x_attr_key()
          != kDstXLiteralInfo.x_attr_key()){
          // the "A" of src x literal is not the
          // same as "A" of dst x literal
          break;
        }
        if (kSrcXLiteralInfo.y_attr_key()
         != kDstXLiteralInfo.y_attr_key()){
          // the "B" of src x literal is not the
          // same as "B" of dst x literal
          break;
        }
        if (match.MapTo(src_pattern.FindVertex(kSrcXLiteralInfo.x_id()))
                     != dst_pattern.FindVertex(kDstXLiteralInfo.x_id())){
          // the vertex "x" of src x literal is not the same under this match
          break;
        }
        if (match.MapTo(src_pattern.FindVertex(kSrcXLiteralInfo.y_id()))
                     != dst_pattern.FindVertex(kDstXLiteralInfo.y_id())){
          // the vertex "x" of src x literal is not the same under this match
          break;
        }
        has_same = true;
        break;
      case gar::LiteralType::kAttrValueLiteral:
        // x.A
        if (kSrcXLiteralInfo.x_attr_key()
         != kDstXLiteralInfo.x_attr_key()){
          // the "A" of src x literal is not the
          // same as "A" of dst x literal
          break;
        }
        if (match.MapTo(src_pattern.FindVertex(kSrcXLiteralInfo.x_id()))
                     != dst_pattern.FindVertex(kDstXLiteralInfo.x_id())){
          // the vertex "x" of src x literal is not the same under this match
          break;
        }
        has_same = true;
        break;
      case gar::LiteralType::kConstantLiteral:
        // x.A = c
        if (kSrcXLiteralInfo.x_attr_key()
         != kDstXLiteralInfo.x_attr_key()){
          // the "A" of src x literal is not the
          // same as "A" of dst x literal
          break;
        }
        if (kSrcXLiteralInfo.c_str()
         != kDstXLiteralInfo.c_str()){
          // the "c" of src x literal is not the
          // same as "c" of dst x literal
          break;
        }
        if (match.MapTo(src_pattern.FindVertex(kSrcXLiteralInfo.x_id()))
                     != dst_pattern.FindVertex(kDstXLiteralInfo.x_id())){
          // the vertex "x" of src x literal is not the same under this match
          break;
        }
        has_same = true;
        break;
      case gar::LiteralType::kEdgeLiteral:
        // x_id_ - edge_label_ -> y_id_
        if (kSrcXLiteralInfo.edge_label()
         != kDstXLiteralInfo.edge_label()){
          // the "edge_label" of src x literal is not the
          // same as "edge_label" of dst x literal
          break;
        }
        assert(src_pattern.FindVertex(kSrcXLiteralInfo.x_id()));
        assert(dst_pattern.FindVertex(kDstXLiteralInfo.x_id()));
        if (match.MapTo(src_pattern.FindVertex(kSrcXLiteralInfo.x_id()))
                     != dst_pattern.FindVertex(kDstXLiteralInfo.x_id())){
          // the vertex "x" of src x literal is not the same under this match
          break;
        }
        assert(src_pattern.FindVertex(kSrcXLiteralInfo.y_id()));
        assert(dst_pattern.FindVertex(kDstXLiteralInfo.y_id()));
        if (match.MapTo(src_pattern.FindVertex(kSrcXLiteralInfo.y_id()))
                     != dst_pattern.FindVertex(kDstXLiteralInfo.y_id())){
          // the vertex "x" of src x literal is not the same under this match
          break;
        }
        has_same = true;
        break;
      case gar::LiteralType::kMlLiteral:
        // x_id_ - (module_url_: edge_label_) -> y_id_
        if (kSrcXLiteralInfo.edge_label()
         != kDstXLiteralInfo.edge_label()){
          // the "edge_label" of src x literal is not the
          // same as "edge_label" of dst x literal
          break;
        }
        if (kSrcXLiteralInfo.module_url()
         != kDstXLiteralInfo.module_url()){
          // the "module_url" of src x literal is not the
          // same as "module_url" of dst x literal
          break;
        }
        assert(src_pattern.FindVertex(kSrcXLiteralInfo.x_id()));
        assert(dst_pattern.FindVertex(kDstXLiteralInfo.x_id()));
        if (match.MapTo(src_pattern.FindVertex(kSrcXLiteralInfo.x_id()))
                     != dst_pattern.FindVertex(kDstXLiteralInfo.x_id())){
          // the vertex "x" of src x literal is not the same under this match
          break;
        }
        assert(src_pattern.FindVertex(kSrcXLiteralInfo.y_id()));
        assert(dst_pattern.FindVertex(kDstXLiteralInfo.y_id()));
        if (match.MapTo(src_pattern.FindVertex(kSrcXLiteralInfo.y_id()))
                     != dst_pattern.FindVertex(kDstXLiteralInfo.y_id())){
          // the vertex "x" of src x literal is not the same under this match
          break;
        }
        has_same = true;
        break;
      case gar::LiteralType::kKeyLiteral:
        // x_id_ == y_id_
        assert(src_pattern.FindVertex(kSrcXLiteralInfo.x_id()));
        assert(dst_pattern.FindVertex(kDstXLiteralInfo.x_id()));
        if (match.MapTo(src_pattern.FindVertex(kSrcXLiteralInfo.x_id()))
                     != dst_pattern.FindVertex(kDstXLiteralInfo.x_id())){
          // the vertex "x" of src x literal is not the same under this match
          break;
        }
        assert(src_pattern.FindVertex(kSrcXLiteralInfo.y_id()));
        assert(dst_pattern.FindVertex(kDstXLiteralInfo.y_id()));
        if (match.MapTo(src_pattern.FindVertex(kSrcXLiteralInfo.y_id()))
                     != dst_pattern.FindVertex(kDstXLiteralInfo.y_id())){
          // the vertex "x" of src x literal is not the same under this match
          break;
        }
        has_same = true;
        break;
      default:
        std::cout<<"unknown x literal type!"<<std::endl;
        return false;
      }
      if (has_same){
        break;
      }
    }
    if (!has_same){
      literal_set_same = false;
      break;
    }
  }
  return literal_set_same;
}

}; // namespace _gar_same

template <typename Pattern0, typename DataGraph0,
          typename Pattern1, typename DataGraph1>
inline bool SameGar(
  const GraphAssociationRule<Pattern0, DataGraph0>& gar0,
  const GraphAssociationRule<Pattern1, DataGraph1>& gar1) {

  using namespace _gar_same;

  if (!GUNDAM::SamePattern(gar0.pattern(),
                           gar1.pattern())){
    return false;
  }

  using MatchType = GUNDAM::Match<const Pattern0,  
                                  const Pattern1>;

  std::function<bool(const MatchType&)> prune_nothing_callback
    = [](const MatchType& match_state) {
      // prune nothing, continue the matching
      return false;
    };

  bool literal_set_same = false;

  std::function<bool(const MatchType&)> verify_literal_callback
    = [&gar0, &gar1,
       &literal_set_same](const MatchType& match) {

      bool literal_set_same_in_this_match 
         = LiteralSetIsSame(gar0.pattern(), gar0.x_literal_set(),
                            gar1.pattern(), gar1.x_literal_set(),
                            match);
      if (!literal_set_same_in_this_match){
        // the x_literal_set of gar0 and gar1 is not the same
        // move to the next matching
        return true;
      }

      literal_set_same_in_this_match 
         = LiteralSetIsSame(gar0.pattern(), gar0.y_literal_set(),
                            gar1.pattern(), gar1.y_literal_set(),
                            match);
      if (!literal_set_same_in_this_match){
        // the y_literal_set of gar0 and gar1 is not the same
        // move to the next matching
        return true;
      }
      // both the x_literal_set and y_literal_set 
      // of gar0 and gar1 are the same
      // does not need to continue matching
      literal_set_same = true;
      return false;
    };

  GUNDAM::MatchUsingMatch(gar0.pattern(), 
                          gar1.pattern(),
                      prune_nothing_callback,
                     verify_literal_callback);

  return literal_set_same;
}

} // namespace gar

#endif  // SAME_GAR_H_