#ifndef GOR_GOR_H_
#define GOR_GOR_H_

#include "gar/gar.h"

namespace gor {

// does not support additional literal now
template <typename Pattern, 
          typename DataGraph>
class GraphOracleRule : public gar::GraphAssociationRule<Pattern,
                                                       DataGraph> {
 public:
   using GraphAssociationRuleType
  = gar::GraphAssociationRule<Pattern, DataGraph>;

 public:
  using GraphAssociationRuleType::GraphAssociationRuleType;
};

}  // namespace gor

#endif