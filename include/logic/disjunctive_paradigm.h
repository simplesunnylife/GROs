#ifndef _LOGIC_DISJUNCTIVE_PARADIGM_H
#define _LOGIC_DISJUNCTIVE_PARADIGM_H

#include "include/util/log.h"

namespace logic {

template <typename AtomType>
class IntersectionSet {
 public:
  IntersectionSet() {
    return;
  }

  ~IntersectionSet() {
    return;
  }

  IntersectionSet(const std::vector<AtomType>& intersection_set)
                            :intersection_set_(intersection_set) {
    std::sort(this->intersection_set_.begin(),
              this->intersection_set_.end());
    assert(std::is_sorted(this->intersection_set_.begin(),
                          this->intersection_set_.end()));
    return;
  }

  IntersectionSet(const AtomType& atom)
               :intersection_set_{atom} {
    assert(std::is_sorted(this->intersection_set_.begin(),
                          this->intersection_set_.end()));
    return;
  }

  IntersectionSet(const IntersectionSet&) = default;
  IntersectionSet(IntersectionSet&&) = default;

  IntersectionSet& operator=(const IntersectionSet&) = default;	  
  IntersectionSet& operator=(IntersectionSet&&) = default;	

  // #######################
  // #  logical intersect  #
  // #######################
  //  return whether have been updated
  //  if intersection_set.intersection_set_ == {a, b}
  //            and this->intersection_set_ == {a, c}
  //  then this->intersection_set_ would be set to {a, b, c} after update
  bool Intersect(const IntersectionSet& intersection_set) {
    assert(std::is_sorted(intersection_set.intersection_set_.begin(),
                          intersection_set.intersection_set_.end()));
    std::vector<AtomType> temp_intersection_set;
    std::set_union(intersection_set .begin(), 
                   intersection_set .end(),
             this->intersection_set_.begin(), 
             this->intersection_set_.end(),
      std::back_inserter(temp_intersection_set));
    assert((temp_intersection_set .size() 
              >= intersection_set .size())
        && (temp_intersection_set .size() 
        >= this->intersection_set_.size()));
    if (this->intersection_set_.size()
      == temp_intersection_set .size()) {
      assert(std::is_sorted(this->intersection_set_.begin(),
                            this->intersection_set_.end()));
      // size are equal, does not modified, return false
      return false;
    }
    // size are not equal, modified, return true
    assert(std::is_sorted(temp_intersection_set.begin(),
                          temp_intersection_set.end()));
    this->intersection_set_.swap(temp_intersection_set);
    return true;
  }

  inline auto begin() const {
    return this->intersection_set_.begin();
  }

  inline auto end() const {
    return this->intersection_set_.end();
  }

  inline auto size() const {
    return this->intersection_set_.size();
  }

  // #####################
  // #  logical contain  #
  // #####################
  //  return true whether the logic of intersection_set is 
  //  contained in this->intersection_set_
  //    e.g. return true if intersection_set 
  //       intersects this->intersection_set_ 
  //               == this->intersection_set_
  //  example:
  //    intersection_set.intersection_set_ = {a, b, c}
  //           and this->intersection_set_ = {a, b}
  //    return false
  //
  //    intersection_set.intersection_set_ = {a, b}
  //           and this->intersection_set_ = {a, b, c}
  //    return true
  inline bool Contain(const IntersectionSet& intersection_set) const {
    if (intersection_set.size() > this->size()) {
      return false;
    }
    assert(intersection_set.size() <= this->size());
    if (!std::is_sorted(this->intersection_set_.begin(),
                        this->intersection_set_.end())) {
      for (const auto& intersection 
               : this->intersection_set_) {
        util::Debug(intersection.ToString());
      }
      for (size_t i = 0; i < this->intersection_set_.size(); i++) {
        for (size_t j = i; j < this->intersection_set_.size(); j++) {
          util::Debug("atom i: " + std::to_string(i) 
         + "is less to atom j: " + std::to_string(j) 
                           + " " + std::to_string(this->intersection_set_[i] 
                                                < this->intersection_set_[j]));
        }
      }
    }
    assert(std::is_sorted(this->intersection_set_.begin(),
                          this->intersection_set_.end()));
    assert(std::is_sorted(intersection_set.intersection_set_.begin(),
                          intersection_set.intersection_set_.end()));
                          
    // wenzhi: optimize me
    IntersectionSet temp_intersection_set = intersection_set;
    
    temp_intersection_set.Intersect(*this);
      
    assert(temp_intersection_set .size() 
             >= intersection_set .size());
    assert(temp_intersection_set .size() 
       >= this->intersection_set_.size());

    return temp_intersection_set.size()
             == intersection_set.size();
  }

  inline bool operator<(const IntersectionSet& intersection_set) const {
    if (this->intersection_set_.size() 
            < intersection_set
             .intersection_set_.size()) {
      return true;
    }
    if (this->intersection_set_.size() 
            > intersection_set
             .intersection_set_.size()) {
      return false;
    }
    assert(intersection_set.intersection_set_.size()
                   == this->intersection_set_.size());
    for (size_t idx = 0; idx < this->intersection_set_.size(); idx++) {
      if (this->intersection_set_[idx]
              < intersection_set
               .intersection_set_[idx]) {
        return true;
      }
      if (intersection_set.intersection_set_[idx]
                   < this->intersection_set_[idx]) {
        return false;
      }
    }
    // are the same
    return false;
  }

  inline std::string ToString() const {
    std::string str = "(";
    for (const auto& atom : this->intersection_set_) {
      str += " ^ " + atom.ToString();
    }
    str += ")";
    return str;
  }
  
 private:
  std::vector<AtomType> intersection_set_;
};

template <typename AtomType>
class DisjunctiveParadigm {
 private:
  /* data */
 public:
  using IntersectionSetType = IntersectionSet<AtomType>;

  DisjunctiveParadigm() = default;

  ~DisjunctiveParadigm() = default;

  DisjunctiveParadigm(const DisjunctiveParadigm&) = default;
  DisjunctiveParadigm(DisjunctiveParadigm&&) = default;

  DisjunctiveParadigm& operator=(const DisjunctiveParadigm&) = default;	  
  DisjunctiveParadigm& operator=(DisjunctiveParadigm&&) = default;	

  // #######################
  // #  logical intersect  #
  // #######################
  //  return whether have been updated
  //  intersect intersection_set with all intersection sets 
  //  contained in this->disjunctive_paradigm_
  bool Intersect(const IntersectionSetType& intersection_set) {
    if (this->disjunctive_paradigm_.empty()) {
      this->disjunctive_paradigm_.emplace(intersection_set);
      return true;
    }
    std::set<IntersectionSetType> modified_intersection_set;
    for (auto this_intersection_set_it  = this->disjunctive_paradigm_.begin();
              this_intersection_set_it != this->disjunctive_paradigm_.end();) {
      IntersectionSetType temp_intersection_set 
                       = *this_intersection_set_it;
      if (!temp_intersection_set.Intersect(intersection_set)) {
        // does not modified, continue
        this_intersection_set_it++;
        continue;
      }
      modified_intersection_set.emplace(temp_intersection_set);
      this_intersection_set_it = this->disjunctive_paradigm_.erase(this_intersection_set_it);
    }
    if (modified_intersection_set.empty()) {
      #ifndef NDEBUG
      for (const auto& this_intersection_set : this->disjunctive_paradigm_) {
        assert(this_intersection_set.Contain(intersection_set));
      }
      #endif // NDEBUG
      return false;
    }
    // modified, remove all duplicated and empty intersection set
    std::set<IntersectionSetType> modified_disjunctive_paradigm;

    std::set_union(modified_intersection_set .begin(), 
                   modified_intersection_set .end(),
                  this->disjunctive_paradigm_.begin(), 
                  this->disjunctive_paradigm_.end(),
          std::inserter(
               modified_disjunctive_paradigm,
               modified_disjunctive_paradigm.begin()));

    this->disjunctive_paradigm_.swap(modified_disjunctive_paradigm);
    return true;
  }

  bool Intersect(const AtomType& atom) {
    IntersectionSetType intersection_set(atom);
    return this->Intersect(intersection_set);
  }

  bool Intersect(const DisjunctiveParadigm& disjunctive_paradigm) {
    if (disjunctive_paradigm.Empty()) {
      return false;
    }
    bool modified = false;
    DisjunctiveParadigm  current_disjunctive_paradigm;
    current_disjunctive_paradigm.disjunctive_paradigm_
                    .swap( this->disjunctive_paradigm_ );

    for (const auto& intersection_set  : disjunctive_paradigm) {
      DisjunctiveParadigm temp_disjunctive_paradigm
                     = current_disjunctive_paradigm;
      temp_disjunctive_paradigm.Intersect(intersection_set);
      this->Union(temp_disjunctive_paradigm);
    }
    return modified;
  }
  
  // first find whether intersection_set is contained in any existed
  // intersection_set
  bool Union(const IntersectionSetType& intersection_set) {
    bool contained_in_an_intersection_set = false;
    for (auto this_intersection_set_it  = this->disjunctive_paradigm_.begin();
              this_intersection_set_it != this->disjunctive_paradigm_.end();) {
      if (!this_intersection_set_it->Contain(intersection_set)) {
        if (!intersection_set.Contain(*this_intersection_set_it)) {
          // this_intersection_set does not needs to be removed
          this_intersection_set_it++;
          continue;
        }
        #ifndef NDEBUG 
        for (const auto& this_intersection_set : this->disjunctive_paradigm_) {
          assert(!this_intersection_set.Contain(intersection_set));
        }
        #endif // NDEBUG 

        // this_intersection_set needs to be removed
        this_intersection_set_it = this->disjunctive_paradigm_.erase(this_intersection_set_it);
        continue;
      }
      contained_in_an_intersection_set = true;
      break;
    }
    if (contained_in_an_intersection_set) {
      // does not need to modify this set
      return false;
    }
    // is not contained in any intersection_set 
    // in this->disjunctive_paradigm_ 
    auto [ disjunctive_paradigm_it,
           disjunctive_paradigm_ret ] = this->disjunctive_paradigm_.emplace(intersection_set);
    assert(disjunctive_paradigm_ret);
    return true;
  }

  bool Union(const DisjunctiveParadigm& disjunctive_paradigm) {
    if (disjunctive_paradigm.Empty()) {
      if (this->Empty()) {
        return false;
      }
      this->Clear();
      return true;
    }
    bool modified = false;
    for (const IntersectionSetType& intersection_set : disjunctive_paradigm) {
      modified |= this->Union(intersection_set);
    }
    return modified;
  }

  inline auto begin() const {
    return this->disjunctive_paradigm_.begin();
  }

  inline auto end() const {
    return this->disjunctive_paradigm_.end();
  }

  inline std::string ToString() const {
    std::string str;
    for (const auto& intersection_set : this->disjunctive_paradigm_) {
      str += " v " + intersection_set.ToString();
    }
    return str;
  }

  inline size_t Size() const {
    return this->disjunctive_paradigm_.size();
  }

  inline bool Empty() const {
    return this->disjunctive_paradigm_.empty();
  }

  inline void Clear() {
    return this->disjunctive_paradigm_.clear();
  }

 private:
  std::set<IntersectionSetType> disjunctive_paradigm_;
};

}; // logic

#endif // _LOGIC_DISJUNCTIVE_PARADIGM_H