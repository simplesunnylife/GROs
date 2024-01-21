
#ifndef EXAMPLES_ANALYTICAL_APPS_TEST_TEST_CONTEXT_H_
#define EXAMPLES_ANALYTICAL_APPS_TEST_TEST_CONTEXT_H_

#include <grape/app/context_base.h>
#include <grape/grape.h>

#include <iomanip>
#include <iostream>
#include <limits>

namespace grape {

/**
 * @brief Context for the parallel version of Test.
 *
 * @tparam FRAG_T
 */
template <typename FRAG_T>
class TestContext : public VertexDataContext<FRAG_T, int64_t> {
 public:
  using oid_t = typename FRAG_T::oid_t;
  using vid_t = typename FRAG_T::vid_t;

  explicit TestContext(const FRAG_T& fragment)
      : VertexDataContext<FRAG_T, int64_t>(fragment, true) {
    return;
  }

  void Init(ParallelMessageManager& messages, std::string yaml_file) {
    this->yaml_file_ = yaml_file;
    // this->supp_.clear();
    this->fid = this->fragment().fid();
#ifdef PROFILING 
    preprocess_time = 0;
    exec_time = 0;
    postprocess_time = 0;
#endif
  }

  void Output(std::ostream& os) {
#ifdef PROFILING 
    VLOG(2) << "preprocess_time: " << preprocess_time << "s.";
    VLOG(2) << "exec_time: " << exec_time << "s.";
    VLOG(2) << "postprocess_time: " << postprocess_time << "s.";
#endif
    // if (this->fid == 0) LOG(INFO) << "last supp = " << this->supp_.size();
  }

#ifdef PROFILING
  double preprocess_time = 0;
  double exec_time = 0;
  double postprocess_time = 0;
#endif
  std::string yaml_file_;
  int fid;
};
}  // namespace grape

#endif  // EXAMPLES_ANALYTICAL_APPS_SSSP_SSSP_CONTEXT_H_
