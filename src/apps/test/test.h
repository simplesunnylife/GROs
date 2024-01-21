
#ifndef EXAMPLES_ANALYTICAL_APPS_TEST_TEST_H_
#define EXAMPLES_ANALYTICAL_APPS_TEST_TEST_H_

#include <grape/grape.h>

#include "../fragment2gundam.h"
#include "../fragment_graph.h"
#include "../fragment_graph_with_index.h"
#include "gar/csv_gar.h"
#include "gar/gar.h"
#include "gundam/algorithm/dp_iso.h"
#include "gundam/algorithm/bfs.h"
#include "gundam/algorithm/simulation/bisimulation.h"
#include "gundam/graph_type/large_graph2.h"
#include "gundam/tool/operator/duplicate_vertex.h"
#include "gundam/tool/topological/tree/to_tree.h"
#include "test/test_context.h"

#include "gundam/io/rapidcsv.h"

namespace grape {

template <typename FRAG_T>
class Test : public ParallelAppBase<FRAG_T, TestContext<FRAG_T>>,
             public ParallelEngine {
 public:
  // specialize the templated worker.
  INSTALL_PARALLEL_WORKER(Test<FRAG_T>, TestContext<FRAG_T>, FRAG_T)
  using vertex_t = typename fragment_t::vertex_t;

  template <class policy = GUNDAM::execution::sequenced_policy,
            typename  QueryGraph, 
            typename TargetGraph,
            class CandidateSetContainer>
  inline bool RefineCandidateSetPartByPart(QueryGraph  &query_graph,
                                          TargetGraph &target_graph,
                                CandidateSetContainer &candidate_set,
                                         const size_t &each_parition_size,
                                    const std::string &outfile_str) {

    static_assert(std::is_same_v<policy, GUNDAM::execution::sequenced_policy>
               || std::is_same_v<policy, GUNDAM::execution:: parallel_policy>, "illegal policy");

    const auto kWayNum = 1000;

    using QueryVertexHandle = typename GUNDAM::VertexHandle<QueryGraph>::type;

    // auto begin_time = std::chrono::high_resolution_clock::now();

    // std::vector<QueryVertexHandle> topo_seq;
    // GUNDAM::_dp_iso::_DAGDP::GetTopoSeq(query_graph, topo_seq);

    // std::cout << "get topo seq time: "
    //           << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin_time).count()
    //           << " (ms)" << std::endl;
    // assert(topo_seq.size() == target_graph.CountVertex());

    auto begin_time = std::chrono::high_resolution_clock::now();

    std::vector<
    std::vector<QueryVertexHandle>> multi_topo_seq;
    GUNDAM::_dp_iso_using_match::_DAGDP::GetMultiWayTopoSeq(query_graph, multi_topo_seq, kWayNum);

    std::cout << "get multi topo seq time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin_time).count()
              << " (ms)" << std::endl;
    assert(multi_topo_seq.size() <= target_graph.CountVertex());

    for (const auto& batch : multi_topo_seq) {
      std::cout << "\t" << batch.size();
    }
    std::cout << std::endl;

    std::map<QueryVertexHandle, size_t> candidate_counter_set;
    std::set<QueryVertexHandle> used_vertex_set;

    const size_t kEachParitionNum = (multi_topo_seq.size() + each_parition_size - 1) 
                                                           / each_parition_size;

    std::vector<uint64_t> candidate_size;

    uint64_t temp_candidate_counter = 0;
    for (const auto& [query_vertex, candidate] : candidate_set) {
      candidate_counter_set[query_vertex] = candidate.size();
      temp_candidate_counter += candidate.size();
    }
    candidate_size.emplace_back(temp_candidate_counter);

    std::vector<QueryVertexHandle> partition_vertex_handle;
    partition_vertex_handle.reserve(each_parition_size + kWayNum);

    size_t   parition_idx = 0;
    size_t last_batch_idx = 0;

    for (size_t batch_idx = 0;
                batch_idx < multi_topo_seq.size();
                batch_idx++) {

      partition_vertex_handle.insert(partition_vertex_handle.end(), 
                                   multi_topo_seq[batch_idx].begin(), 
                                   multi_topo_seq[batch_idx].end());

      const bool kIsLastRound = (batch_idx == (multi_topo_seq.size() - 1));

      if (!kIsLastRound // is not the last batch
        && partition_vertex_handle.size() < each_parition_size) {
        continue;
      }

      const size_t kCurrentPartitionBegin = last_batch_idx;
      const size_t kCurrentPartitionEnd   = batch_idx + 1;

      auto begin_time = std::chrono::high_resolution_clock::now();
      if (!GUNDAM::_dp_iso_using_match::InitCandidateSet<
           GUNDAM::MatchSemantics::kIsomorphism, policy>(query_graph,
                                                         target_graph,
                                               partition_vertex_handle.begin(),
                                               partition_vertex_handle.end(),
                                                      candidate_set)) {
        return false;
      }
      std::cout << "InitCandidateSet time of " << parition_idx << "'th partition: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin_time).count()
                << " (ms)" << std::endl;

      for (const auto& vertex_handle : partition_vertex_handle) {
        used_vertex_set.emplace(vertex_handle);
      }

      uint64_t temp_candidate_counter = 0;
      for (const auto& [query_vertex, candidate] : candidate_set) {
        candidate_counter_set[query_vertex] = candidate.size();
        temp_candidate_counter += candidate.size();
      }
      candidate_size.emplace_back(temp_candidate_counter);

      for (int i = 0; ; i++) {
        auto begin_time = std::chrono::high_resolution_clock::now();

        if (i == 0) {
          if (!GUNDAM::_dp_iso_using_match
                     ::ParDAGDP(query_graph, 
                               target_graph, 
                                multi_topo_seq.begin() + kCurrentPartitionBegin, 
                                multi_topo_seq.begin() + kCurrentPartitionEnd, 
                            candidate_set,
                                 used_vertex_set)) {
            return false;
          }
        }
        else {
          if (!GUNDAM::_dp_iso_using_match
                     ::ParDAGDP(query_graph, 
                               target_graph, 
                                multi_topo_seq.begin(), 
                                multi_topo_seq.begin() + kCurrentPartitionEnd, 
                            candidate_set)) {
            return false;
          }
        }

        std::cout << "        loop_num " << i << " on first " << parition_idx << "'th partition: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin_time).count()
                  << " (ms)" << std::endl;

        bool converged = true;
        uint64_t temp_candidate_counter = 0;
        for (const auto& [query_vertex, candidate] : candidate_set) {
          if (candidate_counter_set[query_vertex] != candidate.size()) {
            converged = false;
          }
          candidate_counter_set[query_vertex] = candidate.size();
          temp_candidate_counter += candidate.size();
        }
        candidate_size.emplace_back(temp_candidate_counter);
        if (converged) {
          break;
        }

        begin_time = std::chrono::high_resolution_clock::now();

        if (!GUNDAM::_dp_iso_using_match::ParDAGDP(query_graph, 
                                                  target_graph, 
                        std::make_reverse_iterator(multi_topo_seq.begin() + kCurrentPartitionEnd - 1), 
                                                   multi_topo_seq.rend(), 
                                              candidate_set)) {
          return false;
        }

        std::cout << "reverse loop_num " << i << " on first " << parition_idx << "'th partition: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin_time).count()
                  << " (ms)" << std::endl;

        converged = true;
        temp_candidate_counter = 0;
        for (const auto& [query_vertex, candidate] : candidate_set) {
          if (candidate_counter_set[query_vertex] != candidate.size()) {
            converged = false;
          }
          candidate_counter_set[query_vertex] = candidate.size();
          temp_candidate_counter += candidate.size();
        }
        candidate_size.emplace_back(temp_candidate_counter);
        if (converged) {
          break;
        }

        std::ofstream outlog(outfile_str);
        for (const auto& candidate_counter : candidate_size) {
          outlog << candidate_counter << std::endl;
        }
      }

      std::ofstream outlog(outfile_str);
      for (const auto& candidate_counter : candidate_size) {
        outlog << candidate_counter << std::endl;
      }

      for (auto& [query_vertex, candidate] : candidate_set) {
        candidate.shrink_to_fit();
      }

      parition_idx++;
      last_batch_idx = kCurrentPartitionEnd;
      partition_vertex_handle.clear();
    }

    for (int i = 0; ; i++) {
      auto begin_time = std::chrono::high_resolution_clock::now();

      if (!GUNDAM::_dp_iso_using_match::ParDAGDP(query_graph, 
                                                target_graph, 
                                                 multi_topo_seq, 
                                                  candidate_set)) {
        return false;
      }

      std::cout << "loop_num " << i << " : "
                << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin_time).count()
                << " (ms)" << std::endl;

      bool converged = true;
      uint64_t temp_candidate_counter = 0;
      for (const auto& [query_vertex, candidate] : candidate_set) {
        if (candidate_counter_set[query_vertex] != candidate.size()) {
          converged = false;
        }
        candidate_counter_set[query_vertex] = candidate.size();
        temp_candidate_counter += candidate.size();
      }
      candidate_size.emplace_back(temp_candidate_counter);
      if (converged) {
        return true;
      }

      std::ofstream outlog(outfile_str);
      for (const auto& candidate_counter : candidate_size) {
        outlog << candidate_counter << std::endl;
      }
      std::reverse(multi_topo_seq.begin(), 
                   multi_topo_seq.end());
    }
    return true;
  }

  template <typename  QueryGraph, 
            typename TargetGraph,
            class CandidateSetContainer>
  inline bool RefineCandidateSetWithTimming(QueryGraph  &query_graph,
                                           TargetGraph &target_graph,
                                 CandidateSetContainer &candidate_set,
                                     const std::string& outfile_str) {
    using QueryVertexHandle = typename GUNDAM::VertexHandle<QueryGraph>::type;
    std::vector<QueryVertexHandle> topo_seq;
    GUNDAM::_dp_iso::_DAGDP::GetTopoSeq(query_graph, topo_seq);
    std::map<QueryVertexHandle, 
        std::vector<size_t>> candidate_counter_set;

    for (const auto& [query_vertex, candidate] : candidate_set) {
      candidate_counter_set[query_vertex].emplace_back(candidate.size());
    }

    constexpr int loop_num = 1000;
    for (int i = 0; i <= loop_num; i++) {
      auto begin_time = std::chrono::high_resolution_clock::now();

      if (!GUNDAM::_dp_iso::DAGDP(query_graph, target_graph, topo_seq, candidate_set))
        return false;

      std::cout << "loop_num " << i << " : "
                << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin_time).count()
                << " (ms)" << std::endl;

      bool converged = true;
      for (const auto& [query_vertex, candidate] : candidate_set) {
        if (candidate_counter_set[query_vertex].back() != candidate.size()) {
          converged = false;
        }
        candidate_counter_set[query_vertex].emplace_back(candidate.size());
      }
      if (converged) {
        return true;
      }

      std::ofstream outlog(outfile_str);
      for (const auto& [query_vertex, candidate_counter] 
                                    : candidate_counter_set) {
        outlog << query_vertex->id();
        for (const auto& counter : candidate_counter) {
          outlog << "," << counter;
        }
        outlog << std::endl;
      }

      std::reverse(topo_seq.begin(), topo_seq.end());
    }
    return true;
  }


  template <typename GraphPatternType,
            typename    DataGraphType>
  void WriteTriple(const std::string kGarSetName){
    std::vector<gar::GraphAssociationRule<GraphPatternType,
                                             DataGraphType>> gar_list;
    using GarType = gar::GraphAssociationRule<GraphPatternType,
                                                 DataGraphType>;
                                                 
    using PatternWithAttr = GUNDAM::LargeGraph2<int, int, std::string, 
                                                int, int, std::string>;

    gar::ReadGARSet(gar_list, kGarSetName + "_v_set.csv",
                              kGarSetName + "_e_set.csv",
                              kGarSetName + "_x_set.csv",
                              kGarSetName + "_y_set.csv");

    const std::string kTripleKey = "triple";

    int counter = 0;
    for (const auto& gar : gar_list){
        PatternWithAttr gar_with_triple(gar.pattern());
        for (auto edge_it = gar_with_triple.EdgeBegin();
                 !edge_it.IsDone();
                  edge_it++){
            const std::string kTriple = GUNDAM::ToString(edge_it->src_handle()->label())
                                + " " + GUNDAM::ToString(edge_it->label())
                                + " " + GUNDAM::ToString(edge_it->dst_handle()->label());
            edge_it->AddAttribute(kTripleKey, kTriple);
        }
        gar::WriteGAR(gar, kGarSetName + "_" + std::to_string(counter) + "_v.csv",
                           kGarSetName + "_" + std::to_string(counter) + "_e.csv",
                           kGarSetName + "_" + std::to_string(counter) + "_x.csv",
                           kGarSetName + "_" + std::to_string(counter) + "_y.csv");
        GUNDAM::WriteCSVGraph(gar_with_triple, 
                             kGarSetName + "_" + std::to_string(counter) + "_v.csv",
                             kGarSetName + "_" + std::to_string(counter) + "_e.csv");
        counter++;
    }  
  }

  template <typename GraphType>
  void MovieLensToGraph(const std::string& rating_file,
                        const std::string&  movie_file,
                                GraphType& graph,
                                      bool make_id_continued = false) {

    rapidcsv::Document  movie_csv_file( movie_file, rapidcsv::LabelParams(0, -1));
    // phase movie column names
    std::vector<std::string> movie_col_name = movie_csv_file.GetColumnNames();

    using VertexIDType    = typename GraphType::VertexType::IDType;
    using VertexLabelType = typename GraphType::VertexType::LabelType;

    using EdgeIDType    = typename GraphType::EdgeType::IDType;
    using EdgeLabelType = typename GraphType::EdgeType::LabelType;

    const VertexLabelType kMovieLabel = 0,
                           kUserLabel = 1;

    const EdgeLabelType kHighReviewLabel = 1,
                         kLowReviewLabel = 0;

    const std::string kMovieIdColName = "movieId",
                        kTitleColName = "title",
                       kGenresColName = "genres";

    std::vector<VertexIDType>  
        movie_id     = movie_csv_file.GetColumn<VertexIDType>(kMovieIdColName);
    std::vector<std::string>  
        movie_title  = movie_csv_file.GetColumn<std::string>(kTitleColName);
    std::vector<std::string>  
        movie_genres = movie_csv_file.GetColumn<std::string>(kGenresColName);

    size_t movie_size = movie_id.size();

    VertexIDType max_movie_id = 0;

    int count_success = 0;
    int count_fail = 0;
    
    for (size_t row = 0; row < movie_size; row++) {
      auto this_movit_id = movie_id[row];
      max_movie_id = max_movie_id > movie_id[row]?
                     max_movie_id : movie_id[row];
      auto [vertex_handle, 
            vertex_ret] = graph.AddVertex(movie_id[row], kMovieLabel);
      if (!vertex_ret) {
        // vertex added fail
        ++count_fail;
        continue;
      }

      // vertex added successfully
      auto [title_handle,
            title_ret] = vertex_handle->AddAttribute(kTitleColName, movie_title[row]);
      if (!title_ret){
        // vertex title added fail
        ++count_fail;
        continue;
      }

      std::istringstream ss( movie_genres[row] );
      std::set <std::string> genres;
      // genres.reserve(20);
      while (ss) {
        std::string s;
        if (!std::getline( ss, s, '|' )) {
          break;
        }
        genres.emplace( s );
      }

      std::string all_genres;
      for (const auto& genre : genres){
        all_genres = all_genres + "|" + genre;
      }
      auto [genre_handle,
            genre_ret] = vertex_handle->AddAttribute(kGenresColName, all_genres);
      if (!genre_ret){
        // vertex title added fail
        ++count_fail;
      }

      // bool add_genre_failed = false;
      // for (const auto& genre : genres){
      //   auto [genre_handle,
      //         genre_ret] = vertex_handle->AddAttribute(genre, 1);
      //   if (!genre_ret){
      //     // vertex title added fail
      //     ++count_fail;
      //     add_genre_failed = true;
      //     break;
      //   }
      // }
    }
    movie_id.clear();
    movie_title.clear();
    movie_genres.clear();

    rapidcsv::Document rating_csv_file(rating_file, rapidcsv::LabelParams(0, -1));

    const std::string kRatingUserIdColName = "userId",
                     kRatingMovieIdColName = "movieId",
                      kRatingRatingColName = "rating",
                   kRatingTimestampColName = "timestamp";

    const std::vector<VertexIDType>  
        rating_user_id  = rating_csv_file.GetColumn<VertexIDType>(kRatingUserIdColName);
    const std::vector<VertexIDType>  
        rating_movie_id = rating_csv_file.GetColumn<VertexIDType>(kRatingMovieIdColName);
    const std::vector<float>  
        rating_rating   = rating_csv_file.GetColumn<float>(kRatingRatingColName);

    size_t rating_size = rating_user_id.size();

    std::set<VertexIDType> user_id_set;

    for (const auto& user_id : rating_user_id) {
      user_id_set.emplace(user_id);
    }

    VertexIDType user_id_offset = max_movie_id + 1 - *(user_id_set.begin());

    for (const auto& user_id : user_id_set) {
      auto [vertex_handle, 
            vertex_ret] = graph.AddVertex(user_id + user_id_offset, kUserLabel);
      if (!vertex_ret) {
        // vertex added fail
        ++count_fail;
        continue;
      }
    }
    user_id_set.clear();

    EdgeIDType edge_id_counter = 0;
    for (size_t row = 0; row < rating_size; row++) {
      constexpr float kRatingThreadhold = 4.0;
      auto [edge_handle, 
            edge_ret] = graph.AddEdge(rating_user_id [row] + user_id_offset,
                                      rating_movie_id[row], 
                                      // rating_rating  [row] * 2,
                                      rating_rating  [row] >= kRatingThreadhold ? 
                                                              kHighReviewLabel 
                                                             : kLowReviewLabel,
                                      edge_id_counter++);
      if (!edge_ret) {
        // edge added fail
        ++count_fail;
        continue;
      }
    }

    if (count_fail > 0) {
      std::cout << "Failed: " << count_fail << std::endl;
    }
    return;
  }

  template <typename GraphType>
  void YelpToGraph(const std::string& business_file,
                   const std::string&     user_file,
                   const std::string&   review_file,
                   const std::string&     tips_file,
                           GraphType& graph,
                   const std::string& city) {

    using VertexIDType    = typename GraphType::VertexType::IDType;
    using VertexLabelType = typename GraphType::VertexType::LabelType;

    using EdgeIDType    = typename GraphType::EdgeType::IDType;
    using EdgeLabelType = typename GraphType::EdgeType::LabelType;

    rapidcsv::Document  business_csv_file( business_file, rapidcsv::LabelParams(0, -1));
    // phase business column names
    std::vector<std::string> business_col_name = business_csv_file.GetColumnNames();


    const VertexLabelType kBusinessLabel = 0,
                              kUserLabel = 1;

    const std::string kBusinessBusinessIdColName = "business_id",
                      kBusinessCategoriesColName = "categories",
                            kBusinessCityColName = "city",
                            kBusinessNameColName = "name",
                           kBusinessStarsColName = "stars";

    std::vector<std::string>  business_id
                            = business_csv_file.GetColumn<std::string>(kBusinessBusinessIdColName);
    std::vector<std::string>  business_categories 
                            = business_csv_file.GetColumn<std::string>(kBusinessCategoriesColName);
    std::vector<std::string>  business_city       
                            = business_csv_file.GetColumn<std::string>(kBusinessCityColName);
    std::vector<std::string>  business_name       
                            = business_csv_file.GetColumn<std::string>(kBusinessNameColName);
    std::vector<float>        business_stars 
                            = business_csv_file.GetColumn<float>(kBusinessStarsColName);

    std::map<std::string, int> business_city_counter;
    for (const auto& city : business_city){
      business_city_counter[city]++;
    }
    for (const auto& city_counter : business_city_counter){
      std::cout << "\t" << city_counter.first << "\t" << city_counter.second << std::endl;
    }

    std::map<std::string, int> business_id_map;
    assert(business_id.size() == business_categories.size());
    assert(business_id.size() == business_city.size());
    assert(business_id.size() == business_name.size());
    assert(business_id.size() == business_stars.size());

    for (size_t row = 0; row < business_id.size(); row ++) {
      if (row % 100000 == 0){
        std::cout << "added business number: " << row << std::endl;
      }

      std::size_t found = business_city[row].find(city);
      if (found == std::string::npos) {
        // this business is not in the specified city
        std::cout << "business id: " << business_id[row] 
                  << " is not selected"
                  << std::endl;
        continue;
      }
      std::cout << "business id: " << business_id[row] 
                << " is selected"
                << std::endl;

      const VertexIDType kVertexId = business_id_map.size(); 
      auto [ vertex_handle,
             vertex_ret ]
           = graph.AddVertex(kVertexId, kBusinessLabel);
      assert(vertex_ret);
      assert(vertex_handle);

      vertex_handle->AddAttribute(kBusinessBusinessIdColName, business_id[row]);
      vertex_handle->AddAttribute(kBusinessCategoriesColName, business_categories[row]);
      vertex_handle->AddAttribute(      kBusinessCityColName, business_city[row]);
      vertex_handle->AddAttribute(      kBusinessNameColName, business_name[row]);
      vertex_handle->AddAttribute(     kBusinessStarsColName, business_stars[row]);
      business_id_map.emplace(business_id[row], kVertexId);
    }
    business_id.clear();
    business_categories.clear();
    business_city.clear();
    business_name.clear();
    business_stars.clear();
    std::cout << "# business added #" << std::endl;

    rapidcsv::Document  review_csv_file( review_file, rapidcsv::LabelParams(0, -1));
    // phase review column names
    std::vector<std::string> review_col_name = review_csv_file.GetColumnNames();
    for (const auto& col_name : review_col_name){
      std::cout << "col_name : " << col_name << std::endl;
    }

    const std::string kReviewBusinessIdColName = "business_id",
                          kReviewUserIdColName =     "user_id",
                        kReviewReviewIdColName =   "review_id",
                           kReviewStarsColName =       "stars";

    const std::string kUserFriendsColName = "friends",
                       kUserUserIdColName = "user_id";

    std::vector<std::string> review_user_id
                           = review_csv_file.GetColumn<std::string>(kReviewUserIdColName);
    std::vector<std::string> review_business_id
                           = review_csv_file.GetColumn<std::string>(kReviewBusinessIdColName);
    std::vector<std::string> review_id
                           = review_csv_file.GetColumn<std::string>(kReviewReviewIdColName);
    std::vector<double> review_stars
                      = review_csv_file.GetColumn<double>(kReviewStarsColName);

    // const EdgeLabelType kReviewEdgeLabel = 1;
    const EdgeLabelType kHighReviewEdgeLabel = 3,
                         kLowReviewEdgeLabel = 4;
    
    EdgeIDType edge_id_counter = 0;

    std::map<std::string, int> user_id_map;
    for (size_t row = 0; row < review_user_id.size(); row ++) {
      if (row % 100000 == 0){
        std::cout << "added review: " << row << std::endl;
      }

      const std::string& kReviewBusinessId = review_business_id[row];
      auto review_business_id_it  = business_id_map.find(kReviewBusinessId);
      if ( review_business_id_it == business_id_map.end()){
        std::cout << "business id: " << kReviewBusinessId 
                  << " is not selected"
                  << std::endl;
        continue;
      }

      const std::string& kReviewUserId = review_user_id[row];
      auto [review_user_id_it,
            review_user_id_ret] = user_id_map.emplace(kReviewUserId, business_id_map.size()
                                                                       + user_id_map.size());
      if (review_user_id_ret){
        // this user has not been considered before
        auto [vertex_handle,
              vertex_ret] = graph.AddVertex(review_user_id_it->second, kUserLabel);
        vertex_handle->AddAttribute(kUserUserIdColName, kReviewUserId);
      }


      constexpr float kReviewStarThreadhold = 4.0;

      auto [ edge_handle,
             edge_ret ] = graph.AddEdge(review_user_id_it->second, 
                                    review_business_id_it->second, 
                                   (review_stars[row] * 2) + kLowReviewEdgeLabel,
                                    // review_stars[row] >= kReviewStarThreadhold ? 
                                    //                       kHighReviewEdgeLabel 
                                    //                      : kLowReviewEdgeLabel, 
                                     edge_id_counter++);

      assert(edge_ret);

      const std::string& kReviewId = review_id[row];
      auto [ attr_handle,
             attr_ret ] = edge_handle->AddAttribute(kReviewReviewIdColName, kReviewId);
    }
    review_user_id.clear();
    review_business_id.clear();
    review_id.clear();

    rapidcsv::Document  tips_csv_file( tips_file, rapidcsv::LabelParams(0, -1));
    // phase tips column names
    std::vector<std::string> tips_col_name = tips_csv_file.GetColumnNames();

    const std::string kTipsBusinessIdColName = "business_id",
                          kTipsUserIdColName =     "user_id";

    std::vector<std::string> tips_user_id
                           = tips_csv_file.GetColumn<std::string>(kTipsBusinessIdColName);
    std::vector<std::string> tips_business_id
                           = tips_csv_file.GetColumn<std::string>(kTipsUserIdColName);

    const EdgeLabelType kTipsEdgeLabel = 2;

    for (size_t row = 0; row < tips_user_id.size(); row ++) {
      if (row % 100000 == 0){
        std::cout << "added tips: " << row << std::endl;
      }

      const std::string& kTipsBusinessId = tips_business_id[row];
      auto tips_business_id_it = business_id_map.find(kTipsBusinessId);
      if (tips_business_id_it == business_id_map.end()){
        std::cout << "business id: " << kTipsBusinessId 
                  << " is not selected"
                  << std::endl;
        continue;
      }

      const std::string& kTipsUserId = tips_user_id[row];
      auto [tips_user_id_it,
            tips_user_id_ret] = user_id_map.emplace(kTipsUserId, business_id_map.size()
                                                                   + user_id_map.size());
      if (tips_user_id_ret){
        // this user has not been considered before
        auto [vertex_handle,
              vertex_ret] = graph.AddVertex(tips_user_id_it->second, kUserLabel);
        vertex_handle->AddAttribute(kUserUserIdColName, kTipsUserId);
      }

      auto [ edge_handle,
             edge_ret ] = graph.AddEdge(tips_user_id_it->second, 
                                    tips_business_id_it->second, kTipsEdgeLabel, edge_id_counter++);
      assert(edge_ret);
    }
    tips_user_id.clear();
    tips_business_id.clear();

    rapidcsv::Document  user_csv_file( user_file, rapidcsv::LabelParams(0, -1));
    // phase user column names
    std::vector<std::string> user_col_name = user_csv_file.GetColumnNames();

    const EdgeLabelType kFriendEdgeLabel = 0;
    
    std::vector<std::string> user_id
                           = user_csv_file.GetColumn<std::string>(kUserUserIdColName);
    std::vector<std::string> user_friends
                           = user_csv_file.GetColumn<std::string>(kUserFriendsColName);

    for (size_t row = 0; row < user_friends.size(); row ++) {
      if (row % 100000 == 0){
        std::cout << "added user friends: " << row << std::endl;
      }
      const std::string& kUserId = user_id     [row];
      const std::string& friends = user_friends[row];

      auto user_id_it =  user_id_map.find(kUserId);
      if ( user_id_it == user_id_map.end()){
        std::cout << "user id: " << kUserId 
                  << " is not review/tips to seleceted business"
                  << std::endl;
        continue;
      }
      const VertexIDType kUserVertexID = user_id_it->second;

      std::istringstream ss( friends );
      while (ss) {
        std::string s;
        if (!std::getline( ss, s, ',' )) {
          break;
        }
        s.erase(std::remove_if(s.begin(), s.end(),
                               [](unsigned char x){return std::isspace(x);}),
                s.end());
        s.erase(std::remove_if(s.begin(), s.end(),
                               [](unsigned char x){return x == '"';}),
                s.end());

        auto [user_id_map_it,
              user_id_map_ret] = user_id_map.emplace(s, user_id_map.size());

        auto friend_id_it =  user_id_map.find(s);
        if ( friend_id_it == user_id_map.end()){
          std::cout << "friend user id: " << kUserId 
                    << " is not review/tips to seleceted business"
                    << std::endl;
          continue;
        }
        
        const VertexIDType kFriendVertexId = friend_id_it->second;

        auto [ edge_handle0,
               edge_ret0 ] = graph.AddEdge(kUserVertexID, kFriendVertexId, kFriendEdgeLabel, edge_id_counter++);
        assert(edge_ret0);
        auto [ edge_handle1,
               edge_ret1 ] = graph.AddEdge(kFriendVertexId, kUserVertexID, kFriendEdgeLabel, edge_id_counter++);
        assert(edge_ret1);
      }
    }
    user_id.clear();
    return;
  }

  inline std::vector<std::string>
    LoadGarPath(const std::string& gar_dir,
                const std::string& rule_prefix,
                const std::string& rule_posfix) const {

    std::set<std::string> dir_files;
    std::vector<std::string> qualified_gar_path_set_in_one_dir;

    util::GetFiles(gar_dir, dir_files);
    for (auto gar_file : dir_files){
      const bool isPrefix = rule_prefix.size() 
                            <= gar_file.size() 
           && std::mismatch(rule_prefix.begin(), 
                            rule_prefix. end (),
                               gar_file.begin(), 
                               gar_file. end ()).first == rule_prefix.end();
      const bool isGarPosfix = rule_posfix.size() 
                               <= gar_file.size() 
              && std::mismatch(rule_posfix.begin(), 
                               rule_posfix. end (),
                                  gar_file.end() - rule_posfix.size(), 
                                  gar_file.end()).first == rule_posfix.end();
      if (!isPrefix || !isGarPosfix){
        continue;
      }

      // load seperated single gar
      gar_file = gar_file.substr(0, gar_file.length()
                                - rule_posfix.size());

      const std::string kGarName = gar_dir + "/" + gar_file;

      qualified_gar_path_set_in_one_dir.emplace_back(kGarName);
    }
    return qualified_gar_path_set_in_one_dir;
  }

  inline std::vector<std::string>
    LoadGarPath(const std::string& gar_dir,
                const std::string& rule_prefix) const {
    return LoadGarPath(gar_dir,
                       rule_prefix, "_v.csv");
  }

  inline std::vector<std::string>
    LoadGarSetPath(const std::string& gar_dir,
                   const std::string& rule_prefix) const {
    return LoadGarPath(gar_dir,
                       rule_prefix, "_v_set.csv");
  }

  template <typename GraphType>
  void RoundAttrToHalf(const std::string&  input_v_file,
                       const std::string&  input_e_file,
                       const std::string& vertex_attr_key,
                       const std::string& output_v_file,
                       const std::string& output_e_file) {
    GraphType data_graph;
    GUNDAM::ReadCSVGraph(data_graph, input_v_file,
                                     input_e_file);
    for (auto vertex_it = data_graph.VertexBegin();
             !vertex_it.IsDone();
              vertex_it++) {
      auto attr_handle = vertex_it->FindAttribute(vertex_attr_key);
      if (!attr_handle) {
        // does not has this attribute
        continue;
      }
      attr_handle->template value<double>() = std::llround( attr_handle->template value<double>() * 2 ) / 2.0 ;
    }
    GUNDAM::WriteCSVGraph(data_graph, output_v_file,
                                      output_e_file);
    return;
  }

  inline void GetData(std::string dir,                  // "/home/fwz/JedAIToolkit/data/cleanCleanErDatasets/DMDG/DBLP/"
                      std::string from_graph_csv_name,  // "q3_DM_DG"
                      std::string      graph_id_name,   // "id_DM"
                      std::string      graph_attr_name, // "volumn"
                      std::string from_table_csv_name,  // "output_article"
                      std::string      table_id_name,   // "id"
                      std::string      table_attr_name) // "volume"
                      {

    rapidcsv::Document q_DM_DG_csv(dir + "/" + from_graph_csv_name + ".csv",
                                   rapidcsv::LabelParams(0, -1));

    std::vector<std::string> q_DM_DG_col_name = q_DM_DG_csv.GetColumnNames();

    std::cout << "q_DM_DG col name: " << std::endl;

    for (auto it  = q_DM_DG_col_name.begin(); 
              it != q_DM_DG_col_name.end();
              it++) {
      std::cout << "\t" << *it;
    }
    std::cout << std::endl;

    std::vector<std::string> 
      id_DM = q_DM_DG_csv.GetColumn<std::string>(graph_id_name);
    std::vector<std::string> 
      volumn = q_DM_DG_csv.GetColumn<std::string>(graph_attr_name);

    for (auto& id_it : id_DM) {
      if (id_it.front() != '"'
       || id_it.back () != '"'){
        continue;
      }
      id_it.erase(0, 1);
      id_it.erase(id_it.size() - 1);
    }

    std::sort(id_DM.begin(), id_DM.end());

    std::map<std::string,
             std::string> id_volumn_map;

    for (int i = 0; i < id_DM.size(); i++) {
      id_volumn_map.emplace(id_DM[i], volumn[i]);
    }

    std::cout << "id_DM.size(): "
              <<  id_DM.size() << std::endl;

    rapidcsv::Document output_article_csv(dir + "/" + from_table_csv_name + ".csv",
                                          rapidcsv::LabelParams(0, -1));

    std::vector<std::string> output_article_col_name 
                           = output_article_csv.GetColumnNames();

    std::cout << "output_article col name: " << std::endl;

    for (auto it  = output_article_col_name.begin(); 
              it != output_article_col_name.end();
              it++) {
      std::cout << *it << std::endl;
    }

    rapidcsv::Document q_DM_DG_output_csv (output_article_csv);

    std::vector<int> idx_to_reserve;

    const std::vector<std::string> 
      id = q_DM_DG_output_csv.GetColumn<std::string>(table_id_name);

    const std::vector<std::string> q_DM_DG_output_csv_column_name
                                 = q_DM_DG_output_csv.GetColumnNames();

    const int kIdColumnIdx = q_DM_DG_output_csv.GetColumnIdx(table_id_name);

    for (int i = 0; i < id.size(); i++) {
      std::string temp_id = id[i];

      if (temp_id.front() == '"'
       && temp_id.back () == '"'){
        temp_id.erase(0, 1);
        temp_id.erase(temp_id.size() - 1);
      }

      if (!std::binary_search(id_DM.begin(), id_DM.end(),
                              temp_id)){
        continue;
      }
      
      idx_to_reserve.emplace_back(i);
      q_DM_DG_output_csv.SetCell(kIdColumnIdx, i, id_volumn_map[temp_id]);
    }

    std::sort(idx_to_reserve.begin(),
              idx_to_reserve.end());

    rapidcsv::Document temp_q_DM_DG_output_csv;

    for (const auto& idx : idx_to_reserve) {
      const auto& row = q_DM_DG_output_csv.GetRow<std::string>(idx);
      temp_q_DM_DG_output_csv.SetRow(temp_q_DM_DG_output_csv.GetRowCount(), row);
    }

    for (int i = 0; i < q_DM_DG_output_csv_column_name.size(); i++) {
      temp_q_DM_DG_output_csv.SetColumnName(i, q_DM_DG_output_csv_column_name[i]);
    }

    // temp_q_DM_DG_output_csv.RemoveColumn(0);
    temp_q_DM_DG_output_csv.RemoveColumn(table_id_name);

    auto temp_q_DM_DG_output_csv_colum_names 
       = temp_q_DM_DG_output_csv.GetColumnNames();

    for (const auto& colum_name : temp_q_DM_DG_output_csv_colum_names) {
      std::cout << "colum_name: " << colum_name << std::endl;
    }

    temp_q_DM_DG_output_csv.Save(dir + "/" + from_graph_csv_name + "_output" + ".csv");

    output_article_csv.RemoveColumn(table_attr_name);
    // output_article_csv.RemoveColumn(table_id_name);

    std::cout << "output_article_csv colum name: " << std::endl;
    for (const auto& colum_name : output_article_csv.GetColumnNames()) {
      std::cout << "\t" << colum_name;
    }
    std::cout  << std::endl;

    output_article_csv.Save(dir + "/" + from_table_csv_name + "_remove_" + table_attr_name + ".csv");
    
    return;
  }

  template <typename DataGraphType>
  void ConvertPokec(std::string      profiles_file,
                    std::string relationships_file,
                    std::string      output_v_file,
                    std::string      output_e_file) {

    using VertexIdType = typename GUNDAM::VertexID<DataGraphType>::type;
    using   EdgeIdType = typename GUNDAM::  EdgeID<DataGraphType>::type;

    using VertexLabelType = typename GUNDAM::VertexLabel<DataGraphType>::type;
    using   EdgeLabelType = typename GUNDAM::  EdgeLabel<DataGraphType>::type;

    const std::string kNullValue = "null";

    const std::vector<std::tuple<std::string, // name
                                 bool,  // as attr
                                 bool>> // as vertex
              kColumnNameSet = {std::tuple("user_id", false, false), // id of vertex
                                std::tuple("public", false , true),
                                std::tuple("completion_percentage", false , true),
                                std::tuple("gender", false , true),
                                std::tuple("region", false , true),
                                std::tuple("last_login", false , true),
                                std::tuple("registration", false , true),
                                std::tuple("AGE", false , true),
                                std::tuple("body", false , true),
                                std::tuple("I_am_working_in_field", false , true),
                                std::tuple("spoken_languages", false , true),
                                std::tuple("hobbies", false , true),
                                std::tuple("I_most_enjoy_good_food", false , true),
                                std::tuple("pets", false , true),
                                std::tuple("body_type", false , true),
                                std::tuple("my_eyesight", false , true),
                                std::tuple("eye_color", false , true),
                                std::tuple("hair_color", false , true),
                                std::tuple("hair_type", false , true),
                                std::tuple("completed_level_of_education", false , true),
                                std::tuple("favourite_color", false , true),
                                std::tuple("relation_to_smoking", false , true),
                                std::tuple("relation_to_alcohol", false , true),
                                std::tuple("sign_in_zodiac", false , true),
                                std::tuple("on_pokec_i_am_looking_for", false , true),
                                std::tuple("love_is_for_me", false , true),
                                std::tuple("relation_to_casual_sex", false , true),
                                std::tuple("my_partner_should_be", false , true),
                                std::tuple("marital_status", false , true),
                                std::tuple("children", false , true),
                                std::tuple("relation_to_children", false , true),
                                std::tuple("I_like_movies", false , true),
                                std::tuple("I_like_watching_movie", false , true),
                                std::tuple("I_like_music", false , true),
                                std::tuple("I_mostly_like_listening_to_music", false , true),
                                std::tuple("the_idea_of_good_evening", false , true),
                                std::tuple("I_like_specialties_from_kitchen", false , true),
                                std::tuple("fun", false , true),
                                std::tuple("I_am_going_to_concerts", false , true),
                                std::tuple("my_active_sports", false , true),
                                std::tuple("my_passive_sports", false , true),
                                std::tuple("profession", false , true),
                                std::tuple("I_like_books", false , true),
                                std::tuple("life_style", false , true),
                                std::tuple("music", false , true),
                                std::tuple("cars", false , true),
                                std::tuple("politics", false , true),
                                std::tuple("relationships", false , true),
                                std::tuple("art_culture", false , true),
                                std::tuple("hobbies_interests", false , true),
                                std::tuple("science_technologies", false , true),
                                std::tuple("computers_internet", false , true),
                                std::tuple("education", false , true),
                                std::tuple("sport", false , true),
                                std::tuple("movies", false , true),
                                std::tuple("travelling", false , true),
                                std::tuple("health", false , true),
                                std::tuple("companies_brands", false , true),
                                std::tuple("more", false , true)};

    std::map<std::string, VertexLabelType> column_name_vertex_label_map;

    for (const auto& [colum_name, as_attr, as_vertex] : kColumnNameSet) {
      assert(!as_attr || !as_vertex);
      column_name_vertex_label_map.emplace(colum_name, 
                                           column_name_vertex_label_map.size());
    }

    std::map<std::string, EdgeLabelType> edge_name_vertex_label_map;

    const std::string kRelationEdge = "has_relation";

    edge_name_vertex_label_map.emplace(kRelationEdge, edge_name_vertex_label_map.size());

    for (size_t idx = 0; idx < kColumnNameSet.size(); idx++) {
      if (!std::get<2>(kColumnNameSet[idx])) {
        // does not extract as vertex
        continue;
      }
      auto [ edge_name_vertex_label_map_it,
             edge_name_vertex_label_map_ret ]
           = edge_name_vertex_label_map.emplace(std::get<0>(kColumnNameSet[idx]), 
                                                edge_name_vertex_label_map.size());
      assert(edge_name_vertex_label_map_ret);
    }

    DataGraphType data_graph;
    
    rapidcsv::Document profiles_file_csv(
                       profiles_file, rapidcsv::LabelParams(-1, -1),
                                      rapidcsv::SeparatorParams('\t'));

    constexpr size_t vertex_id_column_idx = 0,
                        region_column_idx = 4,
                    last_login_column_idx = 5;

    std::vector<VertexIdType> vertex_id_set 
        = profiles_file_csv.GetColumn<VertexIdType>(vertex_id_column_idx);

    VertexIdType vertex_id_counter = 0;

    for (const auto& vertex_id : vertex_id_set) {
      const auto& kUserVertexLabel = column_name_vertex_label_map[std::get<0>(kColumnNameSet[0])];
      vertex_id_counter = vertex_id_counter > vertex_id?
                          vertex_id_counter : vertex_id;
      auto [vertex_handle, ret] = data_graph.AddVertex(vertex_id, kUserVertexLabel);
      assert(ret);
      // wenzhi: can be added as an optimization
      //  add attribute for each vertex here
      //  auto [attr_handle, ret] = vertex_handle->AddAttribute(attr_key, attr_value);
      //  assert(attr_handle);
      //  assert(ret);
    }
    vertex_id_counter++;

    size_t edge_id_counter = 0;

    for (size_t column_idx = 1; 
                column_idx < kColumnNameSet.size(); 
                column_idx++) {
      const std::string& kColumnName = std::get<0>(kColumnNameSet[column_idx]);
      util::Info("processing column: " + kColumnName);
      const bool kAsAttr   = std::get<1>(kColumnNameSet[column_idx]);
      const bool kAsVertex = std::get<2>(kColumnNameSet[column_idx]);
      if (!kAsAttr && !kAsVertex) {
        continue;
      }
      const std::vector<std::string> attr_set 
          = profiles_file_csv.GetColumn<std::string>(column_idx);

      std::map<std::string, VertexIdType> attr_vertex_id_map;

      assert(attr_set.size() == vertex_id_set.size());

      for (size_t row_idx = 0; row_idx < attr_set.size(); row_idx++) {
        const VertexIdType kVertexId  = vertex_id_set[row_idx];
        const std::string  kAttrValue =      attr_set[row_idx];
        if (kAttrValue == kNullValue) {
          continue;
        }
        if (kAsAttr) {
          assert(!kAsVertex);
          auto vertex_handle = data_graph.FindVertex(kVertexId);
          assert(vertex_handle);
          auto [attr_handle, ret] = vertex_handle->AddAttribute(kColumnName, kAttrValue);
          assert(attr_handle);
          assert(ret);
          continue;
        }
        assert(kAsVertex);
        auto attr_vertex_id_map_it 
           = attr_vertex_id_map.find(kAttrValue);
        if ( attr_vertex_id_map_it
          == attr_vertex_id_map.end() ) {
          // does not contain this value before
          attr_vertex_id_map_it = attr_vertex_id_map.emplace_hint(attr_vertex_id_map_it,
                                                                 kAttrValue,
                                                                  vertex_id_counter++);
          assert(column_name_vertex_label_map.find(kColumnName)
              != column_name_vertex_label_map.end());
          auto [vertex_handle, 
                vertex_ret] = data_graph.AddVertex(attr_vertex_id_map_it->second, 
                                            column_name_vertex_label_map[kColumnName]);
          assert(vertex_handle);
          assert(vertex_ret);
          auto [attr_handle, 
                attr_ret] = vertex_handle->AddAttribute("attr", kAttrValue);
          assert(attr_handle);
          assert(attr_ret);
        }
        assert(attr_vertex_id_map_it 
            != attr_vertex_id_map.end());
        assert(data_graph.FindVertex(attr_vertex_id_map_it->second));
        
        auto [edge_handle, ret] = data_graph.AddEdge(kVertexId, 
                                                     attr_vertex_id_map_it->second, 
                                                     edge_name_vertex_label_map[kColumnName],
                                                     edge_id_counter++);
        assert(edge_handle);
        assert(ret);
      }
    }
    vertex_id_set.clear();

    rapidcsv::Document relationships_file_csv(
                       relationships_file, rapidcsv::LabelParams(-1, -1),
                                           rapidcsv::SeparatorParams('\t'));

    constexpr size_t src_id_column_idx = 0,
                     dst_id_column_idx = 0;

    const std::vector<VertexIdType> src_id_set 
        = relationships_file_csv.GetColumn<VertexIdType>(src_id_column_idx);

    const std::vector<VertexIdType> dst_id_set 
        = relationships_file_csv.GetColumn<VertexIdType>(dst_id_column_idx);

    assert(src_id_set.size() == dst_id_set.size());

    for (size_t row = 0; row < src_id_set.size(); row++) {
      assert(column_name_vertex_label_map.find(kRelationEdge)
          != column_name_vertex_label_map.end());
      const auto& kRelationEdgeLabel = column_name_vertex_label_map[kRelationEdge];
      auto [edge_handle, ret] = data_graph.AddEdge(src_id_set[row], 
                                                   dst_id_set[row], 
                                                  kRelationEdgeLabel,
                                                   edge_id_counter++);
      assert(ret);
    }

    GUNDAM::WriteCSVGraph(data_graph, output_v_file,
                                      output_e_file);

    return;
  }

  /**
   * @brief Partial evaluation for SSSP.
   *
   * @param frag
   * @param ctx
   * @param messages
   */
  void PEval(const fragment_t& frag, context_t& ctx,
             message_manager_t& messages) {
    auto inner_vertices = frag.InnerVertices();
    messages.InitChannels(thread_num());

#ifdef PROFILING
    ctx.exec_time -= GetCurrentTime();
#endif
    FragmentGraphWithIndex<const fragment_t> index_graph(frag);
    // index_graph.BuildIndex();

    using VertexIDType = int;
    using VertexLabelType = int;
    using VertexAttributeKeyType = std::string;
    using   EdgeIDType = int;
    using   EdgeLabelType = int;
    using   EdgeAttributeKeyType = std::string;

    using DataGraph =
        GUNDAM::LargeGraph<VertexIDType, VertexLabelType, std::string, 
                             EdgeIDType,   EdgeLabelType, std::string>;
    using Pattern =
        GUNDAM::LargeGraph<VertexIDType, VertexLabelType, std::string, 
                             EdgeIDType,   EdgeLabelType, std::string>;
   
    using GarType = gar::GraphAssociationRule<Pattern, DataGraph>;

    // std::vector<std::string> pattern_path_set = {
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_0",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_1",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_2",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_3",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_4",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_5",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_6",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_7",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_8",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_9",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_10",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_11",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_12",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_13",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_14",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_15",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_16",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_17",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_18",
    //   "../../Experiments/gar/imdb/gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/gar_level_8_19"};

    // for (const auto& pattern_path : pattern_path_set) {
    //   Pattern pattern;
    //   GUNDAM::ReadCSVGraph(pattern, pattern_path + "_v.csv",
    //                                 pattern_path + "_e.csv");
    //   for (auto vertex_it = pattern.VertexBegin();
    //            !vertex_it.IsDone();
    //             vertex_it++) {
    //     auto pattern_tree = GUNDAM::ToTree(pattern, vertex_it);
    //     GUNDAM::WriteCSVGraph(pattern_tree, pattern_path + "_rooted_" + std::to_string(vertex_it->id()) + "_v.csv",
    //                                         pattern_path + "_rooted_" + std::to_string(vertex_it->id()) + "_e.csv");
    //   }
    // }

    // using DataGraph = GUNDAM::Graph<
    //   GUNDAM::SetVertexIDType<VertexIDType>,
    //   GUNDAM::SetVertexAttributeStoreType<GUNDAM::AttributeType::kGrouped>,
    //   // GUNDAM::SetVertexAttributeStoreType<GUNDAM::AttributeType::kSeparated>,
    //   GUNDAM::SetVertexLabelType<VertexLabelType>,
    //   GUNDAM::SetVertexLabelContainerType<GUNDAM::ContainerType::Map>,
    //   GUNDAM::SetVertexIDContainerType<GUNDAM::ContainerType::Map>,
    //   GUNDAM::SetVertexPtrContainerType<GUNDAM::ContainerType::Map>,
    //   GUNDAM::SetEdgeLabelContainerType<GUNDAM::ContainerType::Map>,
    //   GUNDAM::SetVertexAttributeKeyType<VertexAttributeKeyType>,
    //   GUNDAM::SetEdgeIDType<EdgeIDType>,
    //   GUNDAM::SetEdgeAttributeStoreType<GUNDAM::AttributeType::kGrouped>,
    //   // GUNDAM::SetEdgeAttributeStoreType<GUNDAM::AttributeType::kSeparated>,
    //   GUNDAM::SetEdgeLabelType<EdgeLabelType>,
    //   GUNDAM::SetEdgeAttributeKeyType<EdgeAttributeKeyType>>;

    // ConvertPokec<DataGraph>("../dataset/pokec/soc-pokec-profiles.txt",
    //                         "../dataset/pokec/soc-pokec-relationships.txt",
    //                         "../dataset/pokec/pokec_12_11_v.csv",
    //                         "../dataset/pokec/pokec_12_11_e.csv");

    // RoundAttrToHalf<DataGraph>("../dataset/reclogic_expr_data/movielen/train_graph/movielen_v.csv",
    //                            "../dataset/reclogic_expr_data/movielen/train_graph/movielen_e.csv",
    //                            "avgrating",
    //                            "../dataset/reclogic_expr_data/movielen/train_graph/movielen_avgratin_to_half_v.csv",
    //                            "../dataset/reclogic_expr_data/movielen/train_graph/movielen_avgratin_to_half_e.csv");


    using CandidateSetType = std::map<typename GUNDAM::VertexHandle<DataGraph>::type, 
                          std::vector<typename GUNDAM::VertexHandle<DataGraph>::type>>;

    // std::vector<size_t> partition_size_set = {1000, 10000, 100000};

    std::vector<
    std::pair<std::string, std::string>> 
      data_graph_path_set = {std::pair("../dataset/yago/yago",      "yago"),
                             std::pair("../dataset/dblp_1980/dblp", "dblp_1980")};


    std::vector<size_t> partition_size_set = {1000000};

    // std::vector<
    // std::pair<std::string, std::string>> 
    //   data_graph_path_set = {std::pair("../dataset/dblp_1980/dblp", "dblp_1980")};

    for (const auto& [data_graph_path, 
                      data_graph_name] : data_graph_path_set) {
      DataGraph data_graph;
      GUNDAM::ReadCSVGraph(data_graph, data_graph_path + "_v.csv",
                                       data_graph_path + "_e.csv");
      GUNDAM::RemoveIsolateVertex(data_graph);

      std::cout << "data_graph.CountVertex() after remove iso vertex: "
                <<  data_graph.CountVertex() << std::endl;

      auto cc_set = GUNDAM::ConnectedComponent(data_graph);

      std::sort(cc_set.begin(), cc_set.end(),
            [](const auto& a, const auto& b){
                return a.CountVertex() < b.CountVertex();});

      std::cout << "cc_set.size(): "
                <<  cc_set.size() << std::endl;

      const size_t kCcSizeBound = 100;

      for (const auto& cc : cc_set) {
        if (cc.CountVertex() > kCcSizeBound) {
          continue;
        }
        for (auto vertex_it = cc.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          data_graph.EraseVertex(vertex_it->id());
        }
      }

      std::cout << "data_graph.CountVertex() after remove small cc: "
                <<  data_graph.CountVertex() << std::endl;


      using EquivalentClassType = std::vector<
                                  std::vector<typename GUNDAM::VertexHandle<DataGraph>::type>>;

      EquivalentClassType equivalent_class_set;

      GUNDAM::BisimulationGeneralCase(data_graph, equivalent_class_set);

      std::ofstream equivalent_class_file("./" + data_graph_name + "_bisimulation.equivalent_class");
      for (const auto& equivalent_class : equivalent_class_set) {
        for (const auto& vertex_handle : equivalent_class) {
          equivalent_class_file << "\t" << vertex_handle->id();
        }
        equivalent_class_file << std::endl;
      }

      // const VertexLabelType kVertexLabel = 2; // paper
      // const size_t kDuplicateVertexSize = 50000; 

      // std::vector<typename GUNDAM::VertexHandle<DataGraph>::type> vertex_handle_set;
      // for (auto vertex_it = data_graph.VertexBegin();
      //          !vertex_it.IsDone();
      //           vertex_it++) {
      //   if (vertex_it->label() != kVertexLabel) {
      //     continue;
      //   }
      //   vertex_handle_set.emplace_back(vertex_it);
      // }
      // if (vertex_handle_set.size() < kDuplicateVertexSize) {
      //   util::Error("does not have enough vertex to duplicate");
      //   return;
      // }

      // std::cout << "vertex_handle_set.size(): "
      //           <<  vertex_handle_set.size() << std::endl;

      // std::random_shuffle ( vertex_handle_set.begin(), 
      //                       vertex_handle_set.end() );

      // vertex_handle_set.resize(kDuplicateVertexSize);

      // auto vertex_id_map = GUNDAM::DuplicateVertex(data_graph, vertex_handle_set);

      // for (const auto& partition_size : partition_size_set) {
      //   CandidateSetType candidate_set;
      //   for (auto vertex_it = data_graph.VertexBegin();
      //            !vertex_it.IsDone();
      //             vertex_it++) {
      //     candidate_set.emplace(vertex_it, std::vector<typename GUNDAM::VertexHandle<DataGraph>::type>());
      //   }

      //   auto begin_time = std::chrono::high_resolution_clock::now();
      //   if (!RefineCandidateSetPartByPart<GUNDAM::execution::parallel_policy>(
      //                   data_graph,   
      //                   data_graph,  
      //             candidate_set, partition_size, "./" + data_graph_name + "_simulation_parallel_part_by_part_" + std::to_string(partition_size) + ".log2")) {
      //     return;
      //   }
      //   std::cout << "##  Time with partition size " << partition_size << ": "
      //             << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin_time).count()
      //             << " (ms)" << std::endl;

      //   std::ofstream ofile("./" + data_graph_name + "_simulation_parallel_part_by_part_" + std::to_string(partition_size) + ".out2");
      //   for (const auto& [query_handle, candidate] : candidate_set) {
      //     ofile << query_handle->id();
      //     for (const auto& target_handle : candidate) {
      //       ofile << "," << target_handle->id();
      //     }
      //     ofile << std::endl;
      //   }

      //   std::ofstream recall_duplicate_file("./" + data_graph_name + "_simulation_parallel_part_by_part_" + std::to_string(partition_size) + ".recall_duplicate_out");
      //   size_t true_positive_counter = 0;
      //   for (const auto& [ori_vertex_handle, 
      //                     new_vertex_handle] : vertex_id_map) {
      //     std::sort(candidate_set[ori_vertex_handle].begin(),
      //               candidate_set[ori_vertex_handle].end());
      //     std::sort(candidate_set[new_vertex_handle].begin(),
      //               candidate_set[new_vertex_handle].end());

      //     if ((!std::binary_search(candidate_set[ori_vertex_handle].begin(),
      //                              candidate_set[ori_vertex_handle].end(),
      //                                            new_vertex_handle))
      //       &&(!std::binary_search(candidate_set[new_vertex_handle].begin(),
      //                              candidate_set[new_vertex_handle].end(),
      //                                            ori_vertex_handle))) {
      //       continue;
      //     }
      //     true_positive_counter++;
      //     recall_duplicate_file << ori_vertex_handle->id() << ", "
      //                           << new_vertex_handle->id() << std::endl;
      //   }

      //   std::sort(vertex_id_map.begin(), 
      //             vertex_id_map.end());

      //   // construct a negative test size
      //   std::set<std::pair<VertexIDType,
      //                      VertexIDType>> duplciate_set;
      //   std::vector<std::pair<typename GUNDAM::VertexHandle<DataGraph>::type,
      //                         typename GUNDAM::VertexHandle<DataGraph>::type>> negative_set;
      //   for (const auto& [query_handle, candidate] : candidate_set) {
      //     if (query_handle->label() != kVertexLabel) {
      //       continue;
      //     }
      //     for (const auto& target_handle : candidate) {
      //       if (target_handle->label() != kVertexLabel) {
      //         continue;
      //       }
      //       auto [duplciate_set_it,
      //             duplciate_set_ret] 
      //           = duplciate_set.emplace(std::min(query_handle->id(), target_handle->id()),
      //                                   std::max(query_handle->id(), target_handle->id()));
      //       if (!duplciate_set_ret) {
      //         // already existed
      //         continue;
      //       }
      //       if (std::binary_search(vertex_id_map.begin(), 
      //                              vertex_id_map.end(),
      //                             (query_handle->id() < target_handle->id() ?
      //                          std::pair( query_handle, target_handle)
      //                        : std::pair(target_handle,  query_handle)))) {
      //         // is in positive set
      //         continue;
      //       }
      //       negative_set.emplace_back(query_handle,
      //                                target_handle);
      //     }
      //   }
      //   std::random_shuffle ( negative_set.begin(), 
      //                         negative_set.end() );
      //   negative_set.resize( vertex_id_map.size() );

      //   size_t true_negative_counter = 0;
      //   for (const auto& [ori_vertex_handle, 
      //                     new_vertex_handle] : negative_set) {
      //     if (( std::binary_search(candidate_set[ori_vertex_handle].begin(),
      //                              candidate_set[ori_vertex_handle].end(),
      //                                            new_vertex_handle))
      //       ||( std::binary_search(candidate_set[new_vertex_handle].begin(),
      //                              candidate_set[new_vertex_handle].end(),
      //                                            ori_vertex_handle))) {
      //       continue;
      //     }
      //     true_negative_counter++;
      //   }

      //   std::cout << "#######################################"  << std::endl
      //             << "#       positive   size : " << vertex_id_map.size()  << std::endl
      //             << "#  true positive counter: " << true_positive_counter << std::endl
      //             << "#       negative   size : " << vertex_id_map.size()  << std::endl
      //             << "#  true negative counter: " << true_negative_counter << std::endl
      //             << "#######################################"  << std::endl;

      //   std::ofstream all_duplicate_file("./" + data_graph_name + "_simulation_parallel_part_by_part_" + std::to_string(partition_size) + ".all_duplicate_out");
      //   duplciate_set.clear();
      //   for (const auto& [query_handle, candidate] : candidate_set) {
      //     if (query_handle->label() != kVertexLabel) {
      //       continue;
      //     }
      //     for (const auto& target_handle : candidate) {
      //       if (target_handle->label() != kVertexLabel) {
      //         continue;
      //       }
      //       auto [duplciate_set_it,
      //             duplciate_set_ret] 
      //           = duplciate_set.emplace(std::min(query_handle->id(), target_handle->id()),
      //                                   std::max(query_handle->id(), target_handle->id()));
      //       if (!duplciate_set_ret) {
      //         // already existed
      //         continue;
      //       }
      //       all_duplicate_file <<  query_handle->id() << ", "
      //                          << target_handle->id() << std::endl;
      //     }
      //   }
      //   candidate_set.clear();
      // }
    }

    // GUNDAM::ReadCSVGraph(data_graph, "../dataset/imdb/imdb_v.csv",
    //                                  "../dataset/imdb/imdb_e.csv");

    // GUNDAM::RemoveIsolateVertex(data_graph);

    // std::cout << "data_graph.CountVertex(): "
    //           <<  data_graph.CountVertex() << std::endl;
          
    // auto cc_set = GUNDAM::ConnectedComponent(data_graph);

    // std::sort(cc_set.begin(), cc_set.end(), 
    //        [](const auto& a, const auto& b){
    //                return a.CountVertex() 
    //                     < b.CountVertex();});

    // std::cout << "cc_set.size(): "
    //           <<  cc_set.size() << std::endl;

    // const size_t kCcSizeBound = 100;

    // for (const auto& cc : cc_set) {
    //   std::cout << "cc.CountVertex(): "
    //             <<  cc.CountVertex() << std::endl;
    //   if (cc.CountVertex() > kCcSizeBound) {
    //     continue;
    //   }
    //   for (auto vertex_it = cc.VertexBegin();
    //            !vertex_it.IsDone();
    //             vertex_it++) {
    //     data_graph.EraseVertex(vertex_it->id());
    //   }
    // }

    // getchar();

    // for (const auto& partition_size : partition_size_set) {
    //   CandidateSetType candidate_set;
    //   for (auto vertex_it = data_graph.VertexBegin();
    //           !vertex_it.IsDone();
    //             vertex_it++) {
    //     candidate_set.emplace(vertex_it, std::vector<typename GUNDAM::VertexHandle<DataGraph>::type>());
    //   }

    //   auto begin_time = std::chrono::high_resolution_clock::now();
    //   if (!RefineCandidateSetPartByPart<GUNDAM::execution::parallel_policy>(
    //                   data_graph,   
    //                   data_graph,  
    //              candidate_set, partition_size, "./imdb_simulation_parallel_part_by_part_" + std::to_string(partition_size) + ".log2")) {
    //     return;
    //   }
    //   std::cout << "##  Time with partition size " << partition_size << ": "
    //             << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin_time).count()
    //             << " (ms)" << std::endl;

    //   std::ofstream ofile("./imdb_simulation_parallel_part_by_part_" + std::to_string(partition_size) + ".out2");
    //   for (const auto& [query_handle, candidate] : candidate_set) {
    //     ofile << query_handle->id();
    //     for (const auto& target_handle : candidate) {
    //       ofile << "," << target_handle->id();
    //     }
    //     ofile << std::endl;
    //   }

    //   candidate_set.clear();
    // }

    // std::vector<std::string> gar_pattern_str_set 
    //   // = {"./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1/gar_level_0_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1/gar_level_1_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1/gar_level_2_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1/gar_level_3_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1/gar_level_4_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1/gar_level_5_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1/gar_level_6_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1/gar_level_7_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1/gar_level_8_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1/gar_level_9_worker_0"};
    //   // = {"./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.1_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_0_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.1_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_1_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.1_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_2_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.1_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_3_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.1_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_4_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.1_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_5_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.1_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_6_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.1_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_7_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.1_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_8_worker_0",
    //   //    "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.1_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_9_worker_0"};
    //   = {"./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.05_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_0_worker_0",
    //      "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.05_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_1_worker_0",
    //      "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.05_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_2_worker_0",
    //      "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.05_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_3_worker_0",
    //      "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.05_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_4_worker_0",
    //      "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.05_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_5_worker_0",
    //      "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.05_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_6_worker_0",
    //      "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.05_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_7_worker_0",
    //      "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.05_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_8_worker_0",
    //      "./gcr_movielens_data_movielen_10_1_s1_e1_d0_s1_e0_d0_supp_1_conf_0.05_single_direct_path_9a5fdb2b2daa07f34cb6849fbf46fc3df8c048d5/gar_level_9_worker_0"};

    // std::vector<Pattern> gar_pattern_set;
    // for (const auto& gar_pattern_str : gar_pattern_str_set) {
    //   std::vector<Pattern> temp_gar_pattern_set;
    //   GUNDAM::ReadCSVGraphSet(temp_gar_pattern_set, gar_pattern_str + "_v_set.csv",
    //                                                 gar_pattern_str + "_e_set.csv");
    //   for (auto& temp_gar_pattern : temp_gar_pattern_set) {
        
    //     std::vector<Pattern> 
    //           connected_components(GUNDAM::ConnectedComponent<true>(temp_gar_pattern));
    //     bool is_one_direct_path = true;

    //     for (const auto& cc : connected_components) {
    //       if (GUNDAM::IsStar<false>(cc)) {
    //         continue;
    //       }
    //       is_one_direct_path = false;
    //       break;
    //     }

    //     if (!is_one_direct_path) {
    //       continue;
    //     }
    //     gar_pattern_set.emplace_back(std::move(temp_gar_pattern));
    //   }
    // }

    // std::cout << "gar_pattern_set.size(): " 
    //           <<  gar_pattern_set.size() << std::endl;

    // GUNDAM::DeduplicatePatternsWithDfsCode(gar_pattern_set);

    // std::cout << "deduplicated gar_pattern_set.size(): " 
    //                         << gar_pattern_set.size() << std::endl;

    // std::vector<Pattern> graph_set;
    // GUNDAM::ReadCSVGraphSet(graph_set, "./TIE_big/TIE_v.cvs",
    //                                    "./TIE_big/TIE_e.cvs");

    // GUNDAM::DeduplicatePatternsWithDfsCode(graph_set);

    // std::cout << "deduplicated TIE_big.size(): " 
    //                       << graph_set.size() << std::endl;

    // std::sort(graph_set.begin(),
    //           graph_set.end(), [](const Pattern& a,
    //                               const Pattern& b){ return a.CountEdge() < b.CountEdge();});

    // GUNDAM::WriteCSVGraphSet(graph_set, "./TIE_big/TIE_uniqued_v.cvs",
    //                                     "./TIE_big/TIE_uniqued_e.cvs");

    // std::vector<Pattern> gar_pattern_contained_in_graph_set;
    // for (const auto& gar_pattern : gar_pattern_set) {
    //   for (const auto& graph : graph_set) {
    //     if (graph.CountEdge() > gar_pattern.CountEdge()) {
    //       break;
    //     }
    //     if (!GUNDAM::SamePattern(gar_pattern, graph)) {
    //       continue;
    //     }
    //     gar_pattern_contained_in_graph_set.emplace_back(graph);
    //     break;
    //   }
    // }
    
    // size_t gar_pattern_set_counter = 0;
    // for (const auto& gar_pattern : gar_pattern_set) {
    //   for (const auto& gar_pattern_contained : gar_pattern_contained_in_graph_set) {
    //     if (!GUNDAM::SubGraphOf(gar_pattern, 
    //                             gar_pattern_contained)) {
    //       continue;
    //     }
    //     gar_pattern_set_counter++;
    //     break;
    //   }
    // }

    // std::cout << "gar_pattern_set covered in TIE_big: "
    //           << gar_pattern_set_counter
    //           << std::endl;

    // graph_set.clear();
    // GUNDAM::ReadCSVGraphSet(graph_set, "./TIE_Small/TIE_Small_v.cvs",
    //                                    "./TIE_Small/TIE_Small_e.cvs");

    // GUNDAM::DeduplicatePatternsWithDfsCode(graph_set);

    // std::cout << "deduplicated TIE_Small.size(): " 
    //                         << graph_set.size() << std::endl;

    // std::sort(graph_set.begin(),
    //           graph_set.end(), [](const Pattern& a,
    //                               const Pattern& b){ return a.CountEdge() < b.CountEdge();});

    // GUNDAM::WriteCSVGraphSet(graph_set, "./TIE_Small/TIE_uniqued_v.cvs",
    //                                     "./TIE_Small/TIE_uniqued_e.cvs");


    // gar_pattern_contained_in_graph_set.clear();
    // for (const auto& gar_pattern : gar_pattern_set) {
    //   for (const auto& graph : graph_set) {
    //     if (graph.CountEdge() > gar_pattern.CountEdge()) {
    //       break;
    //     }
    //     if (!GUNDAM::SamePattern(gar_pattern, graph)) {
    //       continue;
    //     }
    //     gar_pattern_contained_in_graph_set.emplace_back(graph);
    //     break;
    //   }
    // }
    
    // gar_pattern_set_counter = 0;
    // for (const auto& gar_pattern : gar_pattern_set) {
    //   for (const auto& gar_pattern_contained : gar_pattern_contained_in_graph_set) {
    //     if (!GUNDAM::SubGraphOf(gar_pattern, 
    //                             gar_pattern_contained)) {
    //       continue;
    //     }
    //     gar_pattern_set_counter++;
    //     break;
    //   }
    // }
    
    // std::cout << "gar_pattern_set covered in TIE_Small: "
    //           << gar_pattern_set_counter
    //           << std::endl;

    // RoundAttrToHalf<DataGraph>("../dataset/imdb/folded_imdb_v.csv",
    //                            "../dataset/imdb/folded_imdb_e.csv",
    //                            "rating_is_rating",
    //                            "../dataset/imdb/folded_imdb_rating_to_half_v.csv",
    //                            "../dataset/imdb/folded_imdb_rating_to_half_e.csv");

    // DataGraph data_graph;
    // GUNDAM::ReadCSVGraph(data_graph, "../dataset/dbpedia/dbpedia_v.csv",
    //                                  "../dataset/dbpedia/dbpedia_e.csv");

    // std::ofstream bfs_seq("../dataset/dbpedia/bfs_seq.txt");

    // double relative_vertex_preserve_ratio = 0.1;
     
    // std::vector<typename GUNDAM::VertexHandle<DataGraph>::type> vertex_handle_seq;

    // std::set<typename GUNDAM::VertexHandle<DataGraph>::type> visited_vertex;

    // auto store_vertex_id_set_callback 
    //       = [&vertex_handle_seq,
    //          &visited_vertex](const typename GUNDAM::VertexHandle<DataGraph>::type& vertex_handle) {
    //   vertex_handle_seq.emplace_back(vertex_handle);
    //   visited_vertex.emplace(vertex_handle);
    //   return true; // continue bfs
    // };

    // for (auto vertex_it = data_graph.VertexBegin();
    //          !vertex_it.IsDone();
    //           vertex_it++) {
    //   if (visited_vertex.find(vertex_it)
    //    != visited_vertex.end()) {
    //     continue;
    //   }
    //   GUNDAM::Bfs<true>(data_graph, vertex_it,
    //                     store_vertex_id_set_callback);
    //   assert(vertex_handle_seq.size() == visited_vertex.size());
    //   if (vertex_handle_seq.size() == data_graph.CountVertex()) {
    //     break;
    //   }
    // }

    // assert(vertex_handle_seq.size() == data_graph.CountVertex());

    // std::set<EdgeLabelType> dbpedia_relative_edge_label_set = {
    //       48, // kingdom
    //       49, // phylum
    //       50, // class
    //       51, // order
    //       52};// family

    // for (const auto& vertex_handle : vertex_handle_seq) {
    //   bool is_relative_vertex = false;
    //   for (auto out_edge_it = vertex_handle->OutEdgeBegin();
    //            !out_edge_it.IsDone();
    //             out_edge_it++) {
    //     if (dbpedia_relative_edge_label_set.find(out_edge_it->label())
    //      == dbpedia_relative_edge_label_set.end()) {
    //       // is not relative
    //       continue;
    //     }
    //     // is relative
    //     is_relative_vertex = true;
    //     break;
    //   }
    //   if (!is_relative_vertex){
    //     for (auto in_edge_it = vertex_handle->InEdgeBegin();
    //              !in_edge_it.IsDone();
    //               in_edge_it++) {
    //       if (dbpedia_relative_edge_label_set.find(in_edge_it->label())
    //        == dbpedia_relative_edge_label_set.end()) {
    //         // is not relative
    //         continue;
    //       }
    //       // is relative
    //       is_relative_vertex = true;
    //       break;
    //     }
    //   }
    //   if (is_relative_vertex) {
    //     const double kR = ((double) std::rand() / (RAND_MAX));
    //     if (kR > relative_vertex_preserve_ratio) {
    //       // does not need to preserve
    //       continue;
    //     }
    //   }
    //   bfs_seq << vertex_handle->id() << std::endl;
    // }
    
    // bfs_seq.close();
    
    // DataGraph data_graph;
    // GUNDAM::ReadCSVGraph(data_graph, "../graph_dataset/imdb/imdb_v.csv",
    //                                  "../graph_dataset/imdb/imdb_e.csv");

    // std::ofstream bfs_seq("../graph_dataset/imdb/bfs_seq.txt");

    // // double relative_vertex_preserve_ratio = 0.1;
     
    // std::vector<typename GUNDAM::VertexHandle<DataGraph>::type> vertex_handle_seq;

    // std::set<typename GUNDAM::VertexHandle<DataGraph>::type> visited_vertex;

    // auto store_vertex_id_set_callback 
    //       = [&vertex_handle_seq,
    //          &visited_vertex](const typename GUNDAM::VertexHandle<DataGraph>::type& vertex_handle) {
    //   vertex_handle_seq.emplace_back(vertex_handle);
    //   visited_vertex.emplace(vertex_handle);
    //   return true; // continue bfs
    // };

    // for (auto vertex_it = data_graph.VertexBegin();
    //          !vertex_it.IsDone();
    //           vertex_it++) {
    //   if (visited_vertex.find(vertex_it)
    //    != visited_vertex.end()) {
    //     continue;
    //   }
    //   GUNDAM::Bfs<true>(data_graph, vertex_it,
    //                     store_vertex_id_set_callback);
    //   assert(vertex_handle_seq.size() == visited_vertex.size());
    //   if (vertex_handle_seq.size() == data_graph.CountVertex()) {
    //     break;
    //   }
    // }

    // assert(vertex_handle_seq.size() == data_graph.CountVertex());

    // std::set<EdgeLabelType> dbpedia_relative_edge_label_set = {
    //       48, // kingdom
    //       49, // phylum
    //       50, // class
    //       51, // order
    //       52};// family

    // for (const auto& vertex_handle : vertex_handle_seq) {
    //   // bool is_relative_vertex = false;
    //   // for (auto out_edge_it = vertex_handle->OutEdgeBegin();
    //   //          !out_edge_it.IsDone();
    //   //           out_edge_it++) {
    //   //   if (dbpedia_relative_edge_label_set.find(out_edge_it->label())
    //   //    == dbpedia_relative_edge_label_set.end()) {
    //   //     // is not relative
    //   //     continue;
    //   //   }
    //   //   // is relative
    //   //   is_relative_vertex = true;
    //   //   break;
    //   // }
    //   // if (!is_relative_vertex){
    //   //   for (auto in_edge_it = vertex_handle->InEdgeBegin();
    //   //            !in_edge_it.IsDone();
    //   //             in_edge_it++) {
    //   //     if (dbpedia_relative_edge_label_set.find(in_edge_it->label())
    //   //      == dbpedia_relative_edge_label_set.end()) {
    //   //       // is not relative
    //   //       continue;
    //   //     }
    //   //     // is relative
    //   //     is_relative_vertex = true;
    //   //     break;
    //   //   }
    //   // }
    //   // if (is_relative_vertex) {
    //   //   const double kR = ((double) std::rand() / (RAND_MAX));
    //   //   if (kR > relative_vertex_preserve_ratio) {
    //   //     // does not need to preserve
    //   //     continue;
    //   //   }
    //   // }
    //   bfs_seq << vertex_handle->id() << std::endl;
    // }
    
    // bfs_seq.close();

    // GetData("/home/fwz/project/DMDG/DBpedia/",
    //               "q1_DM_DG",    "id_D",   "birth_place",
    //             "Politician", "\"URI\"", "\"birthPlace\"");

    // GetData("/home/fwz/project/DMDG/DBpedia/",
    //               "q2_DM_DG", "id_D", "birthDate",
    //                "Athlete", "\"URI\"" , "\"birthDate\"");

    // GetData("/home/fwz/project/DMDG/DBLP/",
    //               "q3_DM_DG", "id_D", "volumn",
    //         "output_article", "id"  , "volume");

    // GetData("/home/fwz/project/DMDG/DBLP/",
    //               "q4_DM_DG", "id_D", "affiliation",
    //       "output_phdthesis", "id"  , "school");

    // GetData("/home/fwz/project/DMDG/IMDB/",
    //               "q5_DM_DG",   "id_D", "writers",
    //             "title.crew", "tconst", "writers");

    // GetData("/home/fwz/project/DMDG/IMDB/",
    //               "q6_DM_DG",   "id_D", "primaryProfession",
    //            "name.basics", "nconst", "nconst");
    
    // GetData("/home/fwz/JedAIToolkit/data/cleanCleanErDatasets/DMDG/DBpedia/",
    //               "q1_DM_DG",    "id_D",   "birth_place",
    //             "Politician_slim", "URI",  "birthPlace");

    // GetData("/home/fwz/JedAIToolkit/data/cleanCleanErDatasets/DMDG/DBpedia/",
    //               "q2_DM_DG", "id_D", "birthDate",
    //                 "Athlete_slim", "URI" , "birthDate");

    // GetData("/home/fwz/JedAIToolkit/data/cleanCleanErDatasets/DMDG/DBLP/",
    //               "q3_DM_DG", "id_D", "volumn",
    //         "output_article_slim", "id"  , "volume");

    // GetData("/home/fwz/JedAIToolkit/data/cleanCleanErDatasets/DMDG/DBLP/",
    //               "q4_DM_DG", "id_D", "affiliation",
    //       "output_phdthesis_slim", "id"  , "school");

    // GetData("/home/fwz/JedAIToolkit/data/cleanCleanErDatasets/DMDG/IMDB/",
    //               "q5_DM_DG",   "id_D", "writers",
    //             "title.crew.csv_slim", "tconst", "writers");

    // GetData("/home/fwz/JedAIToolkit/data/cleanCleanErDatasets/DMDG/IMDB/",
    //               "q6_DM_DG",   "id_D", "primaryProfession",
    //             "name_basics_csv_slim", "nconst", "nconst");

    // DataGraph dbpedia;

    // GUNDAM::ReadCSVGraph(dbpedia, "../dataset/dbpedia/dbpedia_v.csv",
    //                               "../dataset/dbpedia/dbpedia_e.csv");

    // GUNDAM::WriteCSVGraph(dbpedia, "../dataset/dbpedia/dbpedia_reordered_v.csv",
    //                                "../dataset/dbpedia/dbpedia_reordered_e.csv");

    // DataGraph yago;

    // GUNDAM::ReadCSVGraph(dbpedia, "../dataset/yago/yago_v.csv",
    //                               "../dataset/yago/yago_e.csv");

    // GUNDAM::WriteCSVGraph(dbpedia, "../dataset/yago/yago_reordered_v.csv",
    //                                "../dataset/yago/yago_reordered_e.csv");

    // DataGraph imdb;

    // GUNDAM::ReadCSVGraph(dbpedia, "../dataset/imdb/imdb_v.csv",
    //                               "../dataset/imdb/imdb_e.csv");

    // GUNDAM::WriteCSVGraph(dbpedia, "../dataset/imdb/imdb_reordered_v.csv",
    //                                "../dataset/imdb/imdb_reordered_e.csv");

    // DataGraph dblp;

    // GUNDAM::ReadCSVGraph(dbpedia, "../dataset/dblp_1980/dblp_v.csv",
    //                               "../dataset/dblp_1980/dblp_e.csv");

    // GUNDAM::WriteCSVGraph(dbpedia, "../dataset/dblp_1980/dblp_reordered_v.csv",
    //                                "../dataset/dblp_1980/dblp_reordered_e.csv");

    // using GarType = gar::GraphAssociationRule<Pattern, DataGraph>;
    // // GarType gar1, gar2;
    // std::vector<GarType> gar_list;
    // // gar_list.push_back(gar1);
    // // gar_list.push_back(gar2);

    // // std::string gar_dir = "./gar_dblp_1980_all_graph_3_1_supp_100/";

    // // std::string gar_dir = "./gar_dblp_1980_locality_s2_e1001_d5_random_walk_d2_s30_0.1_all_graph_15_1_v6_supp_10_app_driven_with_full_restrict_7_edge_type_8_nodes_no_time_limit/";

    // // std::string gar_dir = "./gar_imdb_locality_in_each_cluster_s35_e1_d34_s37_e1_d34_cluster_0_bfs_width_limit_r2_w3_0.100000_10_1_v6_triple_larger_10000_without_edge_literal_2_hop_connected_supp_700_app_driven_7_edge_type/";

    // // std::string gar_dir = "./gar_yago_locality_contain_all_small_cluster_s0_e0_d0_random_walk_d2_s30_0.1_10_1_v6_with_constant_edge_type_larger_10000_supp_200_with_full_restrict_gar_set_app_driven_7_edge_type/";

    // // std::string gar_dir = "./gar_dblp_1980_locality_contain_all_small_cluster_s2_e1006_d10_cluster_0_bfs_width_limit_r2_w3_0.1_supp_10_8_nodes_no_time_limit/";

    // std::string gar_dir = ctx.yaml_file_;

    // gar::ReadGARSet(gar_list, "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/gar_level_2_worker_0_v_set.csv",
    //                           "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/gar_level_2_worker_0_e_set.csv",
    //                           "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/gar_level_2_worker_0_x_set.csv",
    //                           "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/gar_level_2_worker_0_y_set.csv");
    
    // for (int gar_idx = 0; gar_idx < gar_list.size(); gar_idx++){
    //   gar::WriteGAR(gar_list[gar_idx],
    //              "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/_gar_level_2_worker_" + std::to_string(gar_idx) + "_v.csv",
    //              "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/_gar_level_2_worker_" + std::to_string(gar_idx) + "_e.csv",
    //              "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/_gar_level_2_worker_" + std::to_string(gar_idx) + "_x.csv",
    //              "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/_gar_level_2_worker_" + std::to_string(gar_idx) + "_y.csv");
    // }

    // gar_list.clear();
    // gar::ReadGARSet(gar_list, "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt.bk/gar_level_2_worker_0_v_set.csv",
    //                           "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt.bk/gar_level_2_worker_0_e_set.csv",
    //                           "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt.bk/gar_level_2_worker_0_x_set.csv",
    //                           "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt.bk/gar_level_2_worker_0_y_set.csv");
    
    // for (int gar_idx = 0; gar_idx < gar_list.size(); gar_idx++){
    //   gar::WriteGAR(gar_list[gar_idx],
    //              "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt.bk/_gar_level_2_worker_" + std::to_string(gar_idx) + "_v.csv",
    //              "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt.bk/_gar_level_2_worker_" + std::to_string(gar_idx) + "_e.csv",
    //              "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt.bk/_gar_level_2_worker_" + std::to_string(gar_idx) + "_x.csv",
    //              "./gar_dbpedia_5_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt.bk/_gar_level_2_worker_" + std::to_string(gar_idx) + "_y.csv");
    // }

    // gar_list.clear();
    // gar::ReadGARSet(gar_list, "./gar_dbpedia_2_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/gar_level_2_worker_0_v_set.csv",
    //                           "./gar_dbpedia_2_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/gar_level_2_worker_0_e_set.csv",
    //                           "./gar_dbpedia_2_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/gar_level_2_worker_0_x_set.csv",
    //                           "./gar_dbpedia_2_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/gar_level_2_worker_0_y_set.csv");
    
    // for (int gar_idx = 0; gar_idx < gar_list.size(); gar_idx++){
    //   gar::WriteGAR(gar_list[gar_idx],
    //              "./gar_dbpedia_2_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/_gar_level_2_worker_" + std::to_string(gar_idx) + "_v.csv",
    //              "./gar_dbpedia_2_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/_gar_level_2_worker_" + std::to_string(gar_idx) + "_e.csv",
    //              "./gar_dbpedia_2_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/_gar_level_2_worker_" + std::to_string(gar_idx) + "_x.csv",
    //              "./gar_dbpedia_2_1_triple_larger_100000_supp_1000_with_constant_with_full_restrictions_opt/_gar_level_2_worker_" + std::to_string(gar_idx) + "_y.csv");
    // }

    // DataGraph data_graph;

    // MovieLensToGraph("../dataset/ml-latest-small/ratings.csv",
    //                  "../dataset/ml-latest-small/movies.csv", data_graph);

    // GUNDAM::WriteCSVGraph(data_graph, "../dataset/ml-latest/movie_lens_v.csv",
    //                                   "../dataset/ml-latest/movie_lens_e.csv");

    // YelpToGraph("../dataset/yelp_dataset/yelp_academic_dataset_business.csv",
    //             "../dataset/yelp_dataset/yelp_academic_dataset_user.csv",
    //             "../dataset/yelp_dataset/yelp_academic_dataset_review.csv",
    //             "../dataset/yelp_dataset/yelp_academic_dataset_tip.csv",
    //             data_graph, "Austin");

    // GUNDAM::WriteCSVGraph(data_graph, "../dataset/yelp_dataset/yelp_v.csv",
    //                                   "../dataset/yelp_dataset/yelp_e.csv");
    
    // gar::ReadGARSet(gar_list, "/Users/apple/Desktop/input/vertex_list.csv",
    //                           "/Users/apple/Desktop/input/edge_list.csv",
    //                           "/Users/apple/Desktop/input/x_list.csv",
    //                           "/Users/apple/Desktop/input/x_list.csv");
    // gar::WriteGARSet(gar_list, "/Users/apple/Desktop/output/vertex_list.csv",
    //                            "/Users/apple/Desktop/output/edge_list.csv",
    //                            "/Users/apple/Desktop/output/x_list.csv",
    //                            "/Users/apple/Desktop/output/y_list.csv");

    // WriteTriple<Pattern, DataGraph>("./gar_yago_5_2_bfs_1_0.01_edge_type_larger_1000_supp_1000_gar_set/gar_level_1_worker_0");
    // WriteTriple<Pattern, DataGraph>("./gar_yago_5_2_bfs_1_0.01_edge_type_larger_1000_supp_1000_gar_set/gar_level_2_worker_0");
    // WriteTriple<Pattern, DataGraph>("./gar_yago_5_2_bfs_1_0.01_edge_type_larger_1000_supp_1000_gar_set/gar_level_3_worker_0");
    // WriteTriple<Pattern, DataGraph>("./gar_yago_5_2_bfs_1_0.01_edge_type_larger_1000_supp_1000_gar_set/gar_level_4_worker_0");
    // WriteTriple<Pattern, DataGraph>("./gar_yago_5_2_bfs_1_0.01_edge_type_larger_1000_supp_1000_gar_set/gar_level_5_worker_0");
    
    // GUNDAM::ReadCSVGraph(pattern, ctx.pattern_v_file_, ctx.pattern_e_file_);

    // auto match_num = GUNDAM::DPISO(pattern, index_graph, -1);

    // LOG(INFO) << "fid = " << frag.fid() << " match num  = " << match_num;
#ifdef PROFILING
    ctx.exec_time += GetCurrentTime();
    ctx.postprocess_time -= GetCurrentTime();
#endif

    messages.ForceContinue();

#ifdef PROFILING
    ctx.postprocess_time += GetCurrentTime();
#endif
  }

  /**
   * @brief Incremental evaluation for SSSP.
   *
   * @param frag
   * @param ctx
   * @param messages
   */
  void IncEval(const fragment_t& frag, context_t& ctx,
               message_manager_t& messages) {
    auto inner_vertices = frag.InnerVertices();

    auto& channels = messages.Channels();

#ifdef PROFILING
    ctx.preprocess_time -= GetCurrentTime();
#endif
    // ctx.next_modified.parallel_clear(thread_num());

#ifdef PROFILING
    ctx.preprocess_time += GetCurrentTime();
    ctx.exec_time -= GetCurrentTime();
#endif

#ifdef PROFILING
    ctx.exec_time += GetCurrentTime();
    ctx.postprocess_time -= GetCurrentTime();
#endif
#ifdef PROFILING
    ctx.postprocess_time += GetCurrentTime();
#endif
  }
//   using Pattern = GUNDAM::SmallGraph<uint32_t, uint32_t, uint32_t, uint32_t>;
  using Pattern = GUNDAM::LargeGraph2<int, int, std::string, 
                                      int, int, std::string>;
  Pattern pattern;
};

}  // namespace grape

#endif  // EXAMPLES_ANALYTICAL_APPS_SSSP_SSSP_H_
