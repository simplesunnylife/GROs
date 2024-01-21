# add by wangyj
cmake_minimum_required(VERSION 2.8)
# IF(APPLE)
#     SET(CMAKE_C_COMPILER "/usr/local/Cellar/gcc/9.1.0/bin/gcc-9")
#     SET(CMAKE_CXX_COMPILER "/usr/local/Cellar/gcc/9.1.0/bin/g++-9")
# ENDIF()
project(libgrape-lite C CXX)

set( CMAKE_CXX_FLAGS "-std=c++11" )

# ------------------------------------------------------------------------------
# cmake options
# ------------------------------------------------------------------------------
option(USE_JEMALLOC "Whether to use jemalloc." OFF)
option(USE_HUGEPAGES "Whether to use hugepages." OFF)
option(BUILD_SHARED_LIBS "Whether to build libgrape-lite as shared library" ON)
option(PROFILING "Whether to enable profiling" OFF)
option(WITH_ASAN "Build with Address Sanitizer" OFF)

if (USE_HUGEPAGES AND LINUX)
    add_definitions(-DUSE_HUGEPAGES)
endif ()

if (PROFILING)
    message("-- Enable profiling")
    add_definitions(-DPROFILING)
endif ()

include_directories(include/libgrape-lite/thirdparty)
include_directories(include/GUNDAM/include/)
include_directories(include/GUNDAM/)
include_directories(include/libgrape-lite/)
include_directories(include/)

# ------------------------------------------------------------------------------
# setting default cmake type to Release
# ------------------------------------------------------------------------------

SET(DEFAULT_BUILD_TYPE "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
            STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

message("[libgrape-lite] will build in type: ${CMAKE_BUILD_TYPE}")

# ------------------------------------------------------------------------------
# cmake configs
# ------------------------------------------------------------------------------

include(CheckLibraryExists)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include_directories(${CMAKE_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (APPLE)
    SET(CMAKE_MACOSX_RPATH ON)
else ()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Werror -Wl,-rpath,$ORIGIN")
endif ()
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -fprofile-arcs -ftest-coverage")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -g -m64 -flto=auto -march=native")

if (WITH_ASAN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -O1")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

# ------------------------------------------------------------------------------
# find_libraries
# ------------------------------------------------------------------------------
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})

# find Threads------------------------------------------------------------------
set(CMAKE_THREAD_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)

# find glog---------------------------------------------------------------------
include("cmake/FindGlog.cmake")
if (NOT GLOG_FOUND)
    message(FATAL_ERROR "glog not found, please install the glog library")
else ()
    include_directories(SYSTEM ${GLOG_INCLUDE_DIRS})
endif ()

# find gflags-------------------------------------------------------------------
Include("cmake/FindGFlags.cmake")
#find_package(gflags CONFIG REQUIRED)
if (NOT GFLAGS_FOUND)
    message(STATUS "gflags not found, build without gflags")
else ()
    include_directories(SYSTEM ${GFLAGS_INCLUDE_DIRS})
endif ()

# find yaml-------------------------------------------------------------------
Include("cmake/FindYaml.cmake")
if (NOT YAML_FOUND)
    message(STATUS "yaml not found, build without yaml")
else ()
    include_directories(SYSTEM ${YAML_INCLUDES})
endif ()

#openmp
find_package(OpenMP)
if(OPENMP_FOUND)
  message("OpenMP FOUND")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# find jemalloc-----------------------------------------------------------------
if (USE_JEMALLOC)
    include("cmake/FindJemalloc.cmake")
    if (NOT JEMALLOC_FOUND)
        message(STATUS "jemalloc not found, build without jemalloc")
    else ()
        add_definitions(-DUSE_JEMALLOC)
        include_directories(SYSTEM ${JEMALLOC_INCLUDE_DIRS})
    endif ()
endif ()

# find rdkafka---------------------------------------------------------------------
include("cmake/FindRdkafka.cmake")
if (NOT RDKAFKA_FOUND)
    message(STATUS "rdkafka not found, build without rdkafka")
endif ()

# ------------------------------------------------------------------------------
# generate libgrape-lite
# ------------------------------------------------------------------------------
file(GLOB_RECURSE CORE_SRC_FILES "include/libgrape-lite/grape/*.cc")
message("wenzhi CORE_SRC_FILES: ${CORE_SRC_FILES}")
add_library(grape-lite ${CORE_SRC_FILES})

target_link_libraries(grape-lite ${MPI_CXX_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
        ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES})

if (JEMALLOC_FOUND)
    target_link_libraries(grape-lite ${JEMALLOC_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable analytical_apps because gflags not found")
else ()
    add_executable(analytical_apps src/apps/flags.cc src/apps/run_app.cc)
    target_include_directories(analytical_apps PRIVATE src/apps)
    set_target_properties(analytical_apps PROPERTIES OUTPUT_NAME run_app)
    target_link_libraries(analytical_apps grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable gar_discover because gflags not found")
else ()
    add_executable(gar_discover src/apps/flags.cc src/apps/gar_discover.cc)
    target_include_directories(gar_discover PRIVATE src/apps)
    set_target_properties(gar_discover PROPERTIES OUTPUT_NAME gar_discover)
    target_link_libraries(gar_discover grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable gcr_discover because gflags not found")
else ()
    add_executable(gcr_discover src/apps/flags.cc src/apps/gcr_discover.cc)
    target_include_directories(gcr_discover PRIVATE src/apps)
    set_target_properties(gcr_discover PROPERTIES OUTPUT_NAME gcr_discover)
    target_link_libraries(gcr_discover grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable gar_supp because gflags not found")
else ()
    add_executable(gar_supp src/apps/flags.cc src/apps/gar_supp.cc)
    target_include_directories(gar_supp PRIVATE src/apps)
    set_target_properties(gar_supp PROPERTIES OUTPUT_NAME gar_supp)
    target_link_libraries(gar_supp grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable support_sample because gflags not found")
else ()
    add_executable(support_sample src/apps/flags.cc src/apps/support_sample.cc)
    target_include_directories(support_sample PRIVATE src/apps)
    set_target_properties(support_sample PROPERTIES OUTPUT_NAME support_sample)
    target_link_libraries(support_sample grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable gar_overlap because gflags not found")
else ()
    add_executable(gar_overlap src/apps/flags.cc src/apps/gar_overlap.cc)
    target_include_directories(gar_overlap PRIVATE src/apps)
    set_target_properties(gar_overlap PROPERTIES OUTPUT_NAME gar_overlap)
    target_link_libraries(gar_overlap grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable gar_analysis because gflags not found")
else ()
    add_executable(gar_analysis src/apps/flags.cc src/apps/gar_analysis.cc)
    target_include_directories(gar_analysis PRIVATE src/apps)
    set_target_properties(gar_analysis PROPERTIES OUTPUT_NAME gar_analysis)
    target_link_libraries(gar_analysis grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable gar_accuracy_analysis because gflags not found")
else ()
    add_executable(gar_accuracy_analysis src/apps/flags.cc src/apps/gar_accuracy_analysis.cc)
    target_include_directories(gar_accuracy_analysis PRIVATE src/apps)
    set_target_properties(gar_accuracy_analysis PROPERTIES OUTPUT_NAME gar_accuracy_analysis)
    target_link_libraries(gar_accuracy_analysis grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable gar_match because gflags not found")
else ()
    add_executable(gar_match src/apps/flags.cc src/apps/gar_match.cc)
    target_include_directories(gar_match PRIVATE src/apps)
    set_target_properties(gar_match PROPERTIES OUTPUT_NAME gar_match)
    target_link_libraries(gar_match grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable pattern_match because gflags not found")
else ()
    add_executable(pattern_match src/apps/flags.cc src/apps/pattern_match.cc)
    target_include_directories(pattern_match PRIVATE src/apps)
    set_target_properties(pattern_match PROPERTIES OUTPUT_NAME pattern_match)
    target_link_libraries(pattern_match grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable gor_match because gflags not found")
else ()
    add_executable(gor_match src/apps/flags.cc src/apps/gor_match.cc)
    target_include_directories(gor_match PRIVATE src/apps)
    set_target_properties(gor_match PROPERTIES OUTPUT_NAME gor_match)
    target_link_libraries(gor_match grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable rule_match because gflags not found")
else ()
    add_executable(rule_match src/apps/flags.cc src/apps/rule_match.cc)
    target_include_directories(rule_match PRIVATE src/apps)
    set_target_properties(rule_match PROPERTIES OUTPUT_NAME rule_match)
    target_link_libraries(rule_match grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable gor_accuracy_analysis because gflags not found")
else ()
    add_executable(gor_accuracy_analysis src/apps/flags.cc src/apps/gor_accuracy_analysis.cc)
    target_include_directories(gor_accuracy_analysis PRIVATE src/apps)
    set_target_properties(gor_accuracy_analysis PROPERTIES OUTPUT_NAME gor_accuracy_analysis)
    target_link_libraries(gor_accuracy_analysis grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

add_executable(tools "src/apps/tools.cc")
target_include_directories(tools PRIVATE src/apps)
set_target_properties(tools PROPERTIES OUTPUT_NAME tools)
target_link_libraries(tools ${YAML_LIBRARIES})

add_executable(tgraph_discover "src/apps/tgraph_discover.cc")
target_include_directories(tgraph_discover PRIVATE src/apps)
set_target_properties(tgraph_discover PROPERTIES OUTPUT_NAME tgraph_discover)
target_link_libraries(tgraph_discover PRIVATE ${YAML_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(tgraph_apply "src/apps/tgraph_apply.cc")
target_include_directories(tgraph_apply PRIVATE src/apps)
set_target_properties(tgraph_apply PROPERTIES OUTPUT_NAME tgraph_apply)
target_link_libraries(tgraph_apply PRIVATE ${YAML_LIBRARIES} ${GLOG_LIBRARIES} ${MPI_CXX_LIBRARIES})

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable prob_gar_chase because gflags not found")
else ()
    add_executable(prob_gar_chase src/apps/flags.cc src/apps/prob_gar_chase.cc)
    target_include_directories(prob_gar_chase PRIVATE src/apps)
    set_target_properties(prob_gar_chase PROPERTIES OUTPUT_NAME prob_gar_chase)
    target_link_libraries(prob_gar_chase grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable gor_accuracy_analysis because gflags not found")
else ()
    add_executable(rdf_to_graph src/apps/rdf_to_graph.cc src/apps/rdf_to_graph/flags.cc)
    target_include_directories(rdf_to_graph PRIVATE src/apps)
    set_target_properties(rdf_to_graph PROPERTIES OUTPUT_NAME rdf_to_graph)
    target_link_libraries(rdf_to_graph ${GFLAGS_LIBRARIES} ${YAML_LIBRARIES})
endif ()

add_executable(graph_discretize "src/apps/graph_discretize.cc")
target_include_directories(graph_discretize PRIVATE src/apps)
set_target_properties(graph_discretize PROPERTIES OUTPUT_NAME graph_discretize)
target_link_libraries(graph_discretize ${YAML_LIBRARIES})

add_executable(duplicate_vertex "src/apps/duplicate_vertex.cc")
target_include_directories(duplicate_vertex PRIVATE src/apps)
set_target_properties(duplicate_vertex PROPERTIES OUTPUT_NAME duplicate_vertex)
target_link_libraries(duplicate_vertex ${YAML_LIBRARIES})

# add_subdirectory(test)

#
# if (NOT GFLAGS_FOUND)
#     message(WARNING "Disable gnn_sampler because flags not found")
# elseif (NOT RDKAFKA_FOUND)
#     message(WARNING "Disable gnn_sampler because rdkafka not found")
# else ()
#     add_executable(gnn_sampler examples/gnn_sampler/run_sampler.cc)
#     set_target_properties(gnn_sampler PROPERTIES OUTPUT_NAME run_sampler)
#     target_include_directories(gnn_sampler PRIVATE
#             examples/gnn_sampler
#             examples/gnn_sampler/thirdparty
#             ${RDKAFKA_INCLUDE_DIR})
#     target_link_libraries(gnn_sampler grape-lite  ${MPI_CXX_LIBRARIES}
#             ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${RDKAFKA_LIBRARIES})
# endif ()

# ------------------------------------------------------------------------------
# format code
# ------------------------------------------------------------------------------
#file(GLOB_RECURSE FILES_NEED_FORMAT "include/libgrape-lite/*.cc"
#        "include/libgrape-lite/grape/*.h"
#        "include/libgrape-lite/examples/*.h"
#        "include/libgrape-lite/examples/*.cc")
#foreach (file_path ${FILES_NEED_FORMAT})
#    if (${file_path} MATCHES ".*thirdparty.*")
#        list(REMOVE_ITEM FILES_NEED_FORMAT ${file_path})
#    endif ()
#endforeach ()

#add_custom_target(clformat
#        COMMAND clang-format --style=file -i ${FILES_NEED_FORMAT}
#        COMMENT "Running clang-format."
#        VERBATIM)

# ------------------------------------------------------------------------------
# cpplint, check for readability with Google style
# ------------------------------------------------------------------------------
#add_custom_target(cpplint
#        COMMAND ${CMAKE_SOURCE_DIR}/misc/cpplint.py ${FILES_NEED_FORMAT}
#        COMMENT "Running cpplint check."
#        VERBATIM)

# ------------------------------------------------------------------------------
# generate docs
# ------------------------------------------------------------------------------
#add_custom_target(doc COMMAND doxygen "${CMAKE_SOURCE_DIR}/misc/doc-config"
#        COMMENT "Generating docs."
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#        VERBATIM)
