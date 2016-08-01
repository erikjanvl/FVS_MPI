cmake_minimum_required(VERSION 3.5)
project(FeedbackVertexSet)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(Lemon REQUIRED)
include_directories(SYSTEM ${LEMON_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11")

set(SOURCE_FILES
        main.cpp
        bafna_et_al.hpp
        componentFinder.hpp
        fvs_compression.hpp
        fvs_optimization.hpp
        IAlgrotihm.hpp
        inplace.hpp
        newgraph.hpp
        reductions.cpp
        reductions.hpp
        uf.hpp
        )
add_executable(FeedbackVertexSet ${SOURCE_FILES})
target_link_libraries(FeedbackVertexSet ${LEMON_LIBRARIES})