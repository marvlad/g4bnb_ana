cmake_minimum_required(VERSION 3.10)
project(g4bnb_ana)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ROOT
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIRS})

# dk2nu include + lib
include_directories(/exp/annie/app/users/mascenci/g4bnb_detached/dk2nu/install/include)
link_directories(/exp/annie/app/users/mascenci/g4bnb_detached/dk2nu/install/lib)

# includes for your headers
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(analysis
    main.cpp
    src/FileManager.cpp
    src/HistogramManager.cpp
    src/Dk2nuAnalyzer.cpp
)

target_link_libraries(analysis
    ${ROOT_LIBRARIES}
    dk2nuTree    # official dk2nu library
)
