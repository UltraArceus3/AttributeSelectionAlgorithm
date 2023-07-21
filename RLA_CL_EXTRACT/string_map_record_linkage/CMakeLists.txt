cmake_minimum_required(VERSION 3.7...3.22)


if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(StringMap VERSION 1.0
                  DESCRIPTION "String Mapping"
                  LANGUAGES CXX)

include_directories(include)
file(GLOB SOURCES "src/*.cpp")
add_executable(stringMap ${SOURCES})

find_package(LibXml2 REQUIRED)
INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})

find_package(Boost 1.40 COMPONENTS regex REQUIRED)
INCLUDE_DIRECTORIES(${boost_INCLUDE_DIR})

find_package(Boost 1.71 COMPONENTS log REQUIRED)

#specify the c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

target_link_libraries(stringMap ${MPI_CXX_LIBRARIES} ${LIBXML2_LIBRARIES} ${Boost_LIBRARIES} )



