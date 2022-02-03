cmake_minimum_required(VERSION 3.0)

project(EE LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")

# build
set(ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(SRC ${ROOT}/src)
set(DEPS ${ROOT}/deps)
set(BIN ${ROOT}/bin)

add_subdirectory(${DEPS}/json)

set(EXECUTABLE_OUTPUT_PATH ${BIN})
file(GLOB_RECURSE sources "${SRC}/*.h" "${SRC}/*.cpp")
add_executable(EE ${sources})
target_include_directories(EE PRIVATE ${SRC} ${DEPS}/cspice/include ${DEPS}/json/include)
target_link_directories(EE PRIVATE ${DEPS}/cspice/lib)
target_link_libraries(EE PRIVATE -l:cspice.a nlohmann_json::nlohmann_json)