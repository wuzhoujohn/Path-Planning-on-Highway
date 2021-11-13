project(Path_Planning)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(sources src/main.cpp)


message("system name is  ${CMAKE_SYSTEM_NAME}")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 


include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1.11.0/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

set(VCPKG_INSTALL ${_VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET})
    message("Using vcpkg_install ${VCPKG_INSTALL}")
    include_directories(${VCPKG_INSTALL}/include)
    link_directories(${VCPKG_INSTALL}/lib)

add_executable(path_planning ${sources})

target_link_libraries(path_planning C:\\Path-Planning-on-Highway\\lib\\zlib.lib C:\\Path-Planning-on-Highway\\lib\\libssl.lib C:\\Path-Planning-on-Highway\\lib\\libuv.lib C:\\Path-Planning-on-Highway\\lib\\uWS.lib C:\\Path-Planning-on-Highway\\lib\\WS2_32.lib)
