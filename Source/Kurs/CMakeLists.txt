cmake_minimum_required(VERSION 2.6)
set(BIN Kurs)
project(${BIN} CXX)
set(SOURCE_ROOT ../src)
aux_source_directory(${SOURCE_ROOT} SOURCES)


#macro(add_app_executable BUILD_TYPE)
#    if("${BUILD_TYPE}" STREQUAL "debug")
#        set(CMAKE_CXX_FLAGS "-o0 -g")
#        add_executable( ${BIN}-${BUILD_TYPE} ${SOURCES} )
#    endif("${BUILD_TYPE}" STREQUAL "debug")
#    if("${BUILD_TYPE}" STREQUAL "release")
#        set(CMAKE_CXX_FLAGS "-o3 -march='athlon' -funsafe-math-optimizations -fforce-addr -fomit-frame-pointer -foptimize-sibling-calls -fno-stack-protector -finline-small-functions -ftree-parallelize-loops=2 -fprefetch-loop-arrays -funsafe-loop-optimizations -floop-parallelize-all -ftree-loop-linear -funroll-loops -ftree-loop-im -fopenmp -std=c++0x -mpreferred-stack-boundary=2")
#        add_executable( ${BIN}-${BUILD_TYPE} ${SOURCES} )
#    endif("${BUILD_TYPE}" STREQUAL "release")
#endmacro()

#add_app_executable("debug")
#add_app_executable("release")
if(CMAKE_BUILD_TYPE Debug)
    set(CMAKE_CXX_FLAGS "-o0 -g")
endif(CMAKE_BUILD_TYPE Debug)

if(CMAKE_BUILD_TYPE Release)
    set(CMAKE_CXX_FLAGS "-o3 -march='athlon' -funsafe-math-optimizations -fforce-addr -fomit-frame-pointer -foptimize-sibling-calls -fno-stack-protector -finline-small-functions -ftree-parallelize-loops=2 -fprefetch-loop-arrays -funsafe-loop-optimizations -floop-parallelize-all -ftree-loop-linear -funroll-loops -ftree-loop-im -fopenmp -std=c++0x -mpreferred-stack-boundary=2")
endif(CMAKE_BUILD_TYPE Debug)

add_executable( ${BIN}-${CMAKE_BUILD_TYPE} ${SOURCES} )
