CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(PROJECT_NAME mh_lib)
PROJECT(${PROJECT_NAME})

MESSAGE(STATUS "CMAKE_C_COMPILER:   ${CMAKE_C_COMPILER}" )
MESSAGE(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}" )

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
        message(STATUS "C++11 activated by adding -std=gnu++11")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
    elseif(GCC_VERSION VERSION_GREATER 4.3 OR GCC_VERSION VERSION_EQUAL 4.3)
        message(WARNING "C++0x activated. If you get any errors update to a compiler which fully supports C++11")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
    else()
        message(FATAL_ERROR "C++11 needed. Therefore a gcc compiler with a version higher than 4.3 is needed.")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message(STATUS "Clang detected, adding -stdlib=libc++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
else()
    message(WARNING "Unknown compiler ${CMAKE_CXX_COMPILER_ID}, your compiler may not support C++11.")
endif()

SET(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS}         -Wall")
SET(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -Wall -g -pg -DDEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE "Release")
ENDIF()

SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
SET(CMAKE_EXPORT_COMPILE_COMMANDS 1)

INCLUDE_DIRECTORIES(include)

### LIBRARIES ###

### EIGEN
#set(Eigen3_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern )
#FIND_PACKAGE(Eigen3 REQUIRED)
set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern )
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})


### OPENGL
FIND_PACKAGE(OpenGL REQUIRED)
INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_DIR})
LIST(APPEND LIB_LIBRARIES ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )
#LINK_LIBRARIES(${OPENGL_gl_LIBRARY})
#LINK_LIBRARIES(${OPENGL_glu_LIBRARY})

### GLFW
# Assuming glfw was built and installed in ../extern
set( GLFW_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern/glfw-3.1.2/install" )
if(WIN32)
    set(GLFW3_INCLUDE_DIR "${GLFW_INSTALL_DIR}/include")
    set(GLFW3_LIBRARY_DIR "${GLFW_INSTALL_DIR}/lib" )
    set(GLFW3_LIBRARY "glfw3dll")
    set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${GLFW3_LIBRARY_DIR}")
else(WIN32) # unix
    set(glfw3_DIR "${GLFW_INSTALL_DIR}/lib/cmake/glfw" CACHE PATH "glfw source path" )
    FIND_PACKAGE( glfw3 REQUIRED )
    if (NOT GLFW3_INCLUDE_DIR )
        set( GLFW3_INCLUDE_DIR ${GLFW_INSTALL_DIR}/include CACHE PATH "GLFW3_INCLUDE_DIR" FORCE )
    endif()
endif()
INCLUDE_DIRECTORIES(${GLFW3_INCLUDE_DIR})
LIST(APPEND LIB_LIBRARIES ${GLFW3_LIBRARY} )

### GLEW
if(WIN32)
    set(GLEW_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern/glew-1.13.0/")
    set(GLEW_INSTALL_DIR "${GLEW_DIR}")
    set(GLEW_LIBRARIES "glew32")
    set(GLEW_INCLUDE_DIRS "${GLEW_INSTALL_DIR}/include")
    set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${GLFW_LIBRARY_DIR};${GLEW_INSTALL_DIR}/lib/Release/x64")
else(WIN32)
    SET(GLEW_STATIC)
    FIND_PACKAGE(GLEW REQUIRED)
    MESSAGE(STATUS "GLEW_INCLUDE_DIRS: ${GLEW_INCLUDE_DIRS}" )
    MESSAGE(STATUS "GLEW_LIBRARIES: ${GLEW_LIBRARIES}" )
endif(WIN32)
INCLUDE_DIRECTORIES(${GLEW_INCLUDE_DIRS})
LIST(APPEND LIB_LIBRARIES ${GLEW_LIBRARIES})

### ASSIMP
FIND_PACKAGE(assimp)
INCLUDE_DIRECTORIES(${ASSIMP_INCLUDE_DIRS})
IF( assimp_LIBRARIES )
    add_definitions(-DUSE_ASSIMP)
    list(APPEND LIB_LIBRARIES ${assimp_LIBRARIES})
ENDIF( assimp_LIBRARIES )

set(GSL_LIBR ${CMAKE_CURRENT_SOURCE_DIR}/extern/gsl/lib )
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${GSL_LIBR}")
link_directories(${GSL_LIBR})

set(OpenMesh_LIBR ${CMAKE_CURRENT_SOURCE_DIR}/extern/OpenMesh/lib )
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${OpenMesh_LIBR}")
link_directories(${OpenMesh_LIBR})

list(APPEND LIB_LIBRARIES ${GSL_LIBR}/libgsl.a 
  ${GSL_LIBR}/libgslcblas.a 
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/ANN/lib/libANN.a  
 # ${OpenMesh_LIBR}/libOpenMeshCore.a
 # ${OpenMesh_LIBR}/libOpenMeshTools.a
   ${OpenMesh_LIBR}/libOpenMeshCore.6.0.dylib
   ${OpenMesh_LIBR}/libOpenMeshTools.6.0.dylib
  
 
)
MESSAGE(STATUS “om_lib: “ ${OpenMesh_LIBR}/libOpenMeshCore.2.2.dylib)
#list(APPEND LIB_LIBRARIES  ${CMAKE_CURRENT_SOURCE_DIR}/extern/ANN/lib/libANN.a)

### TARGET MH_LIB ###

if(WIN32)
    link_directories("${CMAKE_PREFIX_PATH}")
    MESSAGE(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}" )
endif(WIN32)

### SRC FILES
FILE(GLOB_RECURSE LIB_SRC_FILES    ${PROJECT_SOURCE_DIR}/src/*.cpp)
FILE(GLOB_RECURSE LIB_HEADER_FILES ${PROJECT_SOURCE_DIR}/include/*.h
                                   ${PROJECT_SOURCE_DIR}/include/*.hpp)

#ADD_LIBRARY(mh SHARED ${LIB_SRC_FILES} ${EXTRA_SRC_FILES} ${LIB_HEADER_FILES})
ADD_LIBRARY(mh ${LIB_SRC_FILES} ${EXTRA_SRC_FILES} ${LIB_HEADER_FILES})
TARGET_LINK_LIBRARIES(mh ${LIB_LIBRARIES} )


### TARGET MESHVIEWER ###

SET( MESHVIEWER_SRCS
    example/meshviewer/gui/meshviewer.cpp
    example/meshviewer/gui/meshviewer.h
    example/meshviewer/gui/assignment2.cpp
    #example/meshviewer/gui/assignment2.h
    example/meshviewer/main.cpp
    )

INCLUDE_DIRECTORIES( example/meshviewer )


ADD_EXECUTABLE( meshviewer ${MESHVIEWER_SRCS} )
TARGET_LINK_LIBRARIES( meshviewer mh)
MESSAGE(STATUS CMAKE_PREFIPATH: ${CMAKE_PREFIX_PATH})


### HELPER TARGET FOR HEADERS ###

FILE(GLOB_RECURSE LIB_HEADER_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
FILE(GLOB_RECURSE MESHVIEWER_HEADER_FILES meshviewer/*.h)
add_custom_target( HEADERS_TARGET ${LIB_HEADER_FILES} ${MESHVIEWER_HEADER_FILES} )
set_target_properties( HEADERS_TARGET PROPERTIES LINKER_LANGUAGE CXX )

