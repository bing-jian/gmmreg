#=========================================================================
#Program:   Pointset Registration using Gaussian Mixture Model
#Module:    $RCSfile: CMakeLists.txt,v $
#Language:  C++
#Author:    $Author: bing.jian $
#Date:      $Date: 2008-11-06 01:51:00 -0500 (Thu, 06 Nov 2008) $
#Version:   $Revision: 96 $
#=========================================================================

PROJECT(CVGMI_API)

# The CMake command FIND_PACKAGE(VXL) attempts to find VXL binary 
# installation.  CMake will look in the directory specified by the 
# CMake variable VXL_DIR.  Normally, CMake will initially not find 
# VXL, will warn you that it could not find VXL, and then give you the
# chance to set the variable VXL_DIR and reconfigure.  VXL_DIR now 
# replaces VXL_BINARY_PATH.

# The rest of this CMakeLists.txt file could contain commands as seen 
# in the previous examples, or might just be SUBDIRS commands.

OPTION(BUILD_API_SHARED "Build the CVGMI_API library shared." OFF)
IF(BUILD_API_SHARED)
  SET(API_LIBRARY_TYPE SHARED)
ELSE(BUILD_API_SHARED)
  SET(API_LIBRARY_TYPE STATIC)
ENDIF(BUILD_API_SHARED)


IF(WIN32)
  ADD_DEFINITIONS(-DWIN32 -D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
ENDIF(WIN32)

SET(CVGMI_API_SRCS
         DistanceMatrix.c
         GaussTransform.c
)

ADD_LIBRARY(cvgmi ${API_LIBRARY_TYPE} ${CVGMI_API_SRCS})
IF(UNIX)
  TARGET_LINK_LIBRARIES(cvgmi m)
ENDIF(UNIX)
#SET(LIBRARY_OUTPUT_PATH ${CVGMI_API_SOURCE_DIR}/output/)
#SET(CMAKE_INSTALL_PREFIX ${CVGMI_API_SOURCE_DIR})

INSTALL(TARGETS cvgmi
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
INSTALL(FILES cvgmi_API.h DESTINATION include)

