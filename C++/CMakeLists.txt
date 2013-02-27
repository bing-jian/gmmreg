#=========================================================================
# $Author: bingjian $
# $Date: 2013-01-04 01:39:25 -0500 (Fri, 04 Jan 2013) $
# $Revision: 145 $
#=========================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

PROJECT(GMMREG )

# The CMake command FIND_PACKAGE(VXL) attempts to find VXL binary
# installation.  CMake will look in the directory specified by the
# CMake variable VXL_DIR.  Normally, CMake will initially not find
# VXL, will warn you that it could not find VXL, and then give you the
# chance to set the variable VXL_DIR and reconfigure.  VXL_DIR now
# replaces VXL_BINARY_PATH.

FIND_PACKAGE(VXL)

# Whether FIND_PACKAGE(VXL) worked is stored in the variable
# VXL_FOUND.  If VXL was found we then include `UseVXL.cmake'
# which will set many variables prefixed with VXL_ that your project
# can use to determine what parts of VXL are present and how to use
# them.  `UseVXL.cmake' holds all the necessary definitions that
# CMake needs to use VXL.

IF(VXL_FOUND)
  INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
ENDIF(VXL_FOUND)

# `UseVXL.cmake' will in turn include the file
# `VXLConfig.cmake' from the top of the VXL build tree.  You
# should look at this file to see what CMake variables are imported.

# To set VXL_DIR, you could rely on an environment variable:
# SET( VXL_DIR $ENV{VXLBIN} )
# However, you must then make sure that VXLBIN is set correctly every
# time that cmake is run.  This is often a source of much frustration
# and is not recommended.
#
# Another option is to set the path directly:
#   SET( VXL_DIR/usr/local/vxl-1.0-beta2/bin )
# This has the drawback that the CMakeLists.txt file (this file) has
# to be modified to suit your particular setup, which may cause
# trouble.  If this file is under CVS control, for example, you may
# inadvertently commit it with your local setup hard-coded, thus
# playing havoc with other users.

# The rest of this CMakeLists.txt file could contain commands as seen
# in the previous examples, or might just be SUBDIRS commands.


# Remember that all CMake variables defining the VXL configuration
# begin with "VXL_".  If you want to use VGUI in a client project you
# should use something like this.

#OPTION(BUILD_PYTHON_EXT "Build the Python extension." OFF)
#OPTION(BUILD_MATLAB_MEX "Build the Matlab MEX file." OFF)
OPTION(BUILD_INI_SHARED "Build the port_ini libraries shared." OFF)
OPTION(BUILD_API_SHARED "Build the gmmreg_api libraries shared." OFF)

#OPTION(BUILD_RIGID3D_SHARED "Build the gmmreg_rigid3d libraries shared." OFF)

IF(BUILD_INI_SHARED)
  SET(INI_LIBRARY_TYPE SHARED)
ELSE(BUILD_INI_SHARED)
  SET(INI_LIBRARY_TYPE STATIC)
ENDIF(BUILD_INI_SHARED)

IF(BUILD_API_SHARED)
  SET(API_LIBRARY_TYPE SHARED)
ELSE(BUILD_API_SHARED)
  SET(API_LIBRARY_TYPE STATIC)
ENDIF(BUILD_API_SHARED)


IF(WIN32)
  ADD_DEFINITIONS(-DWIN32 -D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
ELSE(WIN32)
  ADD_LIBRARY(port_ini ${INI_LIBRARY_TYPE} port_ini.c)
ENDIF(WIN32)


SET(GMMREG_API_SRCS
         gmmreg_utils.cpp
         gmmreg_base.cpp
         gmmreg_cpd.cpp
         gmmreg_tps.cpp
         gmmreg_tps_func.cpp
         gmmreg_grbf.cpp
         gmmreg_grbf_func.cpp
         gmmreg_rigid.cpp
         gmmreg_rigid_func.cpp
         gmmreg_api.cpp
)
file(GLOB_RECURSE GMMREG_API_HDRS *.h)

#IF( VXL_VNL_FOUND )
INCLUDE_DIRECTORIES( ${VXL_VNL_INCLUDE_DIR} )

SET(gauss_transform_sources ${gauss_transform_sources}
test_gauss_transform.cpp gmmreg_utils.cpp
)
ADD_EXECUTABLE(test_gauss_transform ${gauss_transform_sources}
)
#TARGET_LINK_LIBRARIES( test_gauss_transform ${VXL_VNL_LIBRARIES} )
TARGET_LINK_LIBRARIES(test_gauss_transform vnl_algo vnl vcl)

ADD_EXECUTABLE(gmmreg_aux gmmreg_aux.cpp
)
TARGET_LINK_LIBRARIES(gmmreg_aux vnl vnl_algo)

ADD_EXECUTABLE(extract_correspondence
    extract_correspondence.cpp gmmreg_utils.cpp
)
TARGET_LINK_LIBRARIES(extract_correspondence vnl vnl_algo)

ADD_LIBRARY(gmmreg_api ${API_LIBRARY_TYPE} ${GMMREG_API_SRCS} ${GMMREG_API_HDRS})
TARGET_LINK_LIBRARIES(gmmreg_api vnl_algo vnl vcl)
IF(UNIX)
  TARGET_LINK_LIBRARIES(gmmreg_api m port_ini)
ENDIF(UNIX)


SET(gmmreg_demo_sources ${gmmreg_demo_sources}
  gmmreg_demo.cpp
)

ADD_EXECUTABLE(gmmreg_demo ${gmmreg_demo_sources}
)
TARGET_LINK_LIBRARIES(gmmreg_demo gmmreg_api)

SET(gmmreg_transform_sources ${gmmreg_transform_sources}
  gmmreg_transform.cpp
  gmmreg_utils.cpp
)

ADD_EXECUTABLE(gmmreg_transform ${gmmreg_transform_sources}
)
TARGET_LINK_LIBRARIES(gmmreg_transform vnl vnl_algo vcl)
IF(UNIX)
  TARGET_LINK_LIBRARIES(gmmreg_transform m port_ini)
ENDIF(UNIX)


#ENDIF( VXL_VNL_FOUND )
