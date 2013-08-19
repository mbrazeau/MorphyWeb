cmake_minimum_required (VERSION 2.6)
project (nui)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${nui_SOURCE_DIR}/cmakemodules")

# Find gsl, libedit, and ncl
FIND_PACKAGE(GSL REQUIRED)
FIND_PACKAGE(Libedit REQUIRED)
FIND_PACKAGE(Ncl REQUIRED)

# Crank up the warning levels
if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

add_subdirectory (mfl)
add_subdirectory (nui)
# things to check for:
# ncl
# histedit.h
