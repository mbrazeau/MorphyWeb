SET(mfl_SRC rearrange.c taxpart.c tree.c readnewick.c compare.c coptim.c randtree.c mfyinterface.c cdata.c)
SET_SOURCE_FILES_PROPERTIES(${mfl_SRC} PROPERTIES LANGUAGE CXX )
add_definitions(-D__STDC_LIMIT_MACROS)
add_library (mfl ${mfl_SRC})
# Uncomment this to set the MFY_DEBUG flag on debug builds...
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DMFY_DEBUG")

