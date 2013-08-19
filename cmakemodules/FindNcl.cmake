find_path(PATH_INC_NCL NAMES "ncl/ncl.h" PATHS)

if( NOT PATH_INC_NCL )
        message(FATAL_ERROR "Unable to locate Nexus Class Library include files" )
endif( NOT PATH_INC_NCL )

find_library(PATH_LIB_NCL NAMES "ncl/libncl.a" PATHS)

if( NOT PATH_LIB_NCL )
        message(FATAL_ERROR "Unable to locate Nexus Class Library file" )
endif( NOT PATH_LIB_NCL )

MESSAGE( STATUS "PATH_INC_NCL = \"${PATH_INC_NCL}\"" )
MESSAGE( STATUS "PATH_LIB_NCL = \"${PATH_LIB_NCL}\"" )
