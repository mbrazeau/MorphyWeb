set(EDIT_INC_FILE "histedit.h")
set(EDIT_LIB_FILE "libedit.a")

find_path(PATH_INC_EDIT NAMES ${EDIT_INC_FILE} PATHS env PATH_INC_EDIT)

macro(editerr filename envvar)
    message(FATAL_ERROR
    "Unable to locate libedit: '${filename}' "
    "Either install libedit to a standard location or set the "
    "env variable '${envvar}' to the location of '${filename}'")
endmacro()

if( NOT PATH_INC_EDIT )
    editerr(${EDIT_INC_FILE} "PATH_INC_EDIT")
endif()

find_library(PATH_LIB_EDIT NAMES ${EDIT_LIB_FILE} PATHS ENV PATH_LIB_EDIT)

if( NOT PATH_LIB_EDIT )
    editerr(${EDIT_LIB_FILE} "PATH_LIB_EDIT")
endif()

message(STATUS "PATH_INC_EDIT = \"${PATH_INC_EDIT}\"")
message(STATUS "PATH_LIB_EDIT = \"${PATH_LIB_EDIT}\"")
