file(GLOB nui_SRC
    "*.h"
    "*.cpp"
)
link_directories(${nui_BINARY_DIR}/mfl)
include_directories(${nui_SOURCE_DIR}/mfl)
add_executable(nui ${nui_SRC})
target_link_libraries(nui mfl ${PATH_LIB_NCL} ${PATH_LIB_EDIT} termcap)

