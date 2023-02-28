find_path(HiGHS_INCLUDE_DIR Highs.h
  PATHS
    ENV HiGHS_ROOT
    ENV HiGHS_INCLUDE_DIR
    ${HiGHS_ROOT}
    /usr
    /usr/local
  PATH_SUFFIXES
    include
    include/highs
  )

find_library(HiGHS_LIBRARY
  NAMES
    highs
  PATHS
    ENV HiGHS_ROOT
    ENV HiGHS_LIB_DIR
    ${GMP_ROOT}
    /usr
    /usr/local
  PATH_SUFFIXES
    lib
  )

mark_as_advanced(HiGHS_INCLUDE_DIR HiGHS_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HiGHS
  REQUIRED_VARS
    HiGHS_INCLUDE_DIR
    HiGHS_LIBRARY
  )

if(HiGHS_FOUND AND NOT TARGET HiGHS::HiGHS)
  add_library(HiGHS::HiGHS UNKNOWN IMPORTED)
  set_target_properties(HiGHS::HiGHS PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
    IMPORTED_LOCATION "${HiGHS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${HiGHS_INCLUDE_DIR}"
    )
endif()