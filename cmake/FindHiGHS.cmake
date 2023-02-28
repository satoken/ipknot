# インクルードディレクトリのパスをgmp.hを頼りに検索する
find_path(HiGHS_INCLUDE_DIR Highs.h
  # 検索するパス
  PATHS
    # 環境変数GMP_ROOTまたはGMP_INCLUDE_DIRが存在したらそこを検索
    ENV HiGHS_ROOT
    ENV HiGHS_INCLUDE_DIR
    # CMakeの変数としてGMP_ROOTが定義されていたらそこを検索
    ${HiGHS_ROOT}
    /usr
    /usr/local
  PATH_SUFFIXES
    include
    include/highs
  )
# ライブラリへのパスをライブラリ名を元に検索する
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
# advancedモードでない限り変数の存在を隠す
mark_as_advanced(HiGHS_INCLUDE_DIR HiGHS_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HiGHS
  REQUIRED_VARS
    HiGHS_INCLUDE_DIR
    HiGHS_LIBRARY
  )

# GMPが見つかり、かつGMP::GMPが定義されていない場合
if(HiGHS_FOUND AND NOT TARGET HiGHS::HiGHS)
  # GMP::GMPというターゲット名でGMPライブラリを定義
  # UNKNOWN = STATIC/SHAREDかはまだ不明
  # IMPORTED = このプロジェクトに属さないターゲット
  add_library(HiGHS::HiGHS UNKNOWN IMPORTED)
  set_target_properties(HiGHS::HiGHS PROPERTIES
    # C言語、C++なら"CXX"とする
    IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
    IMPORTED_LOCATION "${HiGHS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${HiGHS_INCLUDE_DIR}"
    )
endif()