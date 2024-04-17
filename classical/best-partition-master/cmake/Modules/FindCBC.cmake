# The contents of this file have been adapted by Daniëlle Schuman

# - Try to find CBC
# Once done this will define
#  CBC_FOUND           - System has CBC
#  CBC_INCLUDE_DIRS    - The CBC include directories
#  CBC_LIBRARY_PATHS   - The libraries needed to use CBC
#  CBC_TARGETS         - The names of imported targets created for CBC
# User can set CBC_ROOT to the preferred installation prefix

find_path(CBC_FILE_LOC_A NAMES CbcModel.hpp PATHS /usr/local/Cellar/cbc/2.10.7_1/include/cbc/coin)
message(STATUS "CBC: Found file `CbcModel.hpp` at `${CBC_FILE_LOC_A}`")
list(APPEND CBC_INCLUDE_DIRS ${CBC_FILE_LOC_A})

find_path(CBC_FILE_LOC_B NAMES ClpSimplex.hpp PATHS /usr/local/Cellar/clp/1.17.7/include/clp/coin)
message(STATUS "CBC: Found file `ClpSimplex.hpp` at `${CBC_FILE_LOC_B}`")
list(APPEND CBC_INCLUDE_DIRS ${CBC_FILE_LOC_B})

find_path(CBC_FILE_LOC_C NAMES OsiClpSolverInterface.hpp PATHS /usr/local/Cellar/clp/1.17.7/include/clp/coin)
message(STATUS "CBC: Found file `OsiClpSolverInterface.hpp` at `${CBC_FILE_LOC_C}`")
list(APPEND CBC_INCLUDE_DIRS ${CBC_FILE_LOC_C})

find_path(CBC_FILE_LOC_D NAMES OsiSolverInterface.hpp PATHS /usr/local/Cellar/osi/0.108.7_1/include/osi/coin)
message(STATUS "CBC: Found file `OsiSolverInterface.hpp` at `${CBC_FILE_LOC_D}`")
list(APPEND CBC_INCLUDE_DIRS ${CBC_FILE_LOC_D})

find_path(CBC_FILE_LOC_E NAMES CoinPragma.hpp PATHS /usr/local/Cellar/coinutils/2.11.6_1/include/coinutils/coin)
message(STATUS "CBC: Found file `CoinPragma.hpp` at `${CBC_FILE_LOC_E}`")
list(APPEND CBC_INCLUDE_DIRS ${CBC_FILE_LOC_E})

list(REMOVE_DUPLICATES CBC_INCLUDE_DIRS)


find_library(CBC_LIB_LOC_A NAMES libCbcSolver.3.10.7.dylib PATHS /usr/local/Cellar/cbc/2.10.7_1/lib)
message(STATUS "CBC: Found lib CbcSolver at `${CBC_LIB_LOC_A}`")
list(APPEND CBC_LIBRARY_PATHS ${CBC_LIB_LOC_A})
add_library(CbcSolver UNKNOWN IMPORTED)
set_target_properties(CbcSolver PROPERTIES
                        IMPORTED_LOCATION ${CBC_LIB_LOC_A}
                        INTERFACE_INCLUDE_DIRECTORIES "${CBC_INCLUDE_DIRS}")
list(APPEND CBC_TARGETS CbcSolver)

find_library(CBC_LIB_LOC_B NAMES libCbc.3.10.7.dylib PATHS /usr/local/Cellar/cbc/2.10.7_1/lib)
list(APPEND CBC_LIBRARY_PATHS ${CBC_LIB_LOC_B})
add_library(Cbc UNKNOWN IMPORTED)
set_target_properties(Cbc PROPERTIES
                        IMPORTED_LOCATION ${CBC_LIB_LOC_B}
                        INTERFACE_INCLUDE_DIRECTORIES "${CBC_INCLUDE_DIRS}")
list(APPEND CBC_TARGETS Cbc)

find_library(CBC_LIB_LOC_C NAMES libCgl.1.10.5.dylib PATHS /usr/local/Cellar/cgl/0.60.5_1/lib)
list(APPEND CBC_LIBRARY_PATHS ${CBC_LIB_LOC_C})
add_library(Cgl UNKNOWN IMPORTED)
set_target_properties(Cgl PROPERTIES
                        IMPORTED_LOCATION ${CBC_LIB_LOC_C}
                        INTERFACE_INCLUDE_DIRECTORIES "${CBC_INCLUDE_DIRS}")
list(APPEND CBC_TARGETS Cgl)

find_library(CBC_LIB_LOC_D NAMES libOsiClp.1.14.7.dylib PATHS /usr/local/Cellar/clp/1.17.7/lib)
list(APPEND CBC_LIBRARY_PATHS ${CBC_LIB_LOC_D})
add_library(OsiClp UNKNOWN IMPORTED)
set_target_properties(OsiClp PROPERTIES
                        IMPORTED_LOCATION ${CBC_LIB_LOC_D}
                        INTERFACE_INCLUDE_DIRECTORIES "${CBC_INCLUDE_DIRS}")
list(APPEND CBC_TARGETS OsiClp)

find_library(CBC_LIB_LOC_E NAMES libClp.1.14.7.dylib PATHS /usr/local/Cellar/clp/1.17.7/lib)
list(APPEND CBC_LIBRARY_PATHS ${CBC_LIB_LOC_E})
add_library(Clp UNKNOWN IMPORTED)
set_target_properties(Clp PROPERTIES
                        IMPORTED_LOCATION ${CBC_LIB_LOC_E}
                        INTERFACE_INCLUDE_DIRECTORIES "${CBC_INCLUDE_DIRS}")
list(APPEND CBC_TARGETS Clp)

find_library(CBC_LIB_LOC_F NAMES libOsi.1.13.7.dylib PATHS /usr/local/Cellar/osi/0.108.7_1/lib)
list(APPEND CBC_LIBRARY_PATHS ${CBC_LIB_LOC_F})
add_library(Osi UNKNOWN IMPORTED)
set_target_properties(Osi PROPERTIES
                        IMPORTED_LOCATION ${CBC_LIB_LOC_F}
                        INTERFACE_INCLUDE_DIRECTORIES "${CBC_INCLUDE_DIRS}")
list(APPEND CBC_TARGETS Osi)

find_library(CBC_LIB_LOC_G NAMES libCoinUtils.3.11.6.dylib PATHS /usr/local/Cellar/coinutils/2.11.6_1/lib)
list(APPEND CBC_LIBRARY_PATHS ${CBC_LIB_LOC_G})
add_library(CoinUtils UNKNOWN IMPORTED)
set_target_properties(CoinUtils PROPERTIES
                        IMPORTED_LOCATION ${CBC_LIB_LOC_G}
                        INTERFACE_INCLUDE_DIRECTORIES "${CBC_INCLUDE_DIRS}")
list(APPEND CBC_TARGETS CoinUtils)

if(UNIX AND NOT WIN32 AND NOT DEFINED EMSCRIPTEN)
    find_package(ZLIB)
    if(NOT ZLIB_FOUND)
      message(STATUS "CBC: Missing dependency `Zlib`")
      set(CBC_LIBRARY_PATHS "")
    else()
      list(APPEND CBC_LIBRARY_PATHS ${ZLIB_LIBRARIES})
      list(APPEND CBC_TARGETS ${ZLIB_LIBRARIES})
    endif()
endif()

#[[
set(COIN_SUB_DIRS "coin-or" "coin")
set(CBC_FIND_FILES CbcModel.hpp ClpSimplex.hpp OsiClpSolverInterface.hpp OsiSolverInterface.hpp CoinPragma.hpp)

foreach(COIN_SUB_DIR ${COIN_SUB_DIRS})
  foreach(CBC_FILE ${CBC_FIND_FILES})
    set(CBC_FILE_LOC "CBC_LIB_LOC-NOTFOUND")
    message(STATUS "CBC: Looking for file `${COIN_SUB_DIR}/${CBC_FILE}`")
    find_path(CBC_FILE_LOC ${COIN_SUB_DIR}/${CBC_FILE}
              PATH_SUFFIXES cbc cgl clp coinutils osi include)
    if("${CBC_FILE_LOC}" STREQUAL "CBC_FILE_LOC-NOTFOUND")
      message(STATUS "CBC: Could not find file `${CBC_FILE}`")
      set(CBC_INCLUDE_DIRS "")
      break()
    else()
      message(STATUS "CBC: Found file `${CBC_FILE}` at `${CBC_FILE_LOC}`")
    endif()
    list(APPEND CBC_INCLUDE_DIRS ${CBC_FILE_LOC})
    list(APPEND CBC_INCLUDE_DIRS ${CBC_FILE_LOC}/${COIN_SUB_DIR})
  endforeach(CBC_FILE)
  if(NOT "${CBC_FILE_LOC}" STREQUAL "CBC_FILE_LOC-NOTFOUND")
    message(STATUS "CBC: Found CBC_INCLUDE_DIRS `${CBC_INCLUDE_DIRS}`")
    break()
  endif()
endforeach(COIN_SUB_DIR)

if(NOT "${CBC_FILE_LOC}" STREQUAL "CBC_FILE_LOC-NOTFOUND")
  list(REMOVE_DUPLICATES CBC_INCLUDE_DIRS)
  unset(CBC_FIND_FILES)
  unset(CBC_FILE_LOC)

  if(WIN32 AND NOT UNIX)
    set(CBC_REQ_LIBS Osi OsiClp OsiCbc Clp Cgl Cbc CbcSolver CoinUtils)
  else()
    set(CBC_REQ_LIBS CbcSolver Cbc Cgl OsiClp Clp Osi CoinUtils)
  endif()

  foreach(CBC_LIB ${CBC_REQ_LIBS})
    set(CBC_LIB_LOC "CBC_LIB_LOC-NOTFOUND")
    find_library(CBC_LIB_LOC NAMES ${CBC_LIB} lib${CBC_LIB}
                PATH_SUFFIXES lib)
    if("${CBC_LIB_LOC}" STREQUAL "CBC_LIB_LOC-NOTFOUND")
      message(STATUS "CBC: Could not find library `${CBC_LIB}`")
      set(CBC_LIBRARY_PATHS "")
      break()
    endif()
    list(APPEND CBC_LIBRARY_PATHS ${CBC_LIB_LOC})
    add_library(${CBC_LIB} UNKNOWN IMPORTED)
    set_target_properties(${CBC_LIB} PROPERTIES
                          IMPORTED_LOCATION ${CBC_LIB_LOC}
                          INTERFACE_INCLUDE_DIRECTORIES "${CBC_INCLUDE_DIRS}")
    list(APPEND CBC_TARGETS ${CBC_LIB})
  endforeach(CBC_LIB)

  unset(CBC_REQ_LIBS)
  unset(CBC_LIB_LOC)

  if(UNIX AND NOT WIN32 AND NOT DEFINED EMSCRIPTEN)
    find_package(ZLIB)
    if(NOT ZLIB_FOUND)
      message(STATUS "CBC: Missing dependency `Zlib`")
      set(CBC_LIBRARY_PATHS "")
    else()
      list(APPEND CBC_LIBRARY_PATHS ${ZLIB_LIBRARIES})
      list(APPEND CBC_TARGETS ${ZLIB_LIBRARIES})
    endif()
  endif()
endif()
#]]

message(STATUS "CBC: CBC_INCLUDE_DIRS = `${CBC_INCLUDE_DIRS}`")
message(STATUS "CBC: CBC_LIBRARY_PATHS = `${CBC_LIBRARY_PATHS}`")
message(STATUS "CBC: CBC_TARGETS = `${CBC_TARGETS}`")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CBC_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CBC
  FOUND_VAR CBC_FOUND
  REQUIRED_VARS CBC_INCLUDE_DIRS CBC_LIBRARY_PATHS
  FAIL_MESSAGE "Could NOT find CBC, use CBC_ROOT to hint its location"
)

mark_as_advanced(CBC_INCLUDE_DIRS CBC_LIBRARY_PATHS)