cmake_minimum_required(VERSION 3.20)

project(Cminor LANGUAGES Fortran C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# --- ISSA (optional private submodule) ---
set(ISSA_SRC_DIR ${CMAKE_SOURCE_DIR}/EXTERNALS/cminor-issa/SRC)
if(EXISTS ${ISSA_SRC_DIR}/issa_reduce.f90)
  option(CMINOR_ENABLE_ISSA "Build with cminor-issa reduction (-DISSA)" ON)
else()
  option(CMINOR_ENABLE_ISSA "Build with cminor-issa reduction (-DISSA)" OFF)
endif()
if(CMINOR_ENABLE_ISSA AND NOT EXISTS ${ISSA_SRC_DIR}/issa_reduce.f90)
  message(FATAL_ERROR "CMINOR_ENABLE_ISSA=ON but ${ISSA_SRC_DIR}/issa_reduce.f90 not found")
endif()

# --- library sources (no PROGRAM) ---
set(CMINOR_LIB_SOURCES
  SRC/UniRnk_Mod.f90
  SRC/LexicalStringSort.f90
  SRC/Kind_Mod.f90
  SRC/Control_Mod.f90
  SRC/Reac_Mod.f90
  SRC/Meteo_Mod.f90
  SRC/InitRoutines_Mod.f90
  SRC/NetCDF_Mod.f90
  SRC/Sparse_Mod.f90
  SRC/String_Mod.f90
  SRC/HashStr_Mod.f90
  SRC/InputTool_Mod.f90
  SRC/ChemSys_Mod.f90
  SRC/IO_Mod.f90
  SRC/ChemKinInput_Mod.f90
  SRC/fp_parameters.f90
  SRC/fparser.f90
  SRC/Rates_Mod.f90
  SRC/Rosenbrock_Mod.f90
  SRC/Integration_Mod.f90
  SRC/Cminor_host_Mod.f90
  SRC/Cminor_Driver_Mod.f90
)

if(CMINOR_ENABLE_ISSA)
  list(APPEND CMINOR_LIB_SOURCES
    ${ISSA_SRC_DIR}/issa_flux_io.f90
    ${ISSA_SRC_DIR}/issa_graph_tools.f90
    ${ISSA_SRC_DIR}/issa_reduce.f90
  )
endif()

add_library(cminor_lib STATIC ${CMINOR_LIB_SOURCES})
add_executable(cminor SRC/Cminor.f90)
target_link_libraries(cminor PRIVATE cminor_lib)

# Module search path for .mod files and METHODS/
target_include_directories(cminor_lib PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
  ${CMAKE_SOURCE_DIR}/METHODS
  ${CMAKE_SOURCE_DIR}
)
target_include_directories(cminor PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Free-form + cpp (matches Makefile)
target_compile_options(cminor_lib PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:-ffree-form;-ffixed-line-length-none;-ffree-line-length-none;-cpp>
)
target_compile_options(cminor PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:-ffree-form;-ffixed-line-length-none;-ffree-line-length-none;-cpp>
)

if(CMINOR_ENABLE_ISSA)
  target_compile_definitions(cminor_lib PUBLIC ISSA)
  target_compile_definitions(cminor PRIVATE ISSA)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(cminor_lib PRIVATE
    $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-O3>
    $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-O0;-g;-fcheck=all;-fbacktrace;-Wall>
  )
  target_compile_options(cminor PRIVATE
    $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-O3>
    $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-O0;-g>
  )
endif()

# --- NetCDF (mirror Makefile: Fortran + C prefixes; nf-config alone misses -L for libnetcdf) ---
find_program(NF_CONFIG nf-config)
find_program(NC_CONFIG nc-config)
set(NETCDF_DIR "$ENV{NETCDF_DIR}")
set(NETCDF_C_DIR "$ENV{NETCDF_C_DIR}")
set(HDF5_DIR "$ENV{HDF5_DIR}")
if(NOT NETCDF_DIR)
  if(APPLE AND EXISTS /opt/homebrew/opt/netcdf-fortran)
    set(NETCDF_DIR /opt/homebrew/opt/netcdf-fortran)
  elseif(EXISTS /usr/include)
    set(NETCDF_DIR /usr)
  endif()
endif()
if(NOT NETCDF_C_DIR)
  if(APPLE AND EXISTS /opt/homebrew/opt/netcdf)
    set(NETCDF_C_DIR /opt/homebrew/opt/netcdf)
  else()
    set(NETCDF_C_DIR ${NETCDF_DIR})
  endif()
endif()
if(NOT HDF5_DIR)
  if(APPLE AND EXISTS /opt/homebrew/opt/hdf5)
    set(HDF5_DIR /opt/homebrew/opt/hdf5)
  else()
    set(HDF5_DIR /usr)
  endif()
endif()

if(NF_CONFIG)
  execute_process(COMMAND ${NF_CONFIG} --fflags OUTPUT_VARIABLE NETCDF_FFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  separate_arguments(NETCDF_FFLAGS_LIST NATIVE_COMMAND "${NETCDF_FFLAGS}")
  target_compile_options(cminor_lib PUBLIC ${NETCDF_FFLAGS_LIST})
endif()
if(NETCDF_DIR)
  target_include_directories(cminor_lib PUBLIC ${NETCDF_DIR}/include)
endif()

# Explicit -L paths (nf-config --flibs only lists netcdf-fortran's lib dir)
target_link_directories(cminor_lib PUBLIC
  ${NETCDF_DIR}/lib
  ${NETCDF_C_DIR}/lib
  ${HDF5_DIR}/lib
)
if(APPLE)
  target_link_libraries(cminor_lib PUBLIC netcdff netcdf hdf5_hl hdf5 curl)
  target_link_options(cminor_lib PUBLIC
    "LINKER:-rpath,${NETCDF_DIR}/lib"
    "LINKER:-rpath,${NETCDF_C_DIR}/lib"
    "LINKER:-rpath,${HDF5_DIR}/lib"
  )
else()
  # Linux package names vary (hdf5 vs hdf5_serial); prefer nc-config libs when present
  if(NC_CONFIG)
    execute_process(COMMAND ${NC_CONFIG} --libs OUTPUT_VARIABLE NETCDF_CLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
    separate_arguments(NETCDF_CLIBS_LIST NATIVE_COMMAND "${NETCDF_CLIBS}")
    target_link_libraries(cminor_lib PUBLIC netcdff ${NETCDF_CLIBS_LIST})
  else()
    target_link_libraries(cminor_lib PUBLIC netcdff netcdf hdf5_hl hdf5 curl)
  endif()
endif()

# --- BLAS/LAPACK ---
if(APPLE)
  target_link_options(cminor_lib PUBLIC "LINKER:-framework,Accelerate")
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  target_link_libraries(cminor_lib PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

# --- optional issa_reduce standalone ---
if(CMINOR_ENABLE_ISSA)
  add_executable(issa_reduce EXTERNALS/cminor-issa/issa_reduce_main.f90)
  target_link_libraries(issa_reduce PRIVATE cminor_lib)
  target_compile_options(issa_reduce PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:-ffree-form;-ffixed-line-length-none;-ffree-line-length-none;-cpp>
  )
  target_compile_definitions(issa_reduce PRIVATE ISSA)
endif()

# --- install (library for ICON coupling later) ---
install(TARGETS cminor_lib cminor
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
  DESTINATION include/cminor
  FILES_MATCHING PATTERN "*.mod"
)

message(STATUS "Cminor: ISSA=${CMINOR_ENABLE_ISSA}  build_type=${CMAKE_BUILD_TYPE}")
