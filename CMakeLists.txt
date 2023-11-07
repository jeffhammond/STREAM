cmake_minimum_required(VERSION 3.15)

project(
  STREAM
  VERSION 1.0
  DESCRIPTION "STREAM benchmark"
  LANGUAGES C)

enable_testing()

option(fortran "build fortran version" ON)

if(fortran)
  enable_language(Fortran)
endif()

if(CMAKE_C_COMPILER_ID MATCHES "(GNU|Clang)")
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:-O3;-march=native;-Wall>")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:-O3;-Wall>")
  if(WIN32)
    add_compile_options("$<$<COMPILE_LANGUAGE:C>:/QxHost>")
  else()
    add_compile_options("$<$<COMPILE_LANGUAGE:C>:-xHost>")
  endif()
elseif(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:/W3>")
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-O3;-march=native;-Werror=line-truncation;-Wall>")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-O3;-warn>")
  if(WIN32)
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:/QxHost>")
  else()
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-xHost>")
  endif()
endif()

# Look for OpenMP support is found, link it to the executables
# Note that if you are using clang on macOS, you will need to
# install libomp via Homebrew and then set the following
# environment variables:
#   export OpenMP_ROOT=$(brew --prefix)/opt/libomp
# see https://www.scivision.dev/cmake-openmp/ for more details

find_package(OpenMP COMPONENTS C Fortran)

# --- C stream_c

add_executable(stream_c stream.c)
target_link_libraries(stream_c PRIVATE $<$<BOOL:${OpenMP_C_FOUND}>:OpenMP::OpenMP_C>)

add_test(NAME STREAM_C COMMAND stream_c)

# --- Fortran stream_f

if(fortran)
add_executable(stream_f stream.f)
target_link_libraries(stream_f PRIVATE $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>)

add_test(NAME STREAM_Fortran COMMAND stream_f)
endif()

# --- ignore build directory
if(NOT PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  file(GENERATE OUTPUT .gitignore CONTENT "*")
endif()
