cmake_minimum_required(VERSION 3.15)

project(
  STREAM
  VERSION 1.0
  DESCRIPTION "STREAM benchmark"
  LANGUAGES C Fortran)

enable_testing()

if(CMAKE_C_COMPILER_ID MATCHES "(GNU|Clang)")
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:-O3;-march=native>")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:-O3;-xHost>")
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-O3;-march=native;-Werror=line-truncation>")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-O3;-xHost>")
endif()

# The C STREAM benchmark only requires
# stream.c
add_executable(stream_c stream.c)

add_test(NAME STREAM_C COMMAND stream_c)

# The Fortran STREAM benchmark requires
# stream.f and mysecond.o from mysecond.c
add_executable(stream_f stream.f mysecond.c)

add_test(NAME STREAM_Fortran COMMAND stream_f)

# Look for OpenMP support is found, link it to the executables
# Note that if you are using clang on macOS, you will need to
# install libomp via Homebrew and then set the following
# environment variables:
#   export OpenMP_ROOT=$(brew --prefix)/opt/libomp
# see https://www.scivision.dev/cmake-openmp/ for more details

find_package(OpenMP COMPONENTS C Fortran)
target_link_libraries(stream_c PRIVATE $<$<BOOL:${OpenMP_C_FOUND}>:OpenMP::OpenMP_C>)
target_link_libraries(stream_f PRIVATE $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>)

# Per @scivision, ignore the build directory
# (see https://www.scivision.dev/cmake-auto-gitignore-build-dir/)
if(NOT PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  file(GENERATE OUTPUT .gitignore CONTENT "*")
endif()
