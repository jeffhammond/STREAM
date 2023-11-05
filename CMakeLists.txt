cmake_minimum_required(VERSION 3.20)

project(
  STREAM
  VERSION 1.0
  DESCRIPTION "STREAM benchmark"
  LANGUAGES C Fortran)

# The C STREAM benchmark only requires
# stream.c
add_executable(stream_c.exe stream.c)

# The Fortran STREAM benchmark requires
# stream.f and mysecond.o from mysecond.c
add_executable(stream_f.exe stream.f mysecond.c)

# Look for OpenMP support is found, link it to the executables
# Note that if you are using clang on macOS, you will need to
# install libomp via Homebrew and then set the following
# environment variables:
#   export OpenMP_ROOT=$(brew --prefix)/opt/libomp
# see https://www.scivision.dev/cmake-openmp/ for more details

find_package(OpenMP COMPONENTS C Fortran)
# Use a generator expression to link OpenMP to the executables
target_link_libraries(
  stream_c.exe 
  PRIVATE
  $<$<BOOL:${OpenMP_C_FOUND}>:OpenMP::OpenMP_C>
)

target_link_libraries(
  stream_f.exe 
  PRIVATE
  $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>
)

# Per @scivision, ignore the build directory
if(NOT PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  # Git auto-ignore out-of-source build directory
  file(GENERATE OUTPUT .gitignore CONTENT "*")
endif()
