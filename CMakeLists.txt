cmake_minimum_required (VERSION 3.13)

project(stream)

include(CheckLanguage)

check_language(C)
if(CMAKE_C_COMPILER)
    enable_language(C)
endif()
  
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
endif()
  
find_package(OpenMP COMPONENTS C Fortran)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(STREAM_UNIX_RELEASE_FLAGS "-O2")
set(STREAM_MSVC_RELEASE_FLAGS "/O2")

if(CMAKE_C_COMPILER)
    add_executable(stream_c stream.c)
    target_link_libraries(stream_c PRIVATE $<$<BOOL:${OpenMP_C_FOUND}>:OpenMP::OpenMP_C>)
    if(CMAKE_C_COMPILER_ID MATCHES "(GNU|Clang|Intel)")
        target_compile_options(stream_c PUBLIC "$<$<CONFIG:RELEASE>:${STREAM_UNIX_RELEASE_FLAGS}>")
    elseif(CMAKE_C_COMPILER_ID MATCHES "MSVC")
        target_compile_options(stream_c PUBLIC "$<$<CONFIG:RELEASE>:${STREAM_MSVC_RELEASE_FLAGS}>")
    endif()
endif()

if(CMAKE_Fortran_COMPILER)
    add_executable(stream_f stream.f mysecond.c)
    target_link_libraries(stream_f PRIVATE $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "(GNU|Clang|Intel)")
        target_compile_options(stream_f PUBLIC "$<$<CONFIG:RELEASE>:${STREAM_UNIX_RELEASE_FLAGS}>")
    endif()
endif()