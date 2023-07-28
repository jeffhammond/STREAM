# -------------------- Project Information ------------------- #

cmake_minimum_required(VERSION 3.25.1)

project(stream)

# ------------------- Environment Settings ------------------- #

if (APPLE)
    # Need llvm installed.
    set(CMAKE_CXX_COMPILER /usr/local/opt/llvm/bin/clang++)
    set(CMAKE_C_COMPILER /usr/local/opt/llvm/bin/clang)
    # Link llvm libraries.
    link_directories(/usr/local/opt/llvm/lib)
elseif(UNIX AND NOT APPLE)
    # Use g++ on Linux.
    set(CMAKE_CXX_COMPILER g++)
    set(CMAKE_C_COMPILER gcc)
else()
    # TODO: Windows support.
    message("Platform not support.")
endif()

# ------------------ Executable Compilation ------------------ #

file(
    GLOB
    SRC
    stream.c
)
# Build with OpenMP
add_executable(${PROJECT_NAME} ${SRC})
target_link_libraries(${PROJECT_NAME} PRIVATE -fopenmp)
target_compile_options(${PROJECT_NAME} PRIVATE -O2 -fopenmp)