cmake_minimum_required( VERSION 3.31 )
project( oktal LANGUAGES CXX )

set( CMAKE_CXX_STANDARD 23 )
set( CMAKE_CXX_STANDARD_REQUIRED ON ) # Added ON here
set( CMAKE_COMPILE_WARNING_AS_ERROR ON )
add_compile_options(-Wall -Wextra -Wconversion)
option( ADVPT_ENABLE_HTGFILE "" ON )
list(APPEND CMAKE_PREFIX_PATH "/proj/cipdata/da15siwa/advpt-packages")


# 1. Pull Helper Library FIRST so the targets exist
message( STATUS "Downloading AdvPT helper libraries..." )
include( FetchContent )
FetchContent_Declare(
    advpt-helpers
    URL https://teaching.pages.i10git.cs.fau.de/advpt-exercise-ws2526/assets/advpt-helpers.tar.gz
)
FetchContent_MakeAvailable( advpt-helpers )


# 2. Create the library
add_library( oktal )
target_include_directories( oktal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include )

# 3. LINK the helper library to oktal
# This automatically fixes the "No such file or directory" error for VtkHtgFile.hpp
target_link_libraries( oktal PUBLIC advpt_htgfile )

# 4. Now add the subdirectories
add_subdirectory( src )

# 5. Add applications
add_subdirectory( apps )

# 6. Add test suite
option( OKTAL_BUILD_TESTSUITE "Build the Oktal test suite" ON )
if( OKTAL_BUILD_TESTSUITE )
    include( CTest )
    add_subdirectory( tests )
endif()



