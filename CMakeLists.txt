cmake_minimum_required(VERSION 3.15)

project(angle_calibration 
DESCRIPTION "angle calibration library for PSI mythen detector"
LANGUAGES  CXX)

# Read VERSION file into project version
set(VERSION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/VERSION")
file(READ "${VERSION_FILE}" VERSION_CONTENT)
string(STRIP "${VERSION_CONTENT}" PROJECT_VERSION_STRING)
set(PROJECT_VERSION ${PROJECT_VERSION_STRING})
message(STATUS "PROJECT_VERSION = ${PROJECT_VERSION}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#Set default build type if none was specified
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (default Release)" FORCE)
endif()

#maybe add c++ 17 compiler tests
#set_std_fs_lib()

# Set lower / upper case project names
string(TOUPPER "${PROJECT_NAME}" PROJECT_NAME_UPPER)
string(TOLOWER "${PROJECT_NAME}" PROJECT_NAME_LOWER)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set targets export name 
set(TARGETS_EXPORT_NAME "${PROJECT_NAME_LOWER}-targets")

include(FetchContent)

option(ANGCAL_PYTHON_BINDINGS "Build python bindings" OFF)
option(ANGCAL_TESTS "Build tests" OFF)
option(ANGCAL_DOCS "Build documentation" OFF)
option(ANGCAL_PLOT "plot for visualization" OFF)
option(ANGCAL_VERBOSE "print debug info" OFF)
option(ANGCAL_INSTALL_PYTHONEXT "Install the python extension in the install tree under CMAKE_INSTALL_PREFIX/angle_calibration/" OFF)
option(ASAN "Enable AddressSanitizer" OFF)
option(ANGCAL_FETCH_HDF5 "Use FetchContent to download hdf5-devel" OFF)
option(ANGCAL_FETCH_AARE "Use FetchContent to download aare from github" OFF)
option(ANGCAL_FETCH_CATCH "Use FetchContent to download catch2" ON)
option(ANGCAL_FETCH_PYBIND11 "Use FetchContent to download pybind11" ON)


if(ANGCAL_FETCH_AARE)
    FetchContent_Declare(aare GIT_REPOSITORY https://github.com/slsdetectorgroup/aare
    GIT_TAG main) #how to use newest version 
    FetchContent_MakeAvailable(aare)
    install(
        TARGETS aare_core
        EXPORT ${TARGETS_EXPORT_NAME}
    )
    install(TARGETS aare_compiler_flags
    EXPORT ${TARGETS_EXPORT_NAME}
    ) #mmh is this the way to go I think I will define the same compiler options twice now? can directly use aare_compiler_flags instead of compiler_flags
    message(STATUS "target: aare")
else()
    message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
    find_package(aare REQUIRED) 
    message(STATUS "found aare: ${AARE_INCLUDE_DIRS}")

endif()

if(ANGCAL_FETCH_HDF5)
    message(FATAL_ERROR "Fetching HDF5 via FetchContent is not supported here. Please install it via your system. 
    For Ubuntu: sudo apt install libhdf5-dev
    For Red Hat: sudo dnf install hdf5-devel
    For MacOS: brew install hdf5")
else()
    find_package(HDF5 QUIET COMPONENTS CXX)
    if (HDF5_FOUND)
        message(STATUS "Found HDF5: ${HDF5_INCLUDE_DIRS}")
    else()
        message(FATAL_ERROR "HDF5 was NOT found! Please install it via your system  
            For Ubuntu: sudo apt install libhdf5-dev
            For Red Hat: sudo dnf install hdf5-devel
            For MacOS: brew install hdf5")
    endif()
endif()



add_library(compiler_flags INTERFACE) #INTERFACE 
target_compile_features(compiler_flags INTERFACE cxx_std_17)

if(ANGCAL_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

if(MSVC) #Microsoft compiler
    target_compile_options(
        compiler_flags 
        INTERFACE
        /w # disable warnings
    ) #why no warnings?
else() #gcc clang compiler options
    target_compile_options(
        compiler_flags 
        INTERFACE 
            -Wall 
            -Wextra 
            -pedantic 
            -Wshadow 
            -Wold-style-cast
            -Wnon-virtual-dtor
            -Woverloaded-virtual
            -Wdouble-promotion
            -Wformat=2
            -Wredundant-decls
            -Wvla
            -Wdouble-promotion
            -Werror=return-type #important can cause segfault in optimzed build 
    )
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Release build")
    if(MSVC)
        target_compile_options(compiler_flags INTERFACE /O2)
    else()
        target_compile_options(compiler_flags INTERFACE -O3)
    endif()
else()
    message(STATUS "Debug build")
    if(MSVC)
        target_compile_options(
            compiler_flags 
            INTERFACE 
                /Od 
                /Zi 
                /MDd
                /D_ITERATOR_DEBUG_LEVEL=2
        )
        target_link_options(
            compiler_flags
            INTERFACE
                /DEBUG:FULL
        )
    endif()
endif()

if(ASAN)
    message(STATUS "AddressSanitizer enabled")
    target_compile_options(
        compiler_flags 
        INTERFACE 
            -fsanitize=address,undefined,pointer-compare 
            -fno-omit-frame-pointer
    )
    target_link_options(compiler_flags 
    INTERFACE 
        -fsanitize=address,,undefined,pointer-compare 
            -fno-omit-frame-pointer
    )
endif()

set(PUBLICHEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/AngleCalibration.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/CustomFiles.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Flatfield.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Hdf5FileReader.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/MythenFileReader.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/MythenDetectorSpecifications.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/helpers/FileInterface.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/logger.hpp
) 

set(SourceFiles 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/AngleCalibration.cpp
)



add_library(angle_calibration STATIC ${SourceFiles}) 


target_include_directories(angle_calibration PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>" 
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

target_link_libraries(angle_calibration PUBLIC aare_core HDF5::HDF5 PRIVATE compiler_flags) 

get_property(allTargets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
message(STATUS "Available targets: ${allTargets}")

if(ANGCAL_PLOT)
set(TARGET_NAME angle_calibration)
add_subdirectory(plotting_scripts)
endif()

add_subdirectory(examples)


if(ANGCAL_TESTS)
    enable_testing()
    add_subdirectory(tests)
    set(TestSources 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/AngleCalibration.test.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Hdf5FileReader.test.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/MythenFileReader.test.cpp)
    target_sources(tests PRIVATE ${TestSources})
endif()

if(ANGCAL_PYTHON_BINDINGS)
    set_property(TARGET angle_calibration PROPERTY POSITION_INDEPENDENT_CODE ON) #what is this doing
endif()

if(ANGCAL_DOCS)
    add_subdirectory(docs)
endif()

if(ANGCAL_VERBOSE)
    target_compile_definitions(angle_calibration PUBLIC ANGCAL_VERBOSE)
    target_compile_definitions(angle_calibration PUBLIC ANGCAL_LOG_LEVEL=angcal::logDEBUG)
else()
    target_compile_definitions(angle_calibration PUBLIC ANGCAL_LOG_LEVEL=angcal::logINFO)
endif()

set_target_properties(angle_calibration PROPERTIES
    PUBLIC_HEADER "${PUBLICHEADERS}"
)

#only build library if not added in other cmake with add_subdirectory 
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
install(TARGETS angle_calibration compiler_flags EXPORT "${TARGETS_EXPORT_NAME}" 
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/angle_calibration)
        set(CMAKE_INSTALL_DIR "share/cmake/${PROJECT_NAME}")
        set(PROJECT_LIBRARIES angle_calibration compiler_flags )
        include(cmake/package_config.cmake)
endif()

if(NOT ANGCAL_FETCH_AARE)
add_custom_target(
    check-format
    COMMAND find \( -name "*.cpp" -o -name "*.hpp" \) -not -path "./build/*" | xargs -I {} -n 1 -P 10 bash -c  "clang-format -Werror -style=\"file:.clang-format\" {} | diff {} -"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Checking code formatting with clang-format"
    VERBATIM
    
)

add_custom_target(
    format-files
    COMMAND find \( -name "*.cpp" -o -name "*.hpp" \) -not -path "./build/*" | xargs -I {} -n 1 -P 10 bash -c "clang-format -i -style=\"file:.clang-format\" {}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Formatting with clang-format"
    VERBATIM
)

add_custom_target(
    clang-tidy
    COMMAND find \( -path "./src/*" -a -not -path "./src/python/*" -a \( -name "*.cpp"   -not -name "*.test.cpp" \)  \)    -not -name "CircularFifo.hpp" -not -name "ProducerConsumerQueue.hpp" -not -name "VariableSizeClusterFinder.hpp" | xargs -I {} -n 1 -P 10 bash -c "${CLANG_TIDY_COMMAND} --config-file=.clang-tidy -p build {}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "linting with clang-tidy"
    VERBATIM
)
endif()


