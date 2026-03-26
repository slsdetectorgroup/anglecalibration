# TODO: check version
find_package(Python3 3.10 COMPONENTS Development NumPy REQUIRED)

set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
set(PYTHON_INCLUDE_DIR ${Python3_INCLUDE_DIRS})
set(PYTHON_LIBRARY ${Python3_LIBRARIES})

# python environment requires matplotlib 
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import matplotlib; print(matplotlib.__version__)"
    OUTPUT_VARIABLE MATPLOTLIB_VERSION
    RESULT_VARIABLE MATPLOTLIB_RESULT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT MATPLOTLIB_RESULT EQUAL 0)
    message(FATAL_ERROR "matplotlib Python package not found. Install with: conda install matplotlib")
else()
    message(STATUS "Found matplotlib: ${MATPLOTLIB_VERSION}")
endif()

# TODO: need to ensure python version larger than 3.10 for pybind11, maybe add check here?
message(STATUS "Found Python3: ${Python3_EXECUTABLE}, include dirs: ${Python3_INCLUDE_DIRS}, libraries: ${Python3_LIBRARIES}")

FetchContent_Declare(matplotlib_cpp
    GIT_REPOSITORY https://github.com/lava/matplotlib-cpp
    GIT_TAG master)
FetchContent_MakeAvailable(matplotlib_cpp)

if(TARGET matplotlib_cpp)
    message(STATUS "Matplotlib_cpp target is available: matplotlib_cpp")
else()
    message(FATAL_ERROR "Matplotlib_cpp target was not found!")
endif()

target_compile_options(matplotlib_cpp INTERFACE -w) # suppress warnings from matplotlibcpp which we cannot fix and which cause a lot of noise

target_sources(angle_calibration PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PlotCalibrationProcess.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PlotHelpers.cpp
)

target_include_directories(angle_calibration
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_include_directories(angle_calibration SYSTEM PUBLIC 
    $<BUILD_INTERFACE:${matplotlib_cpp_SOURCE_DIR}>
)

target_link_libraries(angle_calibration PUBLIC Python3::Python Python3::NumPy)





