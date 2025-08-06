
if(NOT TARGET ${TARGET_NAME})
    message(FATAL_ERROR "Target '${TARGET_NAME}' not found")
endif()

FetchContent_Declare(
    gnuplot_iostream
    GIT_REPOSITORY https://github.com/dstahlke/gnuplot-iostream.git
    GIT_TAG        master
)


FetchContent_MakeAvailable(gnuplot_iostream)

find_package(Boost REQUIRED COMPONENTS iostreams)

list(APPEND PUBLICHEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/plot_histogram.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PlotHelpers.hpp
)

#add plot_histogram.hpp
target_include_directories(angle_calibration PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

#mmh actually dont know if i have to install it or where to install headers during install process
target_include_directories(angle_calibration PUBLIC $<BUILD_INTERFACE:${gnuplot_iostream_SOURCE_DIR}> 
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")


target_include_directories(angle_calibration PUBLIC ${Boost_INCLUDE_DIRS})

target_link_libraries(angle_calibration PUBLIC ${Boost_LIBRARIES})

target_compile_definitions(angle_calibration PUBLIC ANGCAL_PLOT)




