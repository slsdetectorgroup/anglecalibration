
if(NOT TARGET ${TARGET_NAME})
    message(FATAL_ERROR "Target '${TARGET_NAME}' not found")
endif()

FetchContent_Declare(
    gnuplot_iostream
    GIT_REPOSITORY https://github.com/dstahlke/gnuplot-iostream.git
    GIT_TAG        master
)

FetchContent_MakeAvailable(gnuplot_iostream)

if(ANGCAL_FETCH_BOOST)

# doesnt work !!! 
# a lot is fetched cant only fetch boost_iostream.git too many dependencies 
FetchContent_Declare(
  boost
  GIT_REPOSITORY https://github.com/boostorg/boost.git
  GIT_TAG master
)

set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BOOST_INCLUDE_LIBRARIES "iostreams" CACHE STRING "" FORCE) 
set(BOOST_EXCLUDE_LIBRARIES "mpi,python" CACHE STRING "" FORCE)

FetchContent_MakeAvailable(boost)

if(NOT boost_POPULATED)
    FetchContent_Populate(boost)
    add_subdirectory(${boost_SOURCE_DIR} ${boost_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

install(TARGETS boost
    EXPORT ${project}-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

else() 
    # sudo dnf install boost-devel
    find_package(Boost REQUIRED COMPONENTS iostreams)
endif() 

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


target_link_libraries(angle_calibration PUBLIC Boost::iostreams)
#target_link_libraries(angle_calibration PUBLIC ${Boost_LIBRARIES})


target_compile_definitions(angle_calibration PUBLIC ANGCAL_PLOT)




