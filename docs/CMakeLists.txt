find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

set(TARGET_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

#Doxygen
set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/DoxyFile.in)
set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/DoxyFile)
configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

#Sphinx
set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR})


file(GLOB SPHINX_SOURCE_FILES CONFIGURE_DEPENDS "${SPHINX_SOURCE}/*.rst")

message("SPHINX_SOURCE_FILES: ${SPHINX_SOURCE_FILES}")


if(false)
foreach(filename ${SPHINX_SOURCE_FILES})
    get_filename_component(fname ${filename} NAME)
    message(STATUS "Copying ${filename} to ${SPHINX_BUILD}/src/${fname}")
    configure_file(${filename} "${SPHINX_BUILD}/src/${fname}")  # copy file to build directory
endforeach(filename ${SPHINX_SOURCE_FILES})
endif()


configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
    "${SPHINX_BUILD}/conf.py"
    @ONLY
) 

#creates c++ api using doxygen and breathe 
add_custom_target(
    docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT} #create xml 
    COMMAND ${SPHINX_EXECUTABLE} -a -b html
    -Dbreathe_projects.angle_calibration=${CMAKE_CURRENT_BINARY_DIR}/xml
    -c "${SPHINX_BUILD}"
    ${SPHINX_SOURCE}
    ${SPHINX_BUILD}/html
    COMMENT "Generating documentation with Sphinx"
)

#only creates html using sphinx - no c++ api 
add_custom_target(
    rst
    COMMAND ${SPHINX_EXECUTABLE} -a -b html
    -c "${SPHINX_BUILD}"
    ${SPHINX_SOURCE}
    ${SPHINX_BUILD}/html
    COMMENT "Generating documentation with Sphinx"
)

