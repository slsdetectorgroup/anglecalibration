# SPDX-License-Identifier: MPL-2.0
find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

set(TARGET_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include) # Doxyfile uses this variable to find the include directory

# Doxygen
set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

# Sphinx
set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(
  docs
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
  COMMAND
    ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_BINARY_DIR}" # find angle_calibration module
    ${SPHINX_EXECUTABLE} -W -a -b html -Dversion=${PROJECT_VERSION}
    -Drelease=${PROJECT_VERSION}
    -Dbreathe_projects.angle_calibration=${CMAKE_CURRENT_BINARY_DIR}/xml -c
    "${CMAKE_CURRENT_SOURCE_DIR}" ${SPHINX_SOURCE} ${SPHINX_BUILD}/html
  COMMENT "Generating documentation with Sphinx")