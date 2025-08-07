find_package (Python 3.10 COMPONENTS Interpreter Development.Module REQUIRED)
set(PYBIND11_FINDPYTHON ON) # Needed for RH8

# Download or find pybind11 depending on configuration
if(ANGCAL_FETCH_PYBIND11)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG        v2.13.6
    )
    FetchContent_MakeAvailable(pybind11)
else()
    find_package(pybind11 2.13 REQUIRED)
endif()

# Add the compiled python extension
pybind11_add_module(
    _angcal           # name of the module
    src/module.cpp  # source file
)

message("cmake binary directory: ${CMAKE_BINARY_DIR}")

set_target_properties(_angcal PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/angcal 
)

target_link_libraries(_angcal PRIVATE angle_calibration compiler_flags)


# List of python files to be copied to the build directory
set( PYTHON_FILES
    angcal/__init__.py
)

# Copy the python files to the build directory
foreach(FILE ${PYTHON_FILES})
    configure_file(${FILE} ${CMAKE_BINARY_DIR}/${FILE})
endforeach(FILE ${PYTHON_FILES})

if(ANGCAL_INSTALL_PYTHONEXT)
    install(
        TARGETS _angcal
        EXPORT "${TARGETS_EXPORT_NAME}"
        LIBRARY DESTINATION angle_calibration
        COMPONENT python
    )
    
    install(
        FILES ${PYTHON_FILES} 
        DESTINATION angle_calibration
        COMPONENT python
        )
endif()