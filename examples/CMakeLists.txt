add_executable(example_calibration)

target_link_libraries(example_calibration PRIVATE angle_calibration compiler_flags)
set_target_properties(example_calibration PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
        OUTPUT_NAME run_calibration
)

target_sources(example_calibration PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/example_calibration.cpp )

add_executable(example_conversion)

target_link_libraries(example_conversion PRIVATE angle_calibration compiler_flags)
set_target_properties(example_conversion PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
        OUTPUT_NAME run_conversion
)

target_sources(example_conversion PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/example_conversion.cpp )


