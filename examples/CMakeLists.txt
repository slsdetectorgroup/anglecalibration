add_executable(example)
target_link_libraries(example PRIVATE angle_calibration compiler_flags)
set_target_properties(example PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
        OUTPUT_NAME run_example
)

target_sources(example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/example_calibration.cpp )

