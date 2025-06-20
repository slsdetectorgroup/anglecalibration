# Download catch2 if configured to do so
if (ANGCAL_FETCH_CATCH)
    FetchContent_Declare(
        Catch2
        GIT_SHALLOW    TRUE
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.5.3
    )
    FetchContent_MakeAvailable(Catch2)
else()
    # Otherwise look for installed catch2
    find_package(Catch2 3 REQUIRED)
endif()

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)

add_executable(tests)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain angle_calibration compiler_flags)
set_target_properties(tests PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
        OUTPUT_NAME run_tests
)

include(CTest)
include(Catch)
catch_discover_tests(tests)

target_sources(tests PRIVATE ${TestSources} )

#configure a header to pass test file paths
get_filename_component(TEST_FILE_PATH ${PROJECT_SOURCE_DIR}/data ABSOLUTE)
configure_file(test_config.hpp.in test_config.hpp)
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

