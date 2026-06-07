add_executable(gauss_transform_test gauss_transform_test.cc)
target_include_directories(gauss_transform_test PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${VXL_VNL_INCLUDE_DIR}
)
target_link_libraries(gauss_transform_test
    GTest::gtest_main
    Threads::Threads
    vnl
)
gtest_discover_tests(gauss_transform_test)

add_executable(fgt_utils_test fgt_utils_test.cc)
target_include_directories(fgt_utils_test PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${VXL_VNL_INCLUDE_DIR}
)
target_link_libraries(fgt_utils_test
    GTest::gtest_main
    Threads::Threads
    vnl
)
gtest_discover_tests(fgt_utils_test)

add_executable(em_utils_test em_utils_test.cc)
target_include_directories(em_utils_test PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${VXL_VNL_INCLUDE_DIR}
)
target_link_libraries(em_utils_test
    GTest::gtest_main
    Threads::Threads
    vnl
)
gtest_discover_tests(em_utils_test)

add_executable(rotation_utils_test rotation_utils_test.cc)
target_include_directories(rotation_utils_test PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${VXL_VNL_INCLUDE_DIR}
)
target_link_libraries(rotation_utils_test
    GTest::gtest_main
    Threads::Threads
    vnl
)
gtest_discover_tests(rotation_utils_test)

add_executable(rbf_utils_test rbf_utils_test.cc)
target_include_directories(rbf_utils_test PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${VXL_VNL_INCLUDE_DIR}
)
target_link_libraries(rbf_utils_test
    GTest::gtest_main
    Threads::Threads
    vnl
)
gtest_discover_tests(rbf_utils_test)

add_executable(normalization_utils_test normalization_utils_test.cc)
target_include_directories(normalization_utils_test PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${VXL_VNL_INCLUDE_DIR}
)
target_link_libraries(normalization_utils_test
    GTest::gtest_main
    Threads::Threads
    vnl
)
gtest_discover_tests(normalization_utils_test)

add_executable(downsample_utils_test downsample_utils_test.cc)
target_include_directories(downsample_utils_test PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${VXL_VNL_INCLUDE_DIR}
)
target_link_libraries(downsample_utils_test
    GTest::gtest_main
    Threads::Threads
    vnl
)
gtest_discover_tests(downsample_utils_test)

add_executable(match_utils_test match_utils_test.cc)
target_include_directories(match_utils_test PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${VXL_VNL_INCLUDE_DIR}
)
target_link_libraries(match_utils_test
    GTest::gtest_main
    Threads::Threads
    vnl
)
gtest_discover_tests(match_utils_test)

# misc_utils.cc is a regular compiled file (not template header-only),
# so it is included directly as a source rather than linked via a library.
add_executable(misc_utils_test misc_utils_test.cc misc_utils.cc)
target_include_directories(misc_utils_test PRIVATE
    ${CMAKE_SOURCE_DIR}
)
target_link_libraries(misc_utils_test
    GTest::gtest_main
    Threads::Threads
)
gtest_discover_tests(misc_utils_test)

add_executable(io_utils_test io_utils_test.cc)
target_include_directories(io_utils_test PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${VXL_VNL_INCLUDE_DIR}
)
target_link_libraries(io_utils_test
    GTest::gtest_main
    Threads::Threads
    vnl
)
gtest_discover_tests(io_utils_test)
