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
