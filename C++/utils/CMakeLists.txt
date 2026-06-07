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
