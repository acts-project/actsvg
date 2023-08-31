# Set up the core tests.
actsvg_add_test( web
  "web_builder.cpp"
LINK_LIBRARIES GTest::gtest_main actsvg::core actsvg::web)
