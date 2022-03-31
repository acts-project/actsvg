# Set up the core tests.
actsvg_add_test( meta
  "detector.cpp" 
  "volume.cpp"
  "grid.cpp"
LINK_LIBRARIES GTest::gtest_main actsvg::core actsvg::meta actsvg::data)
