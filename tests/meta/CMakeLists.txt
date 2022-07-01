# Set up the core tests.
actsvg_add_test( meta
  "clusters.cpp"
  "detector.cpp"
  "volume.cpp"
  "grid.cpp"
  "barrel_sheet.cpp"
  "endcap_sheet.cpp"
  "surface_sheet.cpp"
LINK_LIBRARIES GTest::gtest_main actsvg::core actsvg::meta actsvg::data)
