# Set up the core tests.
actsvg_add_test( meta
  "barrel_sheet.cpp"
  "clusters.cpp"
  "detector.cpp"
  "endcap_sheet.cpp"
  "eta_lines.cpp"
  "grid.cpp"
  "portal.cpp"
  "surface_sheet.cpp"
  "surfaces.cpp"
  "volume.cpp"
LINK_LIBRARIES GTest::gtest_main actsvg::core actsvg::meta actsvg::data)
