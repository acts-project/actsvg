# Set up the core tests.
actsvg_add_test( meta
  "barrel_sheet.cpp"
  "clusters.cpp"
  "detector.cpp"
  "endcap_sheet.cpp"
  "eta_lines.cpp"
  "grid.cpp"
  "wire_chamber.cpp"
  "portal.cpp"
  "seeds.cpp"
  "surface_sheet.cpp"
  "surface_materials.cpp"
  "surfaces.cpp"
  "trajectories.cpp"
  "volume.cpp"
LINK_LIBRARIES GTest::gtest_main actsvg::core actsvg::meta actsvg::data)
