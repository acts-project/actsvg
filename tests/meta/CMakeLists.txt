# Set up the core tests.
actsvg_add_test( meta
  "clusters.cpp"
  "detector.cpp"
  "eta_lines.cpp"
  "grid.cpp"
  "wire_chamber.cpp"
  "portal.cpp"
  "seeds.cpp"
  "surface_materials.cpp"
  "surfaces.cpp"
  "trajectories.cpp"
  "volume.cpp"
LINK_LIBRARIES GTest::gtest_main actsvg::core actsvg::meta actsvg::data
)
