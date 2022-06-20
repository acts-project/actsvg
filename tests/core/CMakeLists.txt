# Set up the core tests.
actsvg_add_test( core
  "arc.cpp"
  "barrel.cpp" 
  "circle.cpp"
  #"endcap.cpp"
  "grids.cpp"
  "line.cpp"
  "polygon.cpp"
  "style.cpp"
  "svg.cpp"
  "text.cpp"
  "views.cpp"
LINK_LIBRARIES GTest::gtest_main actsvg::core actsvg::data)