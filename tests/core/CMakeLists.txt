# Set up the core tests.
actsvg_add_test( core
  "barrel.cpp" 
  "endcap.cpp"
  "polygon.cpp"
  "style.cpp"
  "svg.cpp"
  "text.cpp"
  "views.cpp"
LINK_LIBRARIES GTest::gtest_main actsvg::core )