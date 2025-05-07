configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/setup.sh.in
    ${PROJECT_BINARY_DIR}/this_actsvg.sh
    @ONLY
)
install(
    FILES ${PROJECT_BINARY_DIR}/this_actsvg.sh
    DESTINATION ${CMAKE_INSTALL_BINDIR}
)