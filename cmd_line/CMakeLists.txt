add_executable(cmd_line cmd_line.cc)
target_link_libraries(cmd_line PRIVATE voro++)
#cannot have two targets with the same name, so renaming cmd_line
set_target_properties(cmd_line PROPERTIES OUTPUT_NAME voro++ 
RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/src") 
install(TARGETS cmd_line RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
