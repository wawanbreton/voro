file(GLOB EXAMPLE_SOURCES */*.cc) 
foreach(SOURCE ${EXAMPLE_SOURCES})
string(REGEX REPLACE "^.*/([^/]*)\\.cc$" "\\1" PROGNAME "${SOURCE}")
if (NOT PROGNAME STREQUAL ellipsoid) #ellipsoid is broken
	string(REGEX REPLACE "^.*/(.*)/${PROGNAME}\\.cc$" "\\1" DIRNAME "${SOURCE}")
	add_executable(${PROGNAME} ${SOURCE})
	target_link_libraries(${PROGNAME} PRIVATE voro++)
	set_target_properties(${PROGNAME} PROPERTIES 
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DIRNAME}" )
endif()
endforeach(SOURCE)
