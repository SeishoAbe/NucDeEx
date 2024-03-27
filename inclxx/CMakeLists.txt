project(NUCDEEX CXX)
cmake_minimum_required(VERSION 2.6)

if(INCL_DEEXCITATION_NUCDEEX)

  # Make the include directories available to the top level
  set_property(DIRECTORY ../.. APPEND
    PROPERTY INCLUDE_DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/interface/include
    ${CMAKE_CURRENT_SOURCE_DIR}/upstream/include
  )
  include_directories(interface/include upstream/include)

  # NUCDEEX library
  file (GLOB NUCDEEX_SRC upstream/src/*.cc interface/src/*.cc)
  list(REMOVE_ITEM NUCDEEX_SRC ${CMAKE_CURRENT_SOURCE_DIR}/upstream/src/G4NucDeExInterface.cc)
  list(REMOVE_ITEM NUCDEEX_SRC ${CMAKE_CURRENT_SOURCE_DIR}/upstream/src/G4NucDeExInterfaceMessenger.cc)
  list(REMOVE_ITEM NUCDEEX_SRC ${CMAKE_CURRENT_SOURCE_DIR}/upstream/src/NEUTNucDeExInterace.cc)
  list(REMOVE_ITEM NUCDEEX_SRC ${CMAKE_CURRENT_SOURCE_DIR}/upstream/src/ReadTALYS.cc)
  add_library(NUCDEEX ${NUCDEEX_SRC})
  target_link_libraries(NUCDEEX DeExcitation)

endif(INCL_DEEXCITATION_NUCDEEX)
