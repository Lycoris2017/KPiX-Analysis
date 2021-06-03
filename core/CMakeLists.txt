# Required packages
FIND_PACKAGE(LibXml2 REQUIRED)
FIND_PACKAGE(ROOT REQUIRED COMPONENTS vectorDict Minuit Minuit2 Gui GenVector Geom Core Hist RIO NO_MODULE)

INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS})
message("--- ${LIBXML2_INCLUDE_DIR}")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)


# glob the files for the library
file(GLOB_RECURSE kpixData_SOURCES *.cpp)
file(GLOB_RECURSE kpixData_HEADERS *.hpp *.h)
# remove euDataRead and xml parts
list(REMOVE_ITEM kpixData_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/euDataRead.cpp)
list(REMOVE_ITEM kpixData_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/euDataRead.h)
#list(REMOVE_ITEM kpixData_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/XmlVariables.cpp)
#list(REMOVE_ITEM kpixData_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/XmlVariables.h)

# create the library
add_library(kpixData SHARED ${kpixData_SOURCES})

#target_include_directories(kpixData PUBLIC ${LIBXML2_INCLUDE_DIR})
# link external packages
target_link_libraries(kpixData ROOT::Core ROOT::vectorDict ROOT::Tree ROOT::GenVector ROOT::Geom ROOT::RIO ROOT::MathCore ROOT::Hist ROOT::GuiBld ROOT::Minuit ROOT::Minuit2 ${LIBXML2_LIBRARY})

