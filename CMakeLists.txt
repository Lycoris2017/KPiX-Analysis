
#------------------
#- Author: Mengqing Wu <mengqing.wu@desy.de>
#- Generated: September, 2018
#- Modified: 2019-08-02
# * new feature:
#   kpix package is importable from outside via find_package()
#- for LYCORIS project at DESY.
#------------------


##-- Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
include(CMakeDependentOption)

##-- Project's name

project(Lycoris VERSION 1.0)

##-- Options for kpix core packages
option(USE_KPIX 
  "Use local KPiX core libraries to decode data!" ON)

option(MODE_DEBUG
  "Use DEBUG mode, ONLY single module to compile!" OFF)

option(USE_YAML "Turn libyaml-dependency " OFF)

#option(USE_ROOT_INSTALL "Install shared lib and headers to /usr/local/include(or lib)/kpix/" OFF)

option(DO_EXPORT "Export kpix pacakge " ON)

##-- cmake setup
set(CMAKE_BUILD_TYPE Debug )

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}" CACHE PATH "Install prefix" FORCE )

set(LIBRARY_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH ${LIBRARY_INSTALL_RPATH})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

##-- C11 setup
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ADD_DEFINITIONS("-Wunused-value -Wall -pedantic -Wshadow -fdiagnostics-color=auto") 


##-- Flags
find_package(LibXml2 REQUIRED)
find_package(ROOT REQUIRED)
#find_package(ROOT COMPONENTS RooStats)
#find_package(eudaq)
if (eudaq_FOUND)
  message("debug: ${EUDAQ_LIBRARIES_DIRS} ${EUDAQ_LIBRARIES}")
endif(eudaq_FOUND)

##-- everything needed for YAML: 'libyaml-dev' used
if (USE_YAML)
  # w/o Anaconda, use this other than the 2 below lines!
  #FIND_PATH(LIBYAML_INCLUDE_DIR NAMES yaml.h) 
  #FIND_LIBRARY(LIBYAML_LIBRARIES NAMES yaml libyaml)
  set(LIBYAML_INCLUDE_DIR /usr/include)
  set(LIBYAML_LIBRARIES /usr/lib/x86_64-linux-gnu/libyaml.so)
  
  INCLUDE(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(Yaml DEFAULT_MSG LIBYAML_LIBRARIES LIBYAML_INCLUDE_DIR)
  MARK_AS_ADVANCED(LIBYAML_INCLUDE_DIR LIBYAML_LIBRARIES)
  #message("[debug] LIBYAML_INCLUDE_DIR = ${LIBYAML_INCLUDE_DIR}")
  #message("[debug] LIBYAML_LIBRARIES = ${LIBYAML_LIBRARIES}")

endif(USE_YAML)

##-- Library folder to be included
include_directories( ${PROJECT_SOURCE_DIR}/include )
include_directories( ${ROOT_INCLUDE_DIR} )
include_directories( ${LIBXML2_INCLUDE_DIR} )


##-- KPiX Core Packages to install prior
if (USE_KPIX)
  message(" [INFO] local KPiX lite package used.")

  set(KPIX_LOCAL_LIB kpix)
  include_directories( ${PROJECT_SOURCE_DIR}/core/include)

  aux_source_directory( ${PROJECT_SOURCE_DIR}/core/src LIB_SRC)
  
  if (NOT eudaq_FOUND)
    message("debug: remove the euDataRead.cxx")
    list(REMOVE_ITEM LIB_SRC ${PROJECT_SOURCE_DIR}/core/src/euDataRead.cpp)
  endif()
 # message(${LIB_SRC})
  
  add_library(${KPIX_LOCAL_LIB} SHARED ${LIB_SRC})
  set_target_properties(${KPIX_LOCAL_LIB} PROPERTIES SOVERSION 2)

  target_include_directories(${KPIX_LOCAL_LIB}
	PUBLIC
	$<INSTALL_INTERFACE:core/include>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/include>
	PRIVATE
	)

  install(TARGETS ${KPIX_LOCAL_LIB}
	EXPORT kpixConfig
	LIBRARY DESTINATION lib
	)
  # file(GLOB INC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/eudaq/*.hh")

  if (DO_EXPORT)
	
	#--newly added for find_package() outside of kpix
	export(TARGETS
	  ${KPIX_LOCAL_LIB}
      NAMESPACE kpix::
      FILE "${PROJECT_SOURCE_DIR}/cmake/kpixConfig.cmake"
	  )
	install(EXPORT
      kpixConfig
      DESTINATION "${PROJECT_SOURCE_DIR}/cmake"
      NAMESPACE kpix::
	  )
	
  else()
    # start of adding files to usr/local
    file(GLOB KPIX_HEADER "${PROJECT_SOURCE_DIR}/core/include/*.h" )
    install(FILES ${KPIX_HEADER} DESTINATION /usr/local/include/kpix)
    install(TARGETS ${KPIX_LOCAL_LIB} LIBRARY DESTINATION /usr/local/lib/kpix)
    # end of adding files to usr/local
	
  endif(DO_EXPORT)
    
  
else () # not tested, definately not working...
  message(" [INFO] external KPiX packages used.")

  set(CMAKE_PREFIX_PATH "/usr/local/lib/kpix")
  set(KPIX_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/lib")
  set(CMAKE_INSTALL_RPATH ${KPIX_INSTALL_RPATH})
  include_directories( /usr/local/include/kpix )
  find_library(KPIX_EXT_LIB kpix)
  
endif (USE_KPIX)

##-- libraries to link against:
if (DEFINED KPIX_LOCAL_LIB)
  list(APPEND LYCORIS_LIBRARIES ${KPIX_LOCAL_LIB})
endif()

if (DEFINED KPIX_EXT_LIB)
  list(APPEND LYCORIS_LIBRARIES ${KPIX_EXT_LIB})
endif()

##-- LinkDef.h
include(${ROOT_USE_FILE})
#message("\tDebug: ${ROOT_USE_FILE}")
ROOT_GENERATE_DICTIONARY(G__Dict ${PROJECT_SOURCE_DIR}/include/ntupleMaker.h LINKDEF  ${PROJECT_SOURCE_DIR}/include/LinkDef.h)

add_library(Dict SHARED G__Dict.cxx)
target_link_libraries(Dict ${ROOT_LIBRARIES})

install( TARGETS Dict LIBRARY DESTINATION lib )
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/libDict_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/libDict.rootmap DESTINATION lib )

##-- add exe: # maybe cmake can be scripted with a loop using aux_source_directory func?

#add_executable(ymlHandler ${PROJECT_SOURCE_DIR}/src/ymlHandler.cxx)
#target_link_libraries(ymlHandler ${ROOT_LIBRARIES})
#list(APPEND INSTALL_EXE_TARGETS ymlHandler)

# simple evt loop code:

add_executable(printDat ${PROJECT_SOURCE_DIR}/src/printDat.cxx)
target_link_libraries(printDat ${ROOT_LIBRARIES}  ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES})
list(APPEND INSTALL_EXE_TARGETS printDat)

add_executable(count ${PROJECT_SOURCE_DIR}/src/count.cxx)
target_link_libraries(count ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS count)

add_executable(new_analysis ${PROJECT_SOURCE_DIR}/src/new_analysis.cxx)
target_link_libraries(new_analysis ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS new_analysis)

add_executable(new_cluster_analysis ${PROJECT_SOURCE_DIR}/src/new_cluster_analysis.cxx)
target_link_libraries(new_cluster_analysis ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS new_cluster_analysis)

# add_executable(checkTimeVars ${PROJECT_SOURCE_DIR}/src/checkTimeVars.cxx)
# target_link_libraries(checkTimeVars ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
# list(APPEND INSTALL_EXE_TARGETS checkTimeVars)

# calibration ana:

add_executable(new_CalibrationFitter ${PROJECT_SOURCE_DIR}/src/new_CalibrationFitter.cxx)
target_link_libraries(new_CalibrationFitter ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS new_CalibrationFitter)

add_executable(ymlCalibrationFitter ${PROJECT_SOURCE_DIR}/src/ymlCalibrationFitter.cxx)
target_link_libraries(ymlCalibrationFitter ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS ymlCalibrationFitter)

add_executable(DuoKpix_CalibFitter ${PROJECT_SOURCE_DIR}/src/DuoKpix_CalibFitter.cxx)
target_link_libraries(DuoKpix_CalibFitter ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS DuoKpix_CalibFitter)

add_executable(calibrationFitter ${PROJECT_SOURCE_DIR}/src/calibrationFitter.cxx)
target_link_libraries(calibrationFitter ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS calibrationFitter)


# full ana:

add_executable(analysisLite ${PROJECT_SOURCE_DIR}/src/analysisLite.cxx)
target_link_libraries(analysisLite ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS analysisLite)

add_executable(analysis ${PROJECT_SOURCE_DIR}/src/analysis.cxx)
target_link_libraries(analysis ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS analysis)

# add_executable(analysis_newdaq ${PROJECT_SOURCE_DIR}/src/analysis_newdaq.cxx)
# target_link_libraries(analysis_newdaq ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
# list(APPEND INSTALL_EXE_TARGETS analysis_newdaq)

#add_executable(splitFile ${PROJECT_SOURCE_DIR}/src/splitFile.cxx)
#target_link_libraries(splitFile ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
#list(APPEND INSTALL_EXE_TARGETS splitFile)

# add_executable(analysis_newdaq_mq ${PROJECT_SOURCE_DIR}/src/analysis_newdaq_mq.cxx)
# target_link_libraries(analysis_newdaq_mq ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
# list(APPEND INSTALL_EXE_TARGETS analysis_newdaq_mq)

add_executable(cluster_analysis ${PROJECT_SOURCE_DIR}/src/cluster_analysis.cxx)
target_link_libraries(cluster_analysis ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS cluster_analysis)

add_executable(cluster_analysis_step2 ${PROJECT_SOURCE_DIR}/src/cluster_analysis_step2.cxx)
target_link_libraries(cluster_analysis_step2 ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS cluster_analysis_step2)

add_executable(cluster_analysis_Lite ${PROJECT_SOURCE_DIR}/src/cluster_analysis_Lite.cxx)
target_link_libraries(cluster_analysis_Lite ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS cluster_analysis_Lite)

add_executable(external_trig_debug ${PROJECT_SOURCE_DIR}/src/external_trig_debug.cxx)
target_link_libraries(external_trig_debug ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS external_trig_debug)


#add_executable(analysisExternal ${PROJECT_SOURCE_DIR}/src/analysisExternal.cxx)
#target_link_libraries(analysisExternal ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
#list(APPEND INSTALL_EXE_TARGETS analysisExternal)

add_executable(analysisExternal_tree ${PROJECT_SOURCE_DIR}/src/analysisExternal_tree.cxx)
target_link_libraries(analysisExternal_tree ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS analysisExternal_tree)

#add_executable(pedestal_old ${PROJECT_SOURCE_DIR}/src/pedestal_old.cxx)
#target_link_libraries(pedestal_old ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
#list(APPEND INSTALL_EXE_TARGETS pedestal_old)

add_executable(pedestal_tree ${PROJECT_SOURCE_DIR}/src/pedestal_tree.cxx)
target_link_libraries(pedestal_tree ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
list(APPEND INSTALL_EXE_TARGETS pedestal_tree)



# tree builders:

#add_executable(ntupleTest ${PROJECT_SOURCE_DIR}/src/ntupleTest.cxx)
#target_link_libraries(ntupleTest ${ROOT_LIBRARIES} ${LYCORIS_LIBRARIES} ${LIBXML2_LIBRARIES} Dict)
#list(APPEND INSTALL_EXE_TARGETS ntupleTest)




# install 
install(
 TARGETS ${INSTALL_EXE_TARGETS}
 RUNTIME DESTINATION bin
 )
