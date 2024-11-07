cmake_minimum_required(VERSION 3.20)

project(tsmp LANGUAGES CXX)


IF(WIN32)
    message ( "WINDOWS build" )
    # on WINDOWS it is easier to make STATIC libs
    # some libs will by dynamic (dll); add PATH to the lib folder
    add_compile_definitions( _USE_MATH_DEFINES _msvc )
    IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
        message ("debug build ---------------------")
        add_compile_options(-DDEBUG)
        set(CMAKE_VERBOSE_MAKEFILE ON)
    ENDIF(CMAKE_BUILD_TYPE MATCHES DEBUG)
    IF(CMAKE_BUILD_TYPE MATCHES RELEASE)
        message ("release build ---------------------")
        add_compile_options(-DNDEBUG)
    ENDIF(CMAKE_BUILD_TYPE MATCHES RELEASE)
ENDIF(WIN32)


set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt6
    REQUIRED COMPONENTS Core
)

#top tree dir is CMAKE_SOURCE_DIR
include_directories(
${CMAKE_SOURCE_DIR}/fftReal
)

file(GLOB SOURCES ${CMAKE_SOURCE_DIR}/*.cpp)

file(GLOB HEADERS ${CMAKE_SOURCE_DIR}/*.h)

set(PROJECT_SOURCES ${HEADERS} ${SOURCES} )
add_executable(${PROJECT_NAME} ${PROJECT_SOURCES} )

target_link_libraries(${PROJECT_NAME} PRIVATE
    Qt6::Core
)

install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
