cmake_minimum_required(VERSION 3.13)

project(gausskruger CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_CLI "Build CLI application")

add_compile_options(-Wall -Wextra -pedantic -Werror -Wnarrowing)

add_library(gausskruger INTERFACE)

set_target_properties(gausskruger PROPERTIES PUBLIC_HEADER gausskruger.h)
target_include_directories(gausskruger
                           INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                            $<INSTALL_INTERFACE:include>
                           )

if(${BUILD_CLI})
  find_package(Boost REQUIRED COMPONENTS program_options)

  add_executable(gausskruger_cli
      gausskruger_cli.cpp
  )

  target_link_libraries(gausskruger_cli
    PRIVATE
      gausskruger
      Boost::program_options
  )
endif()

install(FILES gausskruger.h DESTINATION include)

install(TARGETS gausskruger
        EXPORT gausskruger)

install(EXPORT gausskruger
  NAMESPACE gausskruger::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake
)
