cmake_minimum_required(VERSION 3.10)

project(clistones C)

find_package(PkgConfig REQUIRED) 

pkg_check_modules(SIGUTILS REQUIRED sigutils)
pkg_check_modules(ALSA REQUIRED alsa)

set(SRCDIR src)
set(INCLUDEDIR include)

set(CLISTONES_HEADERS
  ${INCLUDEDIR}/graves.h)
  
set(CLISTONES_SOURCES
  ${SRCDIR}/graves.c
  ${SRCDIR}/main.c)
  
add_executable(
  clistones
  ${CLISTONES_HEADERS}
  ${CLISTONES_SOURCES})

target_link_libraries(
  clistones 
  ${SIGUTILS_LIBRARIES} 
  ${ALSA_LIBRARIES})
  
target_include_directories(
  clistones PUBLIC 
  ${SIGUTILS_INCLUDE_DIRS}
  ${ALSA_INCLUDE_DIRS}
  ${INCLUDEDIR})
        
target_compile_options(
  clistones PUBLIC
  ${SIGUTILS_CFLAGS_OTHER}
  ${ALSA_CFLAGS_OTHER})
