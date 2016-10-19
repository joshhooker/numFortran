#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(NUMFORTRAN)
enable_language(Fortran)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS"-O3 -std=f2008")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS "-O3 -recursive")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.f90)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(libraryTest libraryTest.f90 ${sources})
#target_link_libraries(libraryTest ${Geant4_LIBRARIES} ${GSL_LIBRARIES} ${ROOT_LIBRARIES})

#
#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(NUMFORTRAN DEPENDS libraryTest)

# # CMake project file for libraryTest

# cmake_minimum_required (VERSION 2.6)
# project (libraryTest)
# enable_language (Fortran)

# # make sure that the default is a RELEASE
# if (NOT CMAKE_BUILD_TYPE)
#   set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
#       "Choose the type of build, options are: None Debug Release."
#       FORCE)
# endif (NOT CMAKE_BUILD_TYPE)

# # default installation
# get_filename_component (default_prefix ".." ABSOLUTE)
# set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
#       "Choose the installation directory; by default it installs in the NORMA directory."
#       FORCE)

# # FFLAGS depend on the compiler
# get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#   # gfortran
#   set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
#   set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
# elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
#   # ifort (untested)
#   set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
#   set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
# elseif (Fortran_COMPILER_NAME MATCHES "g77")
#   # g77
#   set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
#   set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
# else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#   message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
#   message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
#   message ("No optimized Fortran compiler flags are known, we just try -O2...")
#   set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
#   set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
# endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")


# # build executables
# set (PROGRAMS "nfConstants" "nfIntegration")
# set (EXECUTABLES "libraryTest" ${PROGRAMS})

# add_executable ("libraryTest" "libraryTest.f90")
# foreach (p ${PROGRAMS})
#   add_executable (${p} "${p}.f90")
# endforeach (p)

# # install executables and scripts
# install (TARGETS ${EXECUTABLES}
#          RUNTIME DESTINATION "bin")
# install (PROGRAMS ${SCRIPTS}
#          DESTINATION "bin")