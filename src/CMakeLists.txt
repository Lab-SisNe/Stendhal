#[[
filename: CMakeLists.txt

This code simulates the PD microcircuit model (Potjans and Diesmann, 2014)
with the GL (Galves and Locherbach, 2013) neuron in discrete time.

This code is part of STENDHAL package.
A NeuroMat package for simulate neuron networks, and analyse data.

Director:
Antonio Galves

Developers:
Nilton L. Kamiji
Christophe Pouzat
Renan Shimoura

Contributors:
Karine Guimaraes
Aline Duarte
Jorge Stolfi
Antonio Roque

Aug 11, 2020
]]

cmake_minimum_required( VERSION 3.15 )

# Set project name and version
project( stendhal VERSION 0.1 )

# Specify the C++ standard
set( CMAKE_CXX_STANDARD 11)
set( CMAKE_CXX_STANDARD_REQUIRED True )

# Add header file to pass configuration
configure_file( stendhalconfig.h.in stendhalconfig.h )

option( with-optimize "Enable user defined optimization. [defaul ON, when ONdefaults to '-O3']" ON )
# Set optimzation flag. Defaults to -O3
if ( with-optimize )
  if ( ${with-optimize} STREQUAL ON )
    set( with-optimize -O3 )
  endif()
  set( CMAKE_CXX_FLAGS "-Wall -Wextra ${with-optimize}" )
else()
  set( CMAKE_CXX_FLAGS "-Wall -Wextra" )
endif()

# Set "-fp-model strict" if using Intel compiler
if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model strict" )
endif()

# OPTIONS:
# with_pcg option
option(with-pcg "Use PCG as random number generator engine" FALSE)
if ( ${with-pcg} )
  # Check if PCG is present and set compiler flag accordinglly
  find_file(USE_PCG pcg-cpp/pcg_random.hpp "${PROJECT_SOURCE_DIR}")
  if (NOT ${USE_PCG} MATCHES "USE_PCG-NOTFOUND")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_PCG")
  endif()
endif()

# Add executable
add_executable( stendhal_PD main_dPD_GL.cpp dPD_GL.cpp gl_psc_exp.cpp)

# Add directory of config header to the list of paths to search for include files
target_include_directories( stendhal_PD PUBLIC "${PROJECT_BINARY_DIR}" )

