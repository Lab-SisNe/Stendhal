#[[
filename: CMakeLists.txt

This code simulates the PD microcircuit model (Potjans and Diesmann, 2014)
with the GL (Galves and Locherbach, 2013) neuron in discrete time.
The leak function is an exponentially decaying function.
For the pseudo random number generation, the xoroshiro128+ algorithm was used.

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
project( stendhal VERSION 0.2 )

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
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra ${with-optimize}" )
else()
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra" )
endif()

# Set "-fp-model strict" if using Intel compiler
if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model strict" )
endif()

# Add executable
add_executable( dGLPD_delta main_dGLPD_delta.cpp dGLPD_delta.cpp gl_psc_delta.cpp)

# Add directory of config header to the list of paths to search for include files
target_include_directories( dGLPD_delta PUBLIC "${PROJECT_BINARY_DIR}" )

