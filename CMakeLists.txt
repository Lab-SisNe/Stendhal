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
project( Stendhal VERSION 0.1 )

# Specify the C++ standard
set( CMAKE_CXX_STANDARD 11)
set( CMAKE_CXX_STANDARD_REQUIRED True )

# Add header file to pass configuration
configure_file( stendhalconfig.h.in stendhalconfig.h )

# Add directory of config header to the list of paths to search for include files
#target_include_directories( stendhal PUBLIC "${PROJECT_BINARY_DIR}" )

# Add library
#add_library( PD_GL PD_GL.cpp )
add_library( gl_psc_exp gl_psc_exp.cpp )

# Add executable
#add_executable( stendhal_PD_GL main.cpp PD_GL.cpp )
#add_executable( stendhal_PD_GL main.cpp )
add_executable( a.out test_gl_psc_exp.cpp )
