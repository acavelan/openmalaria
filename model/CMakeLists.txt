# CMake configuration for all code in this dir and below

# This file is part of OpenMalaria.
# 
# Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
# Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
# Copyright (C) 2020-2026 University of Basel
# Copyright (C) 2025-2026 The Kids Research Institute Australia
#
# OpenMalaria is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

# Configure this file
file(STRINGS ${CMAKE_SOURCE_DIR}/version.txt OM_VERSION)
# make sure version.h gets rebuilt everytime cmake is run
unset(${CMAKE_CURRENT_BINARY_DIR}/util/version.h CACHE)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/util/version.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/util/version.h
    @ONLY
)

# Plan §7 row 2: GLOB replaces the hand-maintained list (phases 2-7 add/remove
# files freely). main.cpp is the openMalaria executable's entry point (root
# CMakeLists), never part of the lib. StreamValidator is opt-in (below).
file (GLOB_RECURSE Model_CPP CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
list (FILTER Model_CPP EXCLUDE REGEX "/(main|StreamValidator|DescriptiveIPTInfection)\.cpp$")
# DescriptiveIPTInfection.cpp: dead since the IPTDescription schema type was
# removed (census N1); excluded until deleted in this phase, step 2.

if (${OM_STREAM_VALIDATOR})
  list(APPEND Model_CPP util/StreamValidator.cpp)
endif (${OM_STREAM_VALIDATOR})
# Headers - only included so they show up in IDEs:
# This misses loads of headers. Fix if you care.
file (GLOB_RECURSE Model_H "${CMAKE_SOURCE_DIR}/model/*.h")


include_directories (SYSTEM
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${Z_INCLUDE_DIRS}
  # ${SPIRIT_INCLUDE_DIRS}
  
  ${CMAKE_BINARY_DIR}/model
  ${CMAKE_SOURCE_DIR}/contrib
  # ${CMAKE_SOURCE_DIR}/contrib/floating_point_utilities_v3
)
include_directories (
  ${CMAKE_SOURCE_DIR}/model ${CMAKE_BINARY_DIR}
)

add_library (model
  ${Model_CPP}
  ${Model_H}
  ${CMAKE_CURRENT_BINARY_DIR}/util/version.h
)

if (MSVC)
  set_target_properties (model PROPERTIES
    LINK_FLAGS "${OM_LINK_FLAGS}"
    COMPILE_FLAGS "${OM_COMPILE_FLAGS}"
  )
endif (MSVC)
