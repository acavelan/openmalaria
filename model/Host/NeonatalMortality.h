/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2025-2026 The Kids Research Institute Australia
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef Hmod_NeonatalMortality
#define Hmod_NeonatalMortality

#include "Clock.h"
#include "util/random.h"
#include "Host/Human.h"

namespace scnXml {
    class Clinical;
}
namespace OM {
namespace checkpoint { class Archive; }
    class Population;
namespace Host {

using util::LocalRng;

class NeonatalMortality {
public:
  /// Initialisation
  static void init( RunContext& ctx, const scnXml::Clinical& clinical );
  
  /** Called for each birth; returns true if infant dies due to mother's
   * infection. */
  static bool eventNeonatalMortality(LocalRng& rng);
  
  /** Calculate risk of a neonatal mortality based on humans 20-25 years old. */
  static void update (RunContext& ctx, vector<Human>& population);
  
  /// Run state (risk + prevalence ring) to/from checkpoint wire
  /// (model/checkpoint/statics.cpp).
  static void staticCheckpoint (checkpoint::Archive& a);
};

} }
#endif
