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

#ifndef OM_util_checkpoint_containers
#define OM_util_checkpoint_containers

#include "util/checkpoint.h"
#include <list>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

/** Provides some extra functions. See checkpoint.h. */
namespace OM {
namespace util {
namespace checkpoint {
    ///@brief Operator& for smart pointers
    //@{
    template<class T>
    void operator& (std::unique_ptr<T> & x, std::ostream& stream) {
        (*x) & stream;
    }

    template<class T>
    void operator& (std::unique_ptr<T>& x, std::istream& stream) {
        (*x) & stream;
    }
    //@}
    
    ///@brief Operator& for stl containers
    //@{
    template<class U, class V, class S>
    void operator& (std::pair<U,V> x, S& stream) {
        x.first & stream;
        x.second & stream;
    }
    
    template<class T>
    void operator& (std::vector<T>& x, std::ostream& stream) {
        x.size() & stream;
        for (T& y : x) {
            y & stream;
        }
    }

    template<class T>
    void operator& (std::vector<T>& x, std::istream& stream) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.resize (l);
        for (T& y : x) {
            y & stream;
        }
    }

    /// Version of above taking an element to initialize each element from.
    template<class T>
    void checkpoint (std::vector<T>& x, std::istream& stream, T templateInstance) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.resize (l, templateInstance);
        for (T& y : x) {
            y & stream;
        }
    }
    
    template<class T>
    void operator& (std::list<T> x, std::ostream& stream) {
        x.size() & stream;
        for (T& y : x) {
            y & stream;
        }
    }
    template<class T>
    void operator& (std::list<T>& x, std::istream& stream) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.resize (l);
        for (T& y : x) {
            y & stream;
        }
    }

    template<class K, class V>
    void operator& (std::map<K,V> const& m, std::ostream& stream) {
        size_t n = m.size();
        n & stream; // write count
        for (auto const& kv : m) {
            kv.first  & stream;
            kv.second & stream;
        }
    }

    template<class K, class V>
    void operator& (std::map<K,V>& m, std::istream& stream) {
        size_t n;
        n & stream; // read count
        validateListSize(n);
        m.clear();
        for (size_t i = 0; i < n; ++i) {
            K key{stream};
            V val;
            val & stream;
            m.emplace(std::move(key), std::move(val));
        }
    }
    
    void operator& (const std::set<interventions::ComponentId>&x, std::ostream& stream);
    void operator& (std::set<interventions::ComponentId>& x, std::istream& stream);
    
    void operator& (const std::map<std::string,double>& x, std::ostream& stream);
    void operator& (std::map<std::string, double >& x, std::istream& stream);
    
    void operator& (const std::map<double,double>& x, std::ostream& stream);
    void operator& (std::map<double, double>& x, std::istream& stream);
    
    void operator& (const std::map<interventions::ComponentId,SimTime>& x, std::ostream& stream);
    void operator& (std::map<interventions::ComponentId,SimTime>& x, std::istream& stream);
    
    void operator& (const std::multimap<double,double>& x, std::ostream& stream);
    void operator& (std::multimap<double, double>& x, std::istream& stream);
    //@}
    
} } }   // end of namespaces
#endif
