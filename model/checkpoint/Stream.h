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

#pragma once

// The one stream helper of the checkpoint module (model/checkpoint/): one
// Archive class covering both directions, typed primitives only — no
// operator& zoo, no templates over the stream type at call sites.
//
// Wire primitives: fixed-size values (memcpy, native endianness — files are
// read back by the same build), bool as u8, enum as its underlying type,
// string as u32 length + bytes, vector as u64 count + elements, RNG as its
// raw 32-byte state (via RandomBase::checkpoint). Reads are budget-checked:
// a truncated or corrupt file throws checkpoint_error, never UB.

#include "util/errors.h"

#include <cstddef>
#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

namespace OpenMalaria {

/// Bad/incompatible checkpoint file (clean model error, exit code Checkpoint).
class checkpoint_error : public cmd_exception {
  public:
    explicit checkpoint_error(const std::string& msg) : cmd_exception(msg, Error::Checkpoint) {}
};

class Archive {
  public:
    /// Read mode: budget = whole file (stream must be at its start).
    explicit Archive(std::istream& is) : is_(&is) {
        is_->clear();
        is_->seekg(0, std::ios::end);
        end_ = is_->tellg();
        is_->clear();
        is_->seekg(0, std::ios::beg);
        if (is_->fail())
            throw checkpoint_error("cannot size checkpoint input");
    }
    /// Write mode.
    explicit Archive(std::ostream& os) : is_(nullptr), os_(&os) {}

    bool reading() const { return is_ != nullptr; }

    /// Raw block of n bytes (both directions; the escape hatch for POD
    /// payloads — caller owns the byte layout, keep it to simple structs).
    void bytes(void* p, std::size_t n) {
        if (reading()) {
            if (n > remaining())
                throw checkpoint_error("out of data (wanted " + std::to_string(n) + " bytes, " +
                                       std::to_string(remaining()) + " left at pos " + std::to_string(is_->tellg()) +
                                       ")");
            is_->read(static_cast<char*>(p), std::streamsize(n));
            if (is_->gcount() != std::streamsize(n))
                throw checkpoint_error("stream read error");
        } else {
            os_->write(static_cast<const char*>(p), std::streamsize(n));
            if (os_->fail())
                throw checkpoint_error("stream write error");
        }
    }

    /// One fixed-size value (integers, doubles, SimTime, pointers off-wire).
    template <class T> void value(T& v) {
        static_assert(std::is_trivially_copyable<T>::value, "value: POD only");
        static_assert(!std::is_same<T, bool>::value, "bool: use boolean()");
        bytes(&v, sizeof(T));
    }

    /// bool on the wire is exactly one byte (0/1).
    void boolean(bool& v) {
        uint8_t b = v ? 1 : 0;
        value(b);
        v = b != 0;
        if (reading() && b > 1)
            throw checkpoint_error("bool out of range");
    }

    /// Enum as its underlying integer (range-checked on read by the caller
    /// if the domain matters).
    template <class E> void enumerated(E& e) {
        static_assert(std::is_enum<E>::value, "enumerated: enum only");
        using U = typename std::underlying_type<E>::type;
        U u = static_cast<U>(e);
        value(u);
        e = static_cast<E>(u);
    }

    /// string: u32 length + bytes.
    void text(std::string& s) {
        uint32_t n = uint32_t(s.size());
        value(n);
        if (reading()) {
            if (n > remaining())
                throw checkpoint_error("string longer than remaining data");
            s.resize(n);
            if (n)
                bytes(&s[0], n);
        } else {
            if (n)
                bytes(&s[0], n);
        }
    }

    /// vector of fixed-size values: u64 count + elements; on read the count
    /// is rejected unless that many elements could still fit the file.
    template <class T> void vector(std::vector<T>& v) {
        static_assert(std::is_trivially_copyable<T>::value, "vector: POD only");
        uint64_t n = v.size();
        value(n);
        if (reading()) {
            if (n > remaining() / sizeof(T))
                throw checkpoint_error("vector<" + std::to_string(sizeof(T)) + "> claims " + std::to_string(n) +
                                       " elements; file too small");
            v.resize(size_t(n));
            if (n)
                bytes(v.data(), size_t(n * sizeof(T)));
        } else {
            if (n)
                bytes(v.data(), size_t(n * sizeof(T)));
        }
    }

    /// RNG state: the raw 32-byte xoshiro state (binary_checkpoint).
    /// Works for MasterRng, LocalRng — anything with checkpoint(stream).
    template <class R> void rng(R& r) {
        if (reading())
            r.checkpoint(*is_);
        else
            r.checkpoint(*os_);
    }

    /// Bytes left to read (infinite for write mode).
    std::streamoff remaining() const { return reading() ? end_ - is_->tellg() : std::streamoff(-1); }

    /// Fail if anything was left unread (caller bug or file mismatch).
    void checkExhausted() {
        if (!reading())
            return;
        if (remaining() != 0)
            throw checkpoint_error("checkpoint file has " + std::to_string(remaining()) + " bytes remaining");
    }

  private:
    std::istream* is_;
    std::ostream* os_ = nullptr;
    std::streamoff end_ = 0;
};

} // namespace OpenMalaria
