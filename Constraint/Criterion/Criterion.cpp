////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2018 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "Criterion.h"

Criterion::Criterion(const unsigned& T, const unsigned& CT, const unsigned& ST)
    : Tag(T, CT)
    , step_tag(ST) {
    suanpan_debug("Criterion %u ctor() called.\n", T);
}

Criterion::~Criterion() { suanpan_debug("Criterion %u dtor() called.\n", get_tag()); }

void Criterion::set_step_tag(const unsigned& T) { step_tag = T; }

const unsigned& Criterion::get_step_tag() const { return step_tag; }
