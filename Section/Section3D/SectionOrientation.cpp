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

#include "SectionOrientation.h"
#include <utility>

SectionOrientation::SectionOrientation(const unsigned T, const double X, const double Y, const double Z)
    : Tag(T, CT_SECTIONORIENTATION)
    , orientation(std::initializer_list<double>{ X, Y, Z }) {}

SectionOrientation::SectionOrientation(const unsigned T, vec O)
    : Tag(T, CT_SECTIONORIENTATION)
    , orientation(std::move(O)) {}

const vec& SectionOrientation::get_orientation() const { return orientation; }
