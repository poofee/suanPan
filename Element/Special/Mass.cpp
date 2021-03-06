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

#include "Mass.h"

Mass::Mass(const unsigned& T, const unsigned& NT, const double& MA, const uvec& DT)
    : Element(T, ET_MASS, 1, unsigned(DT.max()), uvec{ NT })
    , magnitude(MA)
    , dof_label(DT - 1) {}

void Mass::initialize(const shared_ptr<DomainBase>&) {
    trial_mass.zeros(dof_label.max() + 1, dof_label.max() + 1);
    for(const auto& I : dof_label) trial_mass(I, I) = magnitude;
}

int Mass::update_status() { return 0; }

int Mass::commit_status() { return 0; }

int Mass::clear_status() { return 0; }

int Mass::reset_status() { return 0; }
