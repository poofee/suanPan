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

#include "Frequence.h"
#include <Solver/Solver.h>

Frequence::Frequence(const unsigned T, const unsigned N)
    : Step(T, CT_FREQUENCE, 0.)
    , eigen_number(N) {}

int Frequence::initialize() { return 0; }

int Frequence::analyze() { return get_solver()->analyze(); }

void Frequence::set_eigen_number(const unsigned N) { eigen_number = N; }

unsigned Frequence::get_eigen_number() const { return eigen_number; }
