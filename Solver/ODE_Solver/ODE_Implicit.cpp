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

#include "ODE_Implicit.h"

ODE_Implicit::ODE_Implicit(const unsigned T, const unsigned CT, const unsigned N, const bool C, ODE* O)
    : ODE_Solver(T, CT, O)
    , step_number(N)
    , use_corrector(C) {}

ODE_Implicit::~ODE_Implicit() {}

int ODE_Implicit::analyze() { return update_status(); }
