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

#include "ODE_Solver.h"
#include <Solver/ODE_Solver/ODE.h>

ODE_Solver::ODE_Solver(const unsigned T, const unsigned CT, ODE* O)
    : Tag(T, CT)
    , ode_system(O) {
    suanpan_debug("ODE_Solver $u ctor called.\n", T);
}

ODE_Solver::~ODE_Solver() { suanpan_debug("ODE_Solver $u dtor called.\n", get_tag()); }

int ODE_Solver::initialize() const {
    if(ode_system == nullptr) {
        suanpan_error("initialize() needs a valid ODE.\n");
        return -1;
    }

    return 0;
}

int ODE_Solver::update_status() { throw; }

int ODE_Solver::analyze() { throw; }

void ODE_Solver::set_ode(ODE* E) { ode_system = E; }

ODE* ODE_Solver::get_ode() const { return ode_system; }

unique_ptr<ODE_Solver> suanpan::make_copy(const unique_ptr<ODE_Solver>& old_obj) { return old_obj->get_copy(); }
