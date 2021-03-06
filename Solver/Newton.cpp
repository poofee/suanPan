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

#include "Newton.h"
#include <Converger/Converger.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>

Newton::Newton(const unsigned T, const bool IS)
    : Solver(T, CT_NEWTON)
    , initial_stiffness(IS) {}

int Newton::analyze() {
    auto& C = get_converger();
    auto& G = get_integrator();
    const auto& W = G->get_domain().lock()->get_factory();

    suanpan_info("current analysis time: %.5f.\n", W->get_trial_time());

    auto& max_iteration = C->get_max_iteration();

    // iteration counter
    unsigned counter = 0, flag;

    while(true) {
        // assemble resistance
        G->assemble_resistance();

        if(!initial_stiffness || counter == 0) {
            // assemble stiffness
            G->assemble_matrix();
            // process loads
            G->process_load();
            // process constraints
            G->process_constraint();
            // call solver
            flag = W->get_stiffness()->solve(get_ninja(W), W->get_trial_load() - W->get_sushi());
        } else
            flag = W->get_stiffness()->solve_trs(get_ninja(W), W->get_trial_load() - W->get_sushi());

        suanpan_debug([&]() {
            if(!W->get_ninja().is_finite()) throw;
        });

        // make sure lapack solver succeeds
        if(flag != 0) return flag;

        // avoid machine error accumulation
        G->erase_machine_error();
        // update trial status for factory
        W->update_trial_displacement(W->get_trial_displacement() + W->get_ninja());
        // update for nodes and elements
        if(G->update_trial_status() != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        // exit if converged
        if(C->is_converged()) return SUANPAN_SUCCESS;
        // exit if maximum iteration is hit
        if(++counter > max_iteration) return SUANPAN_FAIL;
    }
}

void Newton::print() { suanpan_info("A solver using Newton--Raphson iteration method.\n"); }
