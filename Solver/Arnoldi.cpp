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

#include "Arnoldi.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Toolbox/arpack_wrapper.h>

Arnoldi::Arnoldi(const unsigned T)
    : Solver(T, CT_ARNOLDI) {}

int Arnoldi::analyze() {
    auto& G = get_integrator();
    auto& W = G->get_domain().lock()->get_factory();

    // assemble resistance
    G->assemble_resistance();
    // assemble stiffness
    G->assemble_matrix();
    // process loads
    G->process_load();
    // process constraints
    G->process_constraint();

    switch(W->get_storage_scheme()) {
    case StorageScheme::BANDSYMM: {
        auto mass = reinterpret_cast<const shared_ptr<BandSymmMat<double>>&>(W->get_mass());
        auto stiffness = reinterpret_cast<const shared_ptr<BandSymmMat<double>>&>(W->get_stiffness());
        eig_solve(get_eigenvalue(W), get_eigenvector(W), stiffness, mass, 10);
        break;
    }
    case StorageScheme::FULL: {
        break;
    }
    case StorageScheme::BAND: {
        break;
    }
    case StorageScheme::SYMMPACK: {
        auto mass = reinterpret_cast<const shared_ptr<SymmPackMat<double>>&>(W->get_mass());
        auto stiffness = reinterpret_cast<const shared_ptr<SymmPackMat<double>>&>(W->get_stiffness());
        eig_solve(get_eigenvalue(W), get_eigenvector(W), stiffness, mass, 10);
        break;
    }
    }

    return 0;
}

void Arnoldi::print() { suanpan_info("A solver using Arnoldi--Raphson iteration method.\n"); }
