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

#include "Bilinear1D.h"
#include <Toolbox/utility.h>

Bilinear1D::Bilinear1D(const unsigned T, const double E, const double Y, const double H, const double B, const double R)
    : Material1D(T, MT_BILINEAR1D, R)
    , elastic_modulus(E)
    , yield_stress(Y)
    , hardening_ratio(H)
    , beta(B)
    , plastic_modulus(elastic_modulus * hardening_ratio / (1. - hardening_ratio))
    , tolerance(1E-12) {}

void Bilinear1D::initialize(const shared_ptr<DomainBase>&) {
    current_history.zeros(2);
    trial_history.zeros(2);

    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;
}

unique_ptr<Material> Bilinear1D::get_copy() { return make_unique<Bilinear1D>(*this); }

int Bilinear1D::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;
    incre_strain = trial_strain - current_strain;

    if(incre_strain(0) == 0.) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& trial_back_stress = trial_history(0);
    auto& trial_plastic_strain = trial_history(1);

    trial_stiffness = initial_stiffness;

    trial_stress = current_stress + elastic_modulus * incre_strain;

    const auto shifted_stress = trial_stress(0) - trial_back_stress;

    const auto yield_func = abs(shifted_stress) - yield_stress - (1. - beta) * plastic_modulus * trial_plastic_strain;

    if(yield_func > tolerance) {
        const auto incre_plastic_strain = yield_func / (elastic_modulus + plastic_modulus);
        const auto tmp_a = suanpan::sign(shifted_stress);
        trial_stress -= tmp_a * elastic_modulus * incre_plastic_strain;
        trial_stiffness *= hardening_ratio;
        trial_back_stress += tmp_a * beta * plastic_modulus * incre_plastic_strain;
        trial_plastic_strain += incre_plastic_strain;
    }

    // keep stiffness as a nonzero number
    if(abs(trial_stiffness(0)) < 1E-10) trial_stiffness(0) = (suanpan::sign(trial_stiffness(0)) == 0. ? 1. : suanpan::sign(trial_stiffness(0))) * 1E-8;

    return SUANPAN_SUCCESS;
}

int Bilinear1D::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history.zeros();
    current_stiffness = initial_stiffness;
    return reset_status();
}

int Bilinear1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return 0;
}

int Bilinear1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return 0;
}

void Bilinear1D::print() {
    suanpan_info("1D bilinear hardening material using J2 plasticity and associated flow rule.\n");
    suanpan_info("Current Strain: %.3E\tCurrent Stress: %.3E\n", current_strain(0), current_stress(0));
}
