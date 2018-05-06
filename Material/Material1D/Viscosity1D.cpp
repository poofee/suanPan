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

#include "Viscosity1D.h"
#include <Toolbox/utility.h>

double Viscosity1D::compute_damping_coefficient(const double strain, const double strain_rate) const {
    auto eta = damping_c;

    const auto denominator_a = 1. + exp(gap_a * strain);
    const auto denominator_b = 1. + exp(gap_b * strain_rate);

    if(factor_b != 0.) eta += factor_b / denominator_a;
    if(factor_c != 0.) eta += factor_c / denominator_b;
    if(factor_d != 0.) eta += factor_d / denominator_a / denominator_b;

    return eta;
}

Viscosity1D::Viscosity1D(const unsigned T, const double A, const double CA, const double CB, const double CC, const double CD, const double GA, const double GB)
    : Material1D(T, MT_VISCOSITY1D, 0.)
    , alpha(abs(A))
    , damping_a(abs(CA))
    , damping_b(abs(CB))
    , damping_c(abs(CC))
    , damping_d(abs(CD))
    , factor_a(damping_c)
    , factor_b(damping_d - damping_c)
    , factor_c(damping_b - damping_c)
    , factor_d(damping_a + damping_c - damping_b - damping_d)
    , gap_a(-abs(GA))
    , gap_b(-abs(GB)) {}

void Viscosity1D::initialize(const shared_ptr<DomainBase>&) {}

unique_ptr<Material> Viscosity1D::get_copy() { return make_unique<Viscosity1D>(*this); }

int Viscosity1D::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
    trial_strain = t_strain;
    trial_strain_rate = t_strain_rate;

    const auto& u = trial_strain(0);
    const auto& v = trial_strain_rate(0);

    trial_stress = suanpan::sign(v) * compute_damping_coefficient(u, v) * pow(fabs(v), alpha);

    return 0;
}

int Viscosity1D::clear_status() {
    current_strain.zeros();
    current_strain_rate.zeros();
    current_stress.zeros();
    trial_strain.zeros();
    trial_strain_rate.zeros();
    trial_stress.zeros();
    return 0;
}

int Viscosity1D::commit_status() {
    current_strain = trial_strain;
    current_strain_rate = trial_strain_rate;
    current_stress = trial_stress;
    return 0;
}

int Viscosity1D::reset_status() {
    trial_strain = current_strain;
    trial_strain_rate = current_strain_rate;
    trial_stress = current_stress;
    return 0;
}

void Viscosity1D::print() { suanpan_info("1-D Vicosity Material %u.\n", get_tag()); }
