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

#include "Concrete22.h"
#include <Domain/DomainBase.h>
#include <Toolbox/tensorToolbox.h>
#include <Toolbox/utility.h>

Concrete22::Concrete22(const unsigned T, const double PS, const BackboneType BB, const double SR, const bool CO, const double FE, const bool PO, const bool SD, const double R, const PlaneType PT)
    : Material2D(T, MT_CONCRETE22, PT, R)
    , concrete_major(0, PS, BB, CO, false, TensionType::LINEAR, FE, R)
    , concrete_minor(0, PS, BB, CO, false, TensionType::LINEAR, FE, R)
    , poisson(PO)
    , degrade(SD)
    , shear_retention(SR)
    , peak_strain(-8.5E-4 * pow(abs(PS), .25) + datum::eps * 1E10) {}

void Concrete22::initialize(const shared_ptr<DomainBase>& D) {
    concrete_major.Material::initialize(D);
    concrete_major.initialize(D);
    concrete_minor.Material::initialize(D);
    concrete_minor.initialize(D);

    elastic_modulus = concrete_major.get_parameter(ParameterType::ELASTICMODULUS);

    initial_stiffness.zeros(3, 3);
    initial_stiffness(0, 1) = initial_stiffness(1, 0) = poissons_ratio * (initial_stiffness(0, 0) = initial_stiffness(1, 1) = elastic_modulus / (1. - poissons_ratio * poissons_ratio));
    initial_stiffness(2, 2) = shear_modulus = .5 * (1. - poissons_ratio) * elastic_modulus;

    trial_stiffness = current_stiffness = initial_stiffness;

    poissons_mat.eye(3, 3);
    if(poisson) poissons_mat(0, 1) = poissons_mat(1, 0) = poissons_ratio * (poissons_mat(0, 0) = poissons_mat(1, 1) = 1. / (1. - poissons_ratio * poissons_ratio));

    trial_history = current_history.zeros(2);
}

unique_ptr<Material> Concrete22::get_copy() { return make_unique<Concrete22>(*this); }

double Concrete22::get_parameter(const ParameterType& P) const {
    switch(P) {
    case ParameterType::E:
    case ParameterType::ELASTICMODULUS:
    case ParameterType::YOUNGSMODULUS:
        return elastic_modulus;
    case ParameterType::POISSONSRATIO:
        return poissons_ratio;
    default:
        return 0.;
    }
}

int Concrete22::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;
    trial_history = current_history;

    auto& principal_angle = trial_history(0);
    auto& max_tensile_strain = trial_history(1);

    if(!concrete_major.is_cracked() && !concrete_minor.is_cracked()) principal_angle = transform::strain::angle(trial_strain);

    const auto trans_mat = transform::strain::trans(principal_angle);

    vec principal_strain = trans_mat * trial_strain;
    if(poisson) principal_strain = poissons_mat * principal_strain;

    const auto& strain_11 = principal_strain(0);
    const auto& strain_22 = principal_strain(1);

    if(strain_11 > max_tensile_strain) max_tensile_strain = strain_11;

    // update status
    if(concrete_major.Material::update_trial_status(strain_11) + concrete_minor.Material::update_trial_status(strain_22) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    // collect principal stress components
    vec principal_stress(3);
    principal_stress(0) = concrete_major.get_trial_stress().at(0);
    principal_stress(1) = concrete_minor.get_trial_stress().at(0);
    principal_stress(2) = shear_retention * shear_modulus * principal_strain(2);

    // collect principal stiffness components
    trial_stiffness(0, 0) = concrete_major.get_trial_stiffness().at(0, 0);
    trial_stiffness(1, 1) = concrete_minor.get_trial_stiffness().at(0, 0);
    trial_stiffness(2, 2) = shear_modulus;

    if(degrade) {
        // doi.org/10.1061/(ASCE)0733-9445(1993)119:12(3590)
        // model b equation (30)
        const auto beta = std::min(1., 1. / (.9 - .27 * max_tensile_strain / peak_strain));

        if(strain_11 < 0.) {
            principal_stress(0) *= beta;
            trial_stiffness(0, 0) *= beta;
        }
        if(strain_22 < 0.) {
            principal_stress(1) *= beta;
            trial_stiffness(1, 1) *= beta;
        }
    }

    // transform back to nominal direction
    trial_stress = trans_mat.t() * principal_stress;

    // transform back to nominal direction
    trial_stiffness = trans_mat.t() * diagmat(trial_stiffness) * poissons_mat * trans_mat;

    return SUANPAN_SUCCESS;
}

int Concrete22::clear_status() {
    trial_strain = current_strain.zeros();
    trial_stress = current_stress.zeros();
    trial_history = current_history.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return concrete_major.clear_status() + concrete_minor.clear_status();
}

int Concrete22::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    current_history = trial_history;
    return concrete_major.commit_status() + concrete_minor.commit_status();
}

int Concrete22::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    trial_history = current_history;
    return concrete_major.reset_status() + concrete_minor.reset_status();
}

void Concrete22::print() {
    suanpan_info("Strain: ");
    get_trial_strain().t().print();
    suanpan_info("Stress: ");
    get_trial_stress().t().print();
}
