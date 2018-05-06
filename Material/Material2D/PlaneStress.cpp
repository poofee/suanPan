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

#include "PlaneStress.h"
#include <Domain/DomainBase.h>

const uvec PlaneStress::F1 = uvec{ std::initializer_list<uword>{ 0, 1, 3 } };
const uvec PlaneStress::F2 = uvec{ std::initializer_list<uword>{ 2, 4, 5 } };

mat PlaneStress::form_stiffness(const mat& full_stiffness) {
    mat t_stiffness = full_stiffness(F1, F1);

    if(full_stiffness(2, 2) == 0.) {
        suanpan_error("K(2,2)=0.\n");
        return t_stiffness;
    }

    for(auto I = 0; I < 3; ++I)
        for(auto J = 0; J < 3; ++J) t_stiffness(I, J) -= full_stiffness(F1[I], 2) * full_stiffness(2, F1[J]) / full_stiffness(2, 2);

    return t_stiffness;
}

PlaneStress::PlaneStress(const unsigned T, const unsigned BT, const unsigned MI, const double TL, const bool FM)
    : Material2D(T, MT_PLANESTRESS, PlaneType::S, 0.)
    , base_tag(BT)
    , max_iteration(MI)
    , tolerance(TL)
    , use_full_matrix(FM) {}

PlaneStress::PlaneStress(const PlaneStress& P)
    : Material2D(P.get_tag(), MT_PLANESTRESS, PlaneType::S, P.density)
    , base_tag(P.base_tag)
    , max_iteration(P.max_iteration)
    , tolerance(P.tolerance)
    , use_full_matrix(P.use_full_matrix) {
    if(P.base != nullptr) base = P.base->get_copy();
    Material::initialize();
    PlaneStress::initialize();
}

void PlaneStress::initialize(const shared_ptr<DomainBase>& D) {
    if(D != nullptr) {
        if(D->find_material(base_tag)) {
            auto& material_proto = D->get_material(base_tag);
            if(!material_proto->initialized) {
                material_proto->Material::initialize(D);
                material_proto->initialize(D);
            }
            base = material_proto->get_copy();
        } else {
            D->disable_material(get_tag());
            return;
        }
    }

    if(base == nullptr) {
        disable();
        return;
    }

    current_full_strain.zeros(6);
    trial_full_strain.zeros(6);

    density = base->get_parameter();

    current_stiffness = trial_stiffness = initial_stiffness = form_stiffness(base->get_initial_stiffness());
}

double PlaneStress::get_parameter(const ParameterType& P) const { return base->get_parameter(P); }

unique_ptr<Material> PlaneStress::get_copy() { return make_unique<PlaneStress>(*this); }

int PlaneStress::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;
    trial_full_strain = current_full_strain;

    for(auto I = 0; I < 3; ++I) trial_full_strain(F1[I]) = trial_strain(I);

    unsigned counter = 0;

    if(use_full_matrix) {
        while(++counter < max_iteration) {
            if(base->update_trial_status(trial_full_strain) != 0) return -1;
            const vec incre_full_strain = solve(mat(base->get_trial_stiffness().submat(F2, F2)), base->get_trial_stress().elem(F2));
            trial_full_strain(F2) -= incre_full_strain;
            const auto error = norm(incre_full_strain);
            if(error < tolerance) break;
            suanpan_extra_debug("PlaneStress state determination error: %.4E.\n", error);
        }
    } else {
        const auto& stress_zz = base->get_trial_stress().at(2);
        const auto& stress_yz = base->get_trial_stress().at(4);
        const auto& stress_zx = base->get_trial_stress().at(5);

        while(++counter < max_iteration) {
            if(base->update_trial_status(trial_full_strain) != 0) return -1;
            trial_full_strain(2) -= stress_zz / base->get_trial_stiffness().at(2, 2);
            trial_full_strain(4) -= stress_yz / base->get_trial_stiffness().at(4, 4);
            trial_full_strain(5) -= stress_zx / base->get_trial_stiffness().at(5, 5);
            const auto error = fabs(stress_zz) + fabs(stress_yz) + fabs(stress_zx);
            if(error < tolerance) break;
            suanpan_extra_debug("PlaneStress state determination error: %.4E.\n", error);
        }
    }

    if(counter == max_iteration) {
        suanpan_error("PlaneStress cannot converge with in %u iterations.\n", counter);
        return -1;
    }

    for(auto I = 0; I < 3; ++I) trial_stress(I) = base->get_trial_stress().at(F1[I]);

    trial_stiffness = form_stiffness(base->get_trial_stiffness());

    return 0;
}

int PlaneStress::clear_status() {
    current_full_strain.zeros(6);
    trial_full_strain.zeros(6);
    Material::clear_status();
    return base->clear_status();
}

int PlaneStress::commit_status() {
    current_full_strain = trial_full_strain;
    Material::commit_status();
    return base->commit_status();
}

int PlaneStress::reset_status() {
    trial_full_strain = current_full_strain;
    Material::reset_status();
    return base->reset_status();
}

vector<vec> PlaneStress::record(const OutputType& P) { return base->record(P); }

void PlaneStress::print() {
    suanpan_info("Strain:\n");
    current_strain.t().print();
    suanpan_info("Stress:\n");
    current_stress.t().print();
}
