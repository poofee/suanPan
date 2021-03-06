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

#include "RC01.h"
#include <Domain/DomainBase.h>

RC01::RC01(const unsigned T, const unsigned ST, const unsigned CT)
    : Material2D(T, MT_RC01)
    , rebar_tag(ST)
    , concrete_tag(CT) {}

RC01::RC01(const RC01& P)
    : Material2D(P.get_tag(), MT_RC01, PlaneType::S, P.density)
    , rebar_tag(P.rebar_tag)
    , concrete_tag(P.concrete_tag) {
    if(P.rebar != nullptr) {
        rebar = P.rebar->get_copy();
        rebar->Material::initialize();
        rebar->initialize();
    }
    if(P.concrete != nullptr) {
        concrete = P.concrete->get_copy();
        concrete->Material::initialize();
        concrete->initialize();
    }
    Material::initialize();
    RC01::initialize();
}

void RC01::initialize(const shared_ptr<DomainBase>& D) {
    if(D != nullptr) {
        if(!D->find_material(rebar_tag) || !D->find_material(concrete_tag)) {
            D->disable_material(get_tag());
            return;
        }

        auto& rebar_proto = D->get_material(rebar_tag);
        rebar_proto->Material::initialize(D);
        rebar_proto->initialize(D);
        rebar = rebar_proto->get_copy();
        rebar->Material::initialize(D);
        rebar->initialize(D);

        auto& concrete_proto = D->get_material(concrete_tag);
        concrete_proto->Material::initialize(D);
        concrete_proto->initialize(D);
        concrete = concrete_proto->get_copy();
        concrete->Material::initialize(D);
        concrete->initialize(D);
    }

    density = concrete->get_parameter(ParameterType::DENSITY) + rebar->get_parameter(ParameterType::DENSITY);

    initial_stiffness = rebar->get_initial_stiffness() + concrete->get_initial_stiffness();
    current_stiffness = initial_stiffness;
    trial_stiffness = initial_stiffness;
}

unique_ptr<Material> RC01::get_copy() { return make_unique<RC01>(*this); }

int RC01::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;

    auto code = 0;
    code += rebar->update_trial_status(trial_strain);
    code += concrete->update_trial_status(trial_strain);
    if(code != 0) return code;

    trial_stress = rebar->get_trial_stress() + concrete->get_trial_stress();
    trial_stiffness = rebar->get_trial_stiffness() + concrete->get_trial_stiffness();

    return 0;
}

int RC01::clear_status() {
    current_strain.zeros();
    trial_strain.zeros();
    current_stress.zeros();
    trial_stress.zeros();
    current_stiffness = initial_stiffness;
    trial_stiffness = initial_stiffness;
    auto code = rebar->clear_status();
    code += concrete->clear_status();
    return code;
}

int RC01::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    auto code = rebar->commit_status();
    code += concrete->commit_status();
    return code;
}

int RC01::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    auto code = rebar->reset_status();
    code += concrete->reset_status();
    return code;
}
