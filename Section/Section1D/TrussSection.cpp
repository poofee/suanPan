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

#include "TrussSection.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

TrussSection::TrussSection(const unsigned T, const double B, const unsigned M)
    : Section1D(T, ST_TRUSSSECTION, M, B) {}

TrussSection::TrussSection(const TrussSection& old_obj)
    : Section1D(old_obj.get_tag(), ST_TRUSSSECTION, old_obj.material_tag)
    , s_material(old_obj.s_material->get_copy()) {}

void TrussSection::initialize(const shared_ptr<DomainBase>& D) {
    s_material = D->get_material(material_tag)->get_copy();

    trial_stiffness = current_stiffness = initial_stiffness = area * s_material->get_initial_stiffness().at(0);
}

unique_ptr<Section> TrussSection::get_copy() { return make_unique<TrussSection>(*this); }

double TrussSection::get_parameter(const ParameterType& P) {
    switch(P) {
    case ParameterType::AREA:
        return area;
    case ParameterType::DENSITY:
        return s_material->get_parameter(ParameterType::DENSITY);
    default:
        return 0.;
    }
}

int TrussSection::update_trial_status(const vec& t_deformation) {
    trial_deformation = t_deformation;

    if(s_material->update_trial_status(t_deformation) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    trial_stiffness = s_material->get_trial_stiffness() * area;

    trial_resistance = s_material->get_trial_stress();

    return SUANPAN_SUCCESS;
}

int TrussSection::clear_status() {
    current_deformation.zeros();
    trial_deformation.zeros();
    current_resistance.zeros();
    trial_resistance.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return s_material->clear_status();
}

int TrussSection::commit_status() {
    current_deformation = trial_deformation;
    current_resistance = trial_resistance;
    current_stiffness = trial_stiffness;
    return s_material->commit_status();
}

int TrussSection::reset_status() {
    trial_deformation = current_deformation;
    trial_resistance = current_resistance;
    trial_stiffness = current_stiffness;
    return s_material->reset_status();
}

void TrussSection::print() { suanpan_info("A TrussSection Section.\n"); }
