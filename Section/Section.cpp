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

#include "Section.h"
#include <Domain/DomainBase.h>
#include <Material/Material.h>

Section::Section(const unsigned T, const unsigned CT, const SectionType ST, const unsigned MT, const double A)
    : Tag(T, CT)
    , material_tag(MT)
    , area(A)
    , section_type(ST) {}

Section::~Section() = default;

void Section::initialize(const shared_ptr<DomainBase>& D) {
    if(material_tag != 0)
        if(!D->find_material(material_tag) || D->get_material(material_tag)->material_type != MaterialType::D1) {
            D->disable_section(get_tag());
            suanpan_error("initialize() cannot find material %u or wrong material type assigned, now disable it.\n", material_tag);
            return;
        }

    const auto size = static_cast<unsigned>(section_type);

    eccentricity.zeros(size - 1);

    if(current_deformation.is_empty()) current_deformation.zeros(size);
    if(trial_deformation.is_empty()) trial_deformation.zeros(size);

    // current_deformation_rate.zeros(size);
    // trial_deformation_rate.zeros(size);

    if(current_resistance.is_empty()) current_resistance.zeros(size);
    if(trial_resistance.is_empty()) trial_resistance.zeros(size);

    if(initial_stiffness.is_empty()) initial_stiffness.zeros(size, size);
    if(trial_stiffness.is_empty()) trial_stiffness.zeros(size, size);
    if(current_stiffness.is_empty()) current_stiffness.zeros(size, size);
}

const vec& Section::get_deformation() const { return trial_deformation; }

const vec& Section::get_deformation_rate() const { return trial_deformation_rate; }

const vec& Section::get_resistance() const { return trial_resistance; }

const mat& Section::get_stiffness() const { return trial_stiffness; }

const mat& Section::get_initial_stiffness() const { return initial_stiffness; }

unique_ptr<Section> Section::get_copy() { throw invalid_argument("hidden method called.\n"); }

double Section::get_parameter(const ParameterType&) { return 0.; }

int Section::update_incre_status(const vec& i_deformation) { return update_trial_status(current_deformation + i_deformation); }

int Section::update_incre_status(const vec& i_deformation, const vec& i_deformation_rate) { return update_trial_status(current_deformation + i_deformation, current_deformation_rate + i_deformation_rate); }

int Section::update_trial_status(const vec&) { throw invalid_argument("hidden method called.\n"); }

int Section::update_trial_status(const vec& t_deformation, const vec&) { return update_trial_status(t_deformation); }

int Section::clear_status() { throw invalid_argument("hidden method called.\n"); }

int Section::commit_status() { throw invalid_argument("hidden method called.\n"); }

int Section::reset_status() { throw invalid_argument("hidden method called.\n"); }

unique_ptr<Section> suanpan::make_copy(const shared_ptr<Section>& S) { return S->get_copy(); }

unique_ptr<Section> suanpan::make_copy(const unique_ptr<Section>& S) { return S->get_copy(); }
