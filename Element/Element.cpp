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

#include "Element.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material.h>
#include <Section/Section.h>
#include <utility>

/**
 * \brief the basic ctor of a typical element
 * \param T element tag
 * \param CT class tag
 * \param NN  number of nodes
 * \param ND  number of dofs
 * \param NT node tags
 * \param MT material tags
 * \param ST section tags
 * \param F nlgeom switch
 */
Element::Element(const unsigned T, const unsigned CT, const unsigned NN, const unsigned ND, uvec NT, uvec MT, uvec ST, const bool F)
    : Tag(T, CT)
    , num_node(NN)
    , num_dof(ND)
    , node_encoding(std::move(NT))
    , material_tag(std::move(MT))
    , section_tag(std::move(ST))
    , nlgeom(F) {
    suanpan_debug("Element %u ctor() called.\n", T);
}

Element::~Element() { suanpan_debug("Element %u dtor() called.\n", get_tag()); }

void Element::initialize(const shared_ptr<DomainBase>& D) {
    // initialized before check node vadality
    if(node_ptr.size() == num_node) {
        for(const auto& I : node_ptr) {
            const auto& t_node = I.lock();
            if(t_node == nullptr || !t_node->is_active()) {
                D->disable_element(get_tag());
                return;
            }
        }
        return;
    }

    // first initiliazation
    const auto total_dof = num_node * num_dof;

    if(total_dof == 0) {
        D->disable_element(get_tag());
        return;
    }

    dof_encoding.set_size(total_dof);

    // check if nodes are still valid
    node_ptr.clear();
    node_ptr.reserve(num_node);
    for(const auto& t_tag : node_encoding) {
        if(!D->find_node(unsigned(t_tag)) || !D->get_node(unsigned(t_tag))->is_active()) {
            suanpan_debug("Element %u finds an invalid node %u, now disable it.\n", get_tag(), t_tag);
            D->disable_element(get_tag());
            return;
        }
        auto& t_node = D->get_node(unsigned(t_tag));
        if(t_node->get_dof_number() < num_dof) t_node->set_dof_number(num_dof);
        node_ptr.emplace_back(t_node);
    }

    // check if material models are valid
    for(const auto& t_material : material_tag) {
        const auto t_tag = unsigned(t_material);
        if(!D->find_material(t_tag) || !D->get_material(t_tag)->is_active()) {
            suanpan_debug("Element %u cannot find valid material %u, now disable it.\n", get_tag(), t_tag);
            D->disable_element(get_tag());
            return;
        }
    }

    // check if section models are valid
    for(const auto& t_section : section_tag) {
        const auto t_tag = unsigned(t_section);
        if(!D->find_section(t_tag) || !D->get_section(t_tag)->is_active()) {
            suanpan_debug("Element %u cannot find valid section %u, now disable it.\n", get_tag(), t_tag);
            D->disable_element(get_tag());
            return;
        }
    }
}

void Element::update_dof_encoding() {
    auto idx = 0;
    for(const auto& tmp_ptr : node_ptr) {
        auto& node_dof = tmp_ptr.lock()->get_reordered_dof();
        for(unsigned i = 0; i < num_dof; ++i) dof_encoding(idx++) = node_dof(i);
    }
}

const unsigned& Element::get_dof_number() const { return num_dof; }

const unsigned& Element::get_node_number() const { return num_node; }

const uvec& Element::get_dof_encoding() const { return dof_encoding; }

const uvec& Element::get_node_encoding() const { return node_encoding; }

const vector<weak_ptr<Node>>& Element::get_node_ptr() const { return node_ptr; }

const vec& Element::get_resistance() const { return trial_resistance; }

const mat& Element::get_mass() const { return trial_mass; }

const mat& Element::get_damping() const { return trial_damping; }

const mat& Element::get_stiffness() const { return trial_stiffness; }

const mat& Element::get_geometry() const { return trial_geometry; }

const mat& Element::get_secant() const { return get_stiffness(); }

const mat& Element::get_initial_mass() const { return initial_mass; }

const mat& Element::get_initial_damping() const { return initial_damping; }

const mat& Element::get_initial_stiffness() const { return initial_stiffness; }

const mat& Element::get_initial_geometry() const { return initial_geometry; }

const mat& Element::get_initial_secant() const { return get_initial_stiffness(); }

Op<vec, op_diagmat> Element::get_diag_mass() const { throw invalid_argument("hidden method called.\n"); }

Op<vec, op_diagmat> Element::get_initial_diag_mass() const { throw invalid_argument("hidden method called.\n"); }

int Element::update_status() { throw invalid_argument("hidden method called.\n"); }

int Element::clear_status() {
    if(!initial_mass.is_empty()) trial_mass = current_mass = initial_mass;
    if(!initial_damping.is_empty()) trial_damping = current_damping = initial_damping;
    if(!initial_stiffness.is_empty()) trial_stiffness = current_stiffness = initial_stiffness;
    if(!initial_geometry.is_empty()) trial_geometry = current_geometry = initial_geometry;

    if(!trial_resistance.is_empty()) trial_resistance.zeros();
    if(!current_resistance.is_empty()) current_resistance.zeros();

    return 0;
}

int Element::commit_status() {
    if(!trial_mass.is_empty()) current_mass = trial_mass;
    if(!trial_damping.is_empty()) current_damping = trial_damping;
    if(!trial_stiffness.is_empty()) current_stiffness = trial_stiffness;
    if(!trial_geometry.is_empty()) current_geometry = trial_geometry;
    if(!trial_resistance.is_empty()) current_resistance = trial_resistance;

    return 0;
}

int Element::reset_status() {
    if(!trial_mass.is_empty()) trial_mass = current_mass;
    if(!trial_damping.is_empty()) trial_damping = current_damping;
    if(!trial_stiffness.is_empty()) trial_stiffness = current_stiffness;
    if(!trial_geometry.is_empty()) trial_geometry = current_geometry;
    if(!trial_resistance.is_empty()) trial_resistance = current_resistance;

    return 0;
}

vector<vec> Element::record(const OutputType&) { return {}; }
