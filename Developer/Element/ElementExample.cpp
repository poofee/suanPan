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

#include "ElementExample.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material.h>

SUANPAN_EXPORT void new_elementexample(unique_ptr<Element>& return_obj, std::istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_error("new_elementexample() needs a tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 3; ++I) {
        if((command >> node).fail()) {
            suanpan_error("new_elementexample() needs 3 nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if((command >> material_tag).fail()) {
        suanpan_error("new_elementexample() needs a material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof())
        suanpan_info("new_elementexample() assumes a unit thickness.\n");
    else if((command >> thickness).fail()) {
        suanpan_error("new_elementexample() needs a valid thickness.\n");
        return;
    }

    return_obj = make_unique<ElementExample>(tag, uvec(node_tag), material_tag, thickness);
}

const unsigned ElementExample::m_node = 3;
const unsigned ElementExample::m_dof = 2;

#ifndef ET_ELEMENTTEMPLATE
#define ET_ELEMENTTEMPLATE 0
#endif

ElementExample::ElementExample(const unsigned& T, const uvec& NT, const unsigned& MT, const double& TH)
    : Element(T, ET_ELEMENTTEMPLATE, m_node, m_dof, NT, uvec{ MT }, {}, false)
    , thickness(TH) {}

void ElementExample::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get_material(unsigned(material_tag(0)));
    m_material = material_proto->get_copy();

    mat ele_coor(m_node, m_node, fill::ones);
    for(auto I = 0; I < m_node; ++I) {
        auto& tmp_coor = node_ptr.at(I).lock()->get_coordinate();
        for(auto J = 0; J < 2; ++J) ele_coor(I, J + 1) = tmp_coor(J);
    }

    area = det(ele_coor) / 2.;

    mat inv_coor = inv(ele_coor);

    strain_mat.zeros(3, m_node * m_dof);
    for(auto I = 0; I < 3; ++I) {
        strain_mat(0, 2 * I) = inv_coor(1, I);
        strain_mat(1, 2 * I + 1) = inv_coor(2, I);
        strain_mat(2, 2 * I) = inv_coor(2, I);
        strain_mat(2, 2 * I + 1) = inv_coor(1, I);
    }

    const auto tmp_density = m_material->get_parameter();
    if(tmp_density != 0.) {
        const vec n = mean(ele_coor) * inv_coor;
        trial_mass = n * n.t() * tmp_density * area * thickness;
    }
}

int ElementExample::update_status() {
    vec trial_disp(m_node * m_dof);
    auto idx = 0;
    for(const auto& I : node_ptr) {
        auto& tmp_disp = I.lock()->get_trial_displacement();
        for(auto J = 0; J < m_dof; ++J) trial_disp(idx++) = tmp_disp(J);
    }

    m_material->update_trial_status(strain_mat * trial_disp);

    trial_stiffness = strain_mat.t() * m_material->get_trial_stiffness() * strain_mat * area * thickness;
    trial_resistance = strain_mat.t() * m_material->get_trial_stress() * area * thickness;

    return 0;
}

int ElementExample::commit_status() { return m_material->commit_status(); }

int ElementExample::clear_status() { return m_material->clear_status(); }

int ElementExample::reset_status() { return m_material->reset_status(); }

void ElementExample::print() { suanpan_info("This is an element example based on CPS3 element.\n"); }
