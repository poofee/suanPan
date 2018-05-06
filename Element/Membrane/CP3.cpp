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

#include "CP3.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/utility.h>

const unsigned CP3::m_node = 3;
const unsigned CP3::m_dof = 2;
const unsigned CP3::m_size = m_dof * m_node;

CP3::CP3(const unsigned T, const uvec& NT, const unsigned MT, const double TH)
    : MaterialElement(T, ET_CP3, m_node, m_dof, NT, uvec{ MT }, false)
    , thickness(TH) {}

void CP3::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get_material(unsigned(material_tag(0)));

    if(material_proto->material_type == MaterialType::D2 && std::dynamic_pointer_cast<Material2D>(material_proto)->plane_type == PlaneType::E) suanpan::hacker(thickness) = 1.;

    m_material = material_proto->get_copy();

    mat ele_coor(m_node, m_node, fill::ones);
    for(unsigned I = 0; I < m_node; ++I) {
        auto& tmp_coor = node_ptr[I].lock()->get_coordinate();
        for(auto J = 0; J < 2; ++J) ele_coor(I, J + 1) = tmp_coor(J);
    }

    area = .5 * det(ele_coor);

    const mat inv_coor = inv(ele_coor);
    pn_pxy = inv_coor.rows(1, 2);

    strain_mat.zeros(3, m_size);
    for(unsigned J = 0; J < m_node; ++J) {
        const auto K = m_dof * J;
        strain_mat(0, K) = strain_mat(2, K + 1) = pn_pxy(0, J);
        strain_mat(2, K) = strain_mat(1, K + 1) = pn_pxy(1, J);
    }
    trial_stiffness = current_stiffness = initial_stiffness = area * thickness * strain_mat.t() * m_material->get_initial_stiffness() * strain_mat;

    initial_mass.zeros(m_size, m_size);
    auto t_density = m_material->get_parameter();
    if(t_density != 0.) {
        t_density *= area * thickness;
        const vec n = mean(ele_coor) * inv_coor;
        for(unsigned I = 0; I < m_node; ++I)
            for(auto J = I; J < m_node; ++J) initial_mass(m_dof * I, m_dof * J) += t_density * n(I) * n(J);
        for(unsigned I = 0; I < m_node * m_dof; I += m_dof) {
            initial_mass(I + 1, I + 1) = initial_mass(I, I);
            for(auto J = I + m_dof; J < m_size; J += m_dof) initial_mass(J, I) = initial_mass(I + 1, J + 1) = initial_mass(J + 1, I + 1) = initial_mass(I, J);
        }
    }
    trial_mass = current_mass = initial_mass;
}

int CP3::update_status() {
    vec t_strain(3, fill::zeros);
    for(unsigned J = 0; J < m_node; ++J) {
        const auto& t_disp = node_ptr[J].lock()->get_trial_displacement();
        t_strain(0) += t_disp(0) * pn_pxy(0, J);
        t_strain(1) += t_disp(1) * pn_pxy(1, J);
        t_strain(2) += t_disp(0) * pn_pxy(1, J) + t_disp(1) * pn_pxy(0, J);
    }

    if(m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    const auto t_factor = area * thickness;

    trial_stiffness = t_factor * strain_mat.t() * m_material->get_trial_stiffness() * strain_mat;
    trial_resistance = t_factor * strain_mat.t() * m_material->get_trial_stress();

    return SUANPAN_SUCCESS;
}

int CP3::commit_status() { return m_material->commit_status(); }

int CP3::clear_status() { return m_material->clear_status(); }

int CP3::reset_status() { return m_material->reset_status(); }

void CP3::print() { suanpan_info("CP3 element.\n"); }
