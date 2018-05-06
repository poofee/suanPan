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

#include "CP6.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/shapeFunction.hpp>
#include <Toolbox/utility.h>
#include <utility>

const unsigned CP6::m_node = 6;
const unsigned CP6::m_dof = 2;
const unsigned CP6::m_size = m_dof * m_node;

CP6::IntegrationPoint::IntegrationPoint(vec C, const double W, unique_ptr<Material>&& M, mat PNPXY)
    : coor(std::move(C))
    , weight(W)
    , m_material(move(M))
    , pn_pxy(std::move(PNPXY))
    , strain_mat(3, m_size, fill::zeros) {}

CP6::CP6(const unsigned T, const uvec& NT, const unsigned MT, const double TH)
    : MaterialElement(T, ET_CP6, m_node, m_dof, NT, uvec{ MT }, false)
    , thickness(TH) {}

void CP6::initialize(const shared_ptr<DomainBase>& D) {
    mat ele_coor(m_node, m_node);
    for(unsigned I = 0; I < m_node; ++I) {
        auto& t_coor = node_ptr[I].lock()->get_coordinate();
        ele_coor(I, 0) = 1.;
        ele_coor(I, 1) = t_coor(0);
        ele_coor(I, 2) = t_coor(1);
        ele_coor(I, 3) = t_coor(0) * t_coor(1);
        ele_coor(I, 4) = t_coor(0) * t_coor(0);
        ele_coor(I, 5) = t_coor(1) * t_coor(1);
    }

    const mat inv_coor = inv(ele_coor);

    area = .5 * det(ele_coor(span(0, 2), span(0, 2)));

    auto& material_proto = D->get_material(unsigned(material_tag(0)));

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    if(material_proto->material_type == MaterialType::D2 && std::dynamic_pointer_cast<Material2D>(material_proto)->plane_type == PlaneType::E) suanpan::hacker(thickness) = 1.;

    initial_stiffness.zeros(m_size, m_size);

    int_pt.clear(), int_pt.reserve(3);
    for(auto I = 0; I < 3; ++I) {
        const vec coor{ ele_coor(I + 3, 1), ele_coor(I + 3, 2) };
        const mat pn_pxy = shape::triangle(coor, 1) * inv_coor;
        int_pt.emplace_back(coor, area * thickness / 3., material_proto->get_copy(), pn_pxy);

        for(unsigned J = 0; J < m_node; ++J) {
            const auto K = m_dof * J;
            int_pt[I].strain_mat(0, K) = int_pt[I].strain_mat(2, K + 1) = pn_pxy(0, J);
            int_pt[I].strain_mat(2, K) = int_pt[I].strain_mat(1, K + 1) = pn_pxy(1, J);
        }
        initial_stiffness += int_pt[I].weight * int_pt[I].strain_mat.t() * ini_stiffness * int_pt[I].strain_mat;
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    initial_mass.zeros(m_size, m_size);
    const auto t_density = material_proto->get_parameter();
    if(t_density != 0.) {
        for(const auto& I : int_pt) {
            const vec n_int = shape::triangle(I.coor, 0) * inv_coor;
            const auto t_factor = t_density * I.weight;
            for(unsigned J = 0; J < m_node; ++J)
                for(auto K = J; K < m_node; ++K) initial_mass(m_dof * J, m_dof * K) += t_factor * n_int(J) * n_int(K);
        }
        for(unsigned I = 0; I < m_size; I += m_dof) {
            initial_mass(I + 1, I + 1) = initial_mass(I, I);
            for(auto J = I + m_dof; J < m_size; J += m_dof) initial_mass(J, I) = initial_mass(I + 1, J + 1) = initial_mass(J + 1, I + 1) = initial_mass(I, J);
        }
    }
    trial_mass = current_mass = initial_mass;
}

int CP6::update_status() {
    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);

    vec t_strain(3);
    for(const auto& I : int_pt) {
        t_strain.zeros();
        for(unsigned J = 0; J < m_node; ++J) {
            const auto& t_disp = node_ptr[J].lock()->get_trial_displacement();
            t_strain(0) += t_disp(0) * I.pn_pxy(0, J);
            t_strain(1) += t_disp(1) * I.pn_pxy(1, J);
            t_strain(2) += t_disp(0) * I.pn_pxy(1, J) + t_disp(1) * I.pn_pxy(0, J);
        }
        if(I.m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        trial_stiffness += I.weight * I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat;
        trial_resistance += I.weight * I.strain_mat.t() * I.m_material->get_trial_stress();
    }

    return SUANPAN_SUCCESS;
}

int CP6::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int CP6::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int CP6::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

void CP6::print() { suanpan_info("CP6 element.\n"); }
