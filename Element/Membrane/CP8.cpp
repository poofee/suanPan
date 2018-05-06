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

#include "CP8.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.hpp>
#include <Toolbox/utility.h>
#include <utility>

const unsigned CP8::m_node = 8;
const unsigned CP8::m_dof = 2;
const unsigned CP8::m_size = m_dof * m_node;

CP8::IntegrationPoint::IntegrationPoint(vec C, const double W, const double J, unique_ptr<Material>&& M, mat PNPXY)
    : coor(std::move(C))
    , weight(W)
    , jacob_det(J)
    , m_material(move(M))
    , pn_pxy(std::move(PNPXY))
    , strain_mat(3, m_size, fill::zeros) {}

CP8::CP8(const unsigned T, const uvec& N, const unsigned M, const double TH, const bool R, const bool F)
    : MaterialElement(T, ET_CP8, m_node, m_dof, N, uvec{ M }, F)
    , thickness(TH)
    , reduced_scheme(R) {}

void CP8::initialize(const shared_ptr<DomainBase>& D) {
    mat ele_coor(m_node, m_dof);
    for(unsigned I = 0; I < m_node; ++I) {
        auto& tmp_coor = node_ptr[I].lock()->get_coordinate();
        for(unsigned J = 0; J < m_dof; ++J) ele_coor(I, J) = tmp_coor(J);
    }

    auto& material_proto = D->get_material(unsigned(material_tag(0)));

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    if(material_proto->material_type == MaterialType::D2 && std::dynamic_pointer_cast<Material2D>(material_proto)->plane_type == PlaneType::E) suanpan::hacker(thickness) = 1.;

    const IntegrationPlan plan(2, 2, reduced_scheme ? IntegrationType::GAUSS : IntegrationType::IRONS);

    initial_stiffness.zeros(m_size, m_size);

    int_pt.clear(), int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        const vec t_vec{ plan(I, 0), plan(I, 1) };
        const auto pn = shape::quad(t_vec, 1, 8);
        const mat jacob = pn * ele_coor;
        const mat pn_pxy = solve(jacob, pn);
        int_pt.emplace_back(t_vec, plan(I, 2), det(jacob), material_proto->get_copy(), pn_pxy);

        for(unsigned J = 0; J < m_node; ++J) {
            const auto K = m_dof * J;
            int_pt[I].strain_mat(0, K) = int_pt[I].strain_mat(2, K + 1) = pn_pxy(0, J);
            int_pt[I].strain_mat(2, K) = int_pt[I].strain_mat(1, K + 1) = pn_pxy(1, J);
        }
        initial_stiffness += int_pt[I].jacob_det * int_pt[I].weight * thickness * int_pt[I].strain_mat.t() * ini_stiffness * int_pt[I].strain_mat;
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    initial_mass.zeros(m_size, m_size);
    const auto tmp_density = material_proto->get_parameter();
    if(tmp_density != 0.) {
        for(const auto& I : int_pt) {
            const auto n_int = shape::quad(I.coor, 0, 8);
            const auto tmp_a = tmp_density * I.jacob_det * I.weight * thickness;
            for(unsigned J = 0; J < m_node; ++J)
                for(auto K = J; K < m_node; ++K) initial_mass(m_dof * J, m_dof * K) += tmp_a * n_int(J) * n_int(K);
        }

        for(unsigned I = 0; I < m_size; I += m_dof) {
            initial_mass(I + 1, I + 1) = initial_mass(I, I);
            for(auto J = I + m_dof; J < m_size; J += m_dof) initial_mass(J, I) = initial_mass(I + 1, J + 1) = initial_mass(J + 1, I + 1) = initial_mass(I, J);
        }
    }
    trial_mass = current_mass = initial_mass;
}

int CP8::update_status() {
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

        const auto t_factor = I.jacob_det * I.weight * thickness;

        trial_stiffness += t_factor * I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat;
        trial_resistance += t_factor * I.strain_mat.t() * I.m_material->get_trial_stress();
    }

    return SUANPAN_SUCCESS;
}

int CP8::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int CP8::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int CP8::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

void CP8::print() { suanpan_info("This is a CP8 element.\n"); }
