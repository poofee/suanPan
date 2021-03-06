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

#include "PS.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.hpp>
#include <Toolbox/utility.h>
#include <utility>

const unsigned PS::m_node = 4;
const unsigned PS::m_dof = 2;
const unsigned PS::m_size = m_dof * m_node;

PS::IntegrationPoint::IntegrationPoint(vec C, const double W, const double J, unique_ptr<Material>&& M)
    : coor(std::move(C))
    , weight(W)
    , jacob_det(J)
    , m_material(move(M))
    , BU(3, m_size, fill::zeros)
    , BS(3, 5, fill::zeros) {}

PS::PS(const unsigned T, const uvec& N, const unsigned M, const double TH)
    : MaterialElement(T, ET_PS, m_node, m_dof, N, uvec{ M }, false)
    , thickness(TH) {}

void PS::initialize(const shared_ptr<DomainBase>& D) {
    mat ele_coor(m_node, m_dof);
    for(unsigned I = 0; I < m_node; ++I) {
        auto& tmp_coor = node_ptr[I].lock()->get_coordinate();
        for(unsigned J = 0; J < m_dof; ++J) ele_coor(I, J) = tmp_coor(J);
    }

    mat jacob_center = shape::quad(vec{ 0., 0. }, 1) * ele_coor;

    const auto jacob_a = jacob_center(0, 0) * jacob_center(0, 0);
    const auto jacob_b = jacob_center(1, 0) * jacob_center(1, 0);
    const auto jacob_c = jacob_center(0, 1) * jacob_center(0, 1);
    const auto jacob_d = jacob_center(1, 1) * jacob_center(1, 1);
    const auto jacob_e = jacob_center(0, 0) * jacob_center(0, 1);
    const auto jacob_f = jacob_center(1, 0) * jacob_center(1, 1);

    auto& material_proto = D->get_material(unsigned(material_tag(0)));

    if(material_proto->material_type == MaterialType::D2 && std::dynamic_pointer_cast<Material2D>(material_proto)->plane_type == PlaneType::E) suanpan::hacker(thickness) = 1.;

    const auto& t_stiffness = material_proto->get_initial_stiffness();

    const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

    A.zeros(5, 5), C.zeros(5, 8);
    int_pt.clear(), int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        const vec t_vec{ plan(I, 0), plan(I, 1) };
        const auto pn = shape::quad(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(t_vec, plan(I, 2), det(jacob), material_proto->get_copy());

        mat pn_pxy = solve(jacob, pn);
        for(unsigned J = 0; J < m_node; ++J) {
            int_pt[I].BU(2, 2 * J + 1) = int_pt[I].BU(0, 2 * J) = pn_pxy(0, J);
            int_pt[I].BU(2, 2 * J) = int_pt[I].BU(1, 2 * J + 1) = pn_pxy(1, J);
        }

        for(auto J = 0; J < 3; ++J) int_pt[I].BS(J, J) = 1.;
        int_pt[I].BS(0, 3) = jacob_a * int_pt[I].coor(1);
        int_pt[I].BS(0, 4) = jacob_b * int_pt[I].coor(0);
        int_pt[I].BS(1, 3) = jacob_c * int_pt[I].coor(1);
        int_pt[I].BS(1, 4) = jacob_d * int_pt[I].coor(0);
        int_pt[I].BS(2, 3) = jacob_e * int_pt[I].coor(1);
        int_pt[I].BS(2, 4) = jacob_f * int_pt[I].coor(0);

        const auto t_factor = thickness * int_pt[I].jacob_det * int_pt[I].weight;
        const mat t_mat = int_pt[I].BS.t() * t_factor;
        C += t_mat * int_pt[I].BU;
        A += t_mat * solve(t_stiffness, int_pt[I].BS);
    }

    trial_stiffness = C.t() * solve(A, C);

    trial_mass.zeros(m_size, m_size);
    const auto tmp_density = material_proto->get_parameter(ParameterType::DENSITY);
    if(tmp_density != 0.) {
        for(const auto& I : int_pt) {
            const auto n_int = shape::quad(I.coor, 0);
            const auto tmp_a = tmp_density * I.jacob_det * I.weight * thickness;
            for(unsigned J = 0; J < m_node; ++J)
                for(auto K = J; K < m_node; ++K) trial_mass(m_dof * J, m_dof * K) += tmp_a * n_int(J) * n_int(K);
        }
        for(unsigned I = 0; I < m_size; I += m_dof) {
            trial_mass(I + 1, I + 1) = trial_mass(I, I);
            for(auto J = I + m_dof; J < m_size; J += m_dof) trial_mass(J, I) = trial_mass(I + 1, J + 1) = trial_mass(J + 1, I + 1) = trial_mass(I, J);
        }
    }
}

int PS::update_status() {
    vec trial_disp(m_size);

    auto idx = 0;
    for(const auto& I : node_ptr) {
        auto& tmp_disp = I.lock()->get_trial_displacement();
        for(unsigned J = 0; J < m_dof; ++J) trial_disp(idx++) = tmp_disp(J);
    }

    trial_resistance = trial_stiffness * trial_disp;

    return 0;
}

int PS::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int PS::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int PS::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}
