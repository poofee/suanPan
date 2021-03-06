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

#include "GQ12.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.hpp>
#include <Toolbox/utility.h>
#include <utility>

const unsigned GQ12::m_node = 4;
const unsigned GQ12::m_dof = 3;
const unsigned GQ12::m_size = m_dof * m_node;

GQ12::IntegrationPoint::IntegrationPoint(vec C, const double W, const double J, unique_ptr<Material>&& M)
    : coor(std::move(C))
    , weight(W)
    , jacob_det(J)
    , m_material(move(M))
    , strain_mat(3, m_size, fill::zeros) {}

GQ12::GQ12(const unsigned T, const uvec& N, const unsigned M, const double TH)
    : MaterialElement(T, ET_GQ12, m_node, m_dof, N, uvec{ M }, false)
    , thickness(TH) {}

void GQ12::initialize(const shared_ptr<DomainBase>& D) {
    mat ele_coor(m_node, 2);
    for(unsigned I = 0; I < m_node; ++I) {
        auto& tmp_coor = node_ptr[I].lock()->get_coordinate();
        for(unsigned J = 0; J < 2; ++J) ele_coor(I, J) = tmp_coor(J);
    }

    const auto LX1 = ele_coor(1, 1) - ele_coor(0, 1);
    const auto LX2 = ele_coor(2, 1) - ele_coor(1, 1);
    const auto LX3 = ele_coor(3, 1) - ele_coor(2, 1);
    const auto LX4 = ele_coor(0, 1) - ele_coor(3, 1);
    const auto LY1 = ele_coor(0, 0) - ele_coor(1, 0);
    const auto LY2 = ele_coor(1, 0) - ele_coor(2, 0);
    const auto LY3 = ele_coor(2, 0) - ele_coor(3, 0);
    const auto LY4 = ele_coor(3, 0) - ele_coor(0, 0);

    auto& material_proto = D->get_material(unsigned(material_tag(0)));

    if(material_proto->material_type == MaterialType::D2 && std::dynamic_pointer_cast<Material2D>(material_proto)->plane_type == PlaneType::E) suanpan::hacker(thickness) = 1.;

    const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

    mat pnt(2, 8);

    int_pt.clear(), int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        const vec t_vec{ plan(I, 0), plan(I, 1) };
        const auto pn = shape::quad(t_vec, 1);
        const mat jacob = pn * ele_coor;
        const mat pn_pxy = solve(jacob, pn);
        int_pt.emplace_back(t_vec, plan(I, 2), det(jacob), material_proto->get_copy());

        const auto X = 2. * plan(I, 0);
        const auto Y = 2. * plan(I, 1);

        const auto AA = plan(I, 0) + 1.;
        const auto BB = plan(I, 0) - 1.;
        const auto CC = plan(I, 1) + 1.;
        const auto DD = plan(I, 1) - 1.;

        pnt(0, 0) = DD * (LX4 * CC - LX1 * X);
        pnt(0, 1) = DD * (LX2 * CC + LX1 * X);
        pnt(0, 2) = CC * (LX3 * X - LX2 * DD);
        pnt(0, 3) = -CC * (LX3 * X + LX4 * DD);
        pnt(0, 4) = DD * (LY4 * CC - LY1 * X);
        pnt(0, 5) = DD * (LY2 * CC + LY1 * X);
        pnt(0, 6) = CC * (LY3 * X - LY2 * DD);
        pnt(0, 7) = -CC * (LY3 * X + LY4 * DD);
        pnt(1, 0) = BB * (LX4 * Y - LX1 * AA);
        pnt(1, 1) = AA * (LX1 * BB + LX2 * Y);
        pnt(1, 2) = AA * (LX3 * BB - LX2 * Y);
        pnt(1, 3) = -BB * (LX3 * AA + LX4 * Y);
        pnt(1, 4) = BB * (LY4 * Y - LY1 * AA);
        pnt(1, 5) = AA * (LY1 * BB + LY2 * Y);
        pnt(1, 6) = AA * (LY3 * BB - LY2 * Y);
        pnt(1, 7) = -BB * (LY3 * AA + LY4 * Y);

        const mat pnt_pxy = solve(jacob, pnt / 16.);

        for(unsigned J = 0; J < m_node; ++J) {
            int_pt[I].strain_mat(0, m_dof * J) = pn_pxy(0, J);
            int_pt[I].strain_mat(2, m_dof * J) = pn_pxy(1, J);
            int_pt[I].strain_mat(1, m_dof * J + 1) = pn_pxy(1, J);
            int_pt[I].strain_mat(2, m_dof * J + 1) = pn_pxy(0, J);
            int_pt[I].strain_mat(0, m_dof * J + 2) = pnt_pxy(0, J);
            int_pt[I].strain_mat(1, m_dof * J + 2) = pnt_pxy(1, J + 4);
            int_pt[I].strain_mat(2, m_dof * J + 2) = pnt_pxy(0, J + 4) + pnt_pxy(1, J);
        }
    }

    initial_stiffness.zeros(m_size, m_size);
    for(const auto& I : int_pt) initial_stiffness += I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat * I.jacob_det * I.weight * thickness;
}

int GQ12::update_status() {
    auto code = 0, idx = 0;

    vec trial_disp(m_size);
    for(const auto& I : node_ptr) {
        auto& tmp_disp = I.lock()->get_trial_displacement();
        for(unsigned J = 0; J < m_dof; ++J) trial_disp(idx++) = tmp_disp(J);
    }

    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);
    for(const auto& I : int_pt) {
        code += I.m_material->update_trial_status(I.strain_mat * trial_disp);
        const mat tmp_mat = I.strain_mat.t() * I.jacob_det * I.weight * thickness;
        trial_stiffness += tmp_mat * I.m_material->get_trial_stiffness() * I.strain_mat;
        trial_resistance += tmp_mat * I.m_material->get_trial_stress();
    }

    return code;
}

int GQ12::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int GQ12::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int GQ12::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

void GQ12::print() {
    suanpan_info("Material model response:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("Integration Point %u:\n", I + 1);
        int_pt[I].m_material->print();
    }
}
