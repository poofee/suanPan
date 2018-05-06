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

#include "RebarLayer.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.hpp>
#include <Toolbox/tensorToolbox.h>
#include <Toolbox/utility.h>
#include <utility>

const unsigned RebarLayer::m_node = 4;
const unsigned RebarLayer::m_dof = 2;
const unsigned RebarLayer::m_size = m_dof * m_node;

RebarLayer::IntegrationPoint::IntegrationPoint(vec C, const double W, const double J, unique_ptr<Material>&& M, mat PNPXY)
    : coor(std::move(C))
    , weight(W)
    , jacob_det(J)
    , m_material(std::move(M))
    , pn_pxy(std::move(PNPXY)) {}

RebarLayer::RebarLayer(const unsigned T, const uvec& N, const unsigned M, const double TH, const bool RS)
    : MaterialElement(T, ET_REBARLAYER, m_node, m_dof, N, uvec{ M })
    , d_pos_a(0)
    , d_pos_b(0)
    , d_rho_a(0)
    , d_rho_b(0)
    , d_area_a(0)
    , d_area_b(0) {}

void RebarLayer::initialize(const shared_ptr<DomainBase>& D) {
    mat ele_coor(m_node, m_dof);
    for(unsigned I = 0; I < m_node; ++I) {
        auto& tmp_coor = node_ptr[I].lock()->get_coordinate();
        for(unsigned J = 0; J < m_dof; ++J) ele_coor(I, J) = tmp_coor(J);
    }

    auto& material_proto = D->get_material(unsigned(material_tag(0)));

    const IntegrationPlan plan(1, 2, IntegrationType::GAUSS);

    int_pt.clear(), int_pt.reserve(2 * plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        const vec t_vec{ d_pos_a, plan(I, 0) };
        const auto pn = shape::quad(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(t_vec, plan(I, 1), det(jacob), material_proto->get_copy(), solve(jacob, pn));
    }
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        const vec t_vec{ d_pos_b, plan(I, 0) };
        const auto pn = shape::quad(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(t_vec, plan(I, 1), det(jacob), material_proto->get_copy(), solve(jacob, pn));
    }
}

int RebarLayer::update_status() {
    auto code = 0;
    code++;
    return code;
}

int RebarLayer::commit_status() { return 0; }

int RebarLayer::clear_status() { return 0; }

int RebarLayer::reset_status() { return 0; }

vector<vec> RebarLayer::record(const OutputType& T) {
    vector<vec> data;
    switch(T) {
    default:
        break;
    }
    return data;
}

void RebarLayer::print() { suanpan_info("Three part rebar layer.\n"); }
