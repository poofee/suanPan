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

#include "S4.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.hpp>
#include <Toolbox/tensorToolbox.h>
#include <Toolbox/utility.h>
#include <utility>

const unsigned S4::s_node = 4;
const unsigned S4::s_dof = 6;
const unsigned S4::s_size = s_dof * s_node;
const uvec S4::m_dof = uvec(std::initializer_list<uword>{ 0, 1 });
const uvec S4::p_dof = uvec(std::initializer_list<uword>{ 2, 3, 4 });
mat S4::iso_mapping;

S4::IntegrationPoint::IntegrationPoint(vec C, const double F, unique_ptr<Material>&& M)
    : coor(std::move(C))
    , factor(F)
    , s_material(std::move(M))
    , BM(3, 8, fill::zeros)
    , BP(3, 12, fill::zeros) {}

mat S4::direction_cosine() {
    const auto& coor_i = node_ptr[0].lock()->get_coordinate();
    const auto& coor_j = node_ptr[1].lock()->get_coordinate();
    const auto& coor_k = node_ptr[2].lock()->get_coordinate();

    const vec vector_a = coor_j(span(0, 2)) - coor_i(span(0, 2));
    const vec vector_b = coor_k(span(0, 2)) - coor_i(span(0, 2));

    mat trans(3, 3);

    trans.col(0) = normalise(vector_a);
    trans.col(2) = normalise(cross(vector_a, vector_b));
    trans.col(1) = cross(trans.col(2), trans.col(0));

    return trans;
}

void S4::transform_from_local_to_global(vec& resistance, const mat& trans) {
    for(unsigned I = 0; I < 8; ++I) {
        const auto J = 3 * I;
        auto resistance_part = resistance(span(J, J + 2));
        resistance_part = trans * resistance_part;
    }
}

void S4::transform_from_local_to_global(mat& stiffness, const mat& trans) {
    for(unsigned I = 0; I < 8; ++I) {
        const auto K = 3 * I;
        for(unsigned J = 0; J < 8; ++J) {
            const auto L = 3 * J;
            auto stiffness_part = stiffness(span(K, K + 2), span(L, L + 2));
            stiffness_part = trans * stiffness_part * trans.t();
        }
    }
}

vec S4::reshuffle(const vec& membrane_resistance, const vec& plate_resistance) {
    vec total_resistance(s_size, fill::zeros);

    for(unsigned I = 0; I < s_node; ++I) {
        const auto J = I * s_dof, K = 2 * I, L = 3 * I;
        total_resistance(J + m_dof) = membrane_resistance(span(K, K + 1));
        total_resistance(J + p_dof) = plate_resistance(span(L, L + 2));
    }

    return total_resistance;
}

mat S4::reshuffle(const mat& membrane_stiffness, const mat& plate_stiffness) {
    mat total_stiffness(s_size, s_size, fill::zeros);

    for(unsigned I = 0; I < s_node; ++I) {
        const auto K = I * s_dof, M = 2 * I, N = 3 * I;
        for(unsigned J = 0; J < s_node; ++J) {
            const auto L = J * s_dof, O = 2 * J, P = 3 * J;
            total_stiffness(K + m_dof, L + m_dof) = membrane_stiffness(span(M, M + 1), span(O, O + 1));
            total_stiffness(K + p_dof, L + p_dof) = plate_stiffness(span(N, N + 2), span(P, P + 2));
            total_stiffness(K + 3, L + 3) = 1.;
        }
    }

    return total_stiffness;
}

S4::S4(const unsigned T, const uvec& N, const unsigned M, const double TH, const char IS)
    : MaterialElement(T, ET_S4, s_node, s_dof, N, uvec{ M })
    , thickness(TH)
    , int_scheme(IS) {}

void S4::initialize(const shared_ptr<DomainBase>& D) {
    trans_mat = direction_cosine();

    auto idx = 0;

    auto z_norm = 0.;

    mat ele_coor(4, 2);
    for(const auto& I : node_ptr) {
        const vec local_coor = trans_mat * I.lock()->get_coordinate()(span(0, 2));
        ele_coor.row(idx++) = local_coor(span(0, 1)).t();
        z_norm += local_coor(2);
    }

    if(z_norm > 1E-6) suanpan_warning("non-planar shell geomtry detected.\n");

    // in-plane
    const IntegrationPlan m_plan(2, 3, IntegrationType::GAUSS);
    // along thickness
    const IntegrationPlan t_plan(1, 5, IntegrationType::GAUSS);

    const auto& material_proto = D->get_material(unsigned(material_tag(0)));

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    idx = 0;

    mat m_stiffness(8, 8, fill::zeros), p_stiffness(12, 12, fill::zeros);

    int_pt.clear(), int_pt.reserve(t_plan.n_rows * m_plan.n_rows);
    for(unsigned I = 0; I < m_plan.n_rows; ++I) {
        const vec t_vec{ m_plan(I, 0), m_plan(I, 1) };
        const auto pnm = shape::quad(t_vec, 1);
        const mat jacob = pnm * ele_coor;
        const mat pnm_pxy = solve(jacob, pnm);
        auto pnp = shape::plate(t_vec, 1);
        for(auto J = 0; J < 12; ++J) {
            const auto K = span(2 * J, 2 * J + 1);
            auto pnp_part = pnp(K, K);
            pnp_part = solve(jacob, pnp_part) * inv(jacob).t();
        }
        for(unsigned J = 0; J < t_plan.n_rows; ++J) {
            int_pt.emplace_back(vec(std::initializer_list<double>{ m_plan(I, 0), m_plan(I, 1), .5 * thickness * t_plan(J, 0) }), .5 * thickness * det(jacob) * m_plan(I, 2) * t_plan(J, 1), material_proto->get_copy());
            for(unsigned K = 0; K < s_node; ++K) {
                const auto L = s_dof * K;
                int_pt[idx].BM(0, L) = int_pt[idx].BM(2, L + 1) = pnm_pxy(0, K);
                int_pt[idx].BM(2, L) = int_pt[idx].BM(1, L + 1) = pnm_pxy(1, K);
                int_pt[idx].BP(0, L) = int_pt[idx].BP(2, L + 1) = pnp(0, K);
                int_pt[idx].BP(2, L) = int_pt[idx].BP(1, L + 1) = pnp(1, K);
            }
            m_stiffness += int_pt[I].factor * int_pt[I].BM.t() * ini_stiffness * int_pt[I].BM;
        }
    }
}

int S4::update_status() {
    // seperate displacement vector
    auto idx = 0;
    vec m_disp(8), p_disp(12);
    for(const auto& t_ptr : node_ptr) {
        auto& t_disp = t_ptr.lock()->get_incre_displacement();
        const vec tm_disp = trans_mat * t_disp(span(0, 2));
        const vec tp_disp = trans_mat * t_disp(span(3, 5));
        const auto I = 2 * idx, J = 3 * idx++;
        m_disp(I) = tm_disp(0);
        m_disp(I + 1) = tm_disp(1);
        p_disp(J) = tm_disp(2);
        p_disp(J + 1) = tp_disp(1);
        p_disp(J + 2) = tp_disp(2);
    }

    vec m_resistance(8, fill::zeros), p_resistance(12, fill::zeros);
    mat m_stiffness(8, 8, fill::zeros), p_stiffness(12, 12, fill::zeros);

    for(const auto& I : int_pt) {
        const vec t_strain = I.BM * m_disp + I.coor(2) * I.BP * p_disp;
        if(!I.s_material->update_trial_status(t_strain) == SUANPAN_SUCCESS) return SUANPAN_FAIL;
        auto& t_stress = I.s_material->get_trial_stress();
        auto& t_stiffness = I.s_material->get_trial_stiffness();
        m_resistance += I.factor * I.BM.t() * t_stress;
        p_resistance += I.coor(2) * I.factor * I.BP.t() * t_stress;
        m_stiffness += I.factor * I.BM.t() * t_stiffness * I.BM;
        p_stiffness += I.coor(2) * I.coor(2) * I.factor * I.BP.t() * t_stiffness * I.BP;
    }

    trial_resistance = reshuffle(m_resistance, p_resistance);
    trial_stiffness = reshuffle(m_stiffness, p_stiffness);

    transform_from_local_to_global(trial_resistance, trans_mat);
    transform_from_local_to_global(trial_stiffness, trans_mat);

    return SUANPAN_SUCCESS;
}

int S4::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->commit_status();
    return code;
}

int S4::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->clear_status();
    return code;
}

int S4::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->reset_status();
    return code;
}

vector<vec> S4::record(const OutputType&) {
    vector<vec> data;

    return data;
}

void S4::print() {}
