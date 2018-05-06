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

#include "QE2.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.hpp>
#include <Toolbox/tensorToolbox.h>
#include <Toolbox/utility.h>
#include <utility>

const unsigned QE2::m_node = 4;
const unsigned QE2::m_dof = 2;
const unsigned QE2::m_size = m_dof * m_node;
mat QE2::mapping;

QE2::IntegrationPoint::IntegrationPoint(vec C, const double F, unique_ptr<Material>&& M)
    : coor(std::move(C))
    , factor(F)
    , m_material(move(M))
    , B(3, m_size, fill::zeros)
    , BI(3, 2, fill::zeros) {}

QE2::QE2(const unsigned T, const uvec& N, const unsigned M, const double TH)
    : MaterialElement(T, ET_QE2, m_node, m_dof, N, uvec{ M })
    , thickness(TH) {}

void QE2::initialize(const shared_ptr<DomainBase>& D) {
    if(mapping.is_empty()) {
        mat t_mapping(4, 4);
        t_mapping.fill(.25);
        t_mapping(1, 0) = t_mapping(1, 3) = t_mapping(2, 0) = t_mapping(2, 1) = t_mapping(3, 1) = t_mapping(3, 3) = -.25;
        mapping = t_mapping;
    }

    mat ele_coor(m_node, 2);
    for(unsigned I = 0; I < m_node; ++I) {
        auto& tmp_coor = node_ptr[I].lock()->get_coordinate();
        for(auto J = 0; J < 2; ++J) ele_coor(I, J) = tmp_coor(J);
    }

    trans_mat = trans(mapping * ele_coor);

    auto& material_proto = D->get_material(unsigned(material_tag(0)));

    if(material_proto->material_type == MaterialType::D2 && std::dynamic_pointer_cast<Material2D>(material_proto)->plane_type == PlaneType::E) suanpan::hacker(thickness) = 1.;

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

    vec disp_mode(4);
    disp_mode(0) = 1.;

    mat H = zeros(7, 7);

    L = zeros(7, 8), LI = zeros(7, 2);

    int_pt.clear(), int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        const vec t_vec{ plan(I, 0), plan(I, 1) };
        const auto pn = shape::quad(t_vec, 1);
        const mat jacob = pn * ele_coor;
        const auto det_jacob = det(jacob);
        int_pt.emplace_back(t_vec, plan(I, 2) * det_jacob * thickness, material_proto->get_copy());

        disp_mode(1) = int_pt[I].coor(0);
        disp_mode(2) = int_pt[I].coor(1);
        disp_mode(3) = int_pt[I].coor(0) * int_pt[I].coor(1);

        int_pt[I].P = shape::stress7(trans_mat * disp_mode);

        int_pt[I].A = solve(ini_stiffness, int_pt[I].P);

        const mat pn_pxy = solve(jacob, pn);
        for(unsigned K = 0; K < m_node; ++K) {
            int_pt[I].B(0, m_dof * K) = int_pt[I].B(2, m_dof * K + 1) = pn_pxy(0, K);
            int_pt[I].B(2, m_dof * K) = int_pt[I].B(1, m_dof * K + 1) = pn_pxy(1, K);
        }

        const vec pe_pxy = int_pt[I].coor / det_jacob;
        int_pt[I].BI(2, 1) = int_pt[I].BI(0, 0) = -trans_mat(1, 2) * pe_pxy(0) - trans_mat(1, 1) * pe_pxy(1);
        int_pt[I].BI(2, 0) = int_pt[I].BI(1, 1) = trans_mat(0, 2) * pe_pxy(0) + trans_mat(0, 1) * pe_pxy(1);

        const mat tmp_mat = int_pt[I].P.t() * int_pt[I].factor;
        H += tmp_mat * int_pt[I].A;
        L += tmp_mat * int_pt[I].B;
        LI += tmp_mat * int_pt[I].BI;
    }

    initial_mass.zeros(m_size, m_size);
    const auto tmp_density = material_proto->get_parameter();
    if(tmp_density != 0.) {
        for(const auto& I : int_pt) {
            const auto n_int = shape::quad(I.coor, 0);
            const auto tmp_a = tmp_density * I.factor;
            for(unsigned J = 0; J < m_node; ++J)
                for(auto K = J; K < m_node; ++K) initial_mass(m_dof * J, m_dof * K) += tmp_a * n_int(J) * n_int(K);
        }
        for(unsigned I = 0; I < m_node * m_dof; I += m_dof) {
            initial_mass(I + 1, I + 1) = initial_mass(I, I);
            for(auto J = I + m_dof; J < m_node * m_dof; J += m_dof) initial_mass(J, I) = initial_mass(I + 1, J + 1) = initial_mass(J + 1, I + 1) = initial_mass(I, J);
        }
    }
    trial_mass = current_mass = initial_mass;

    HT = trans(H);

    solve(HIL, H, L), solve(HILI, H, LI);

    const mat TT = LI.t() * HIL;

    trial_stiffness = current_stiffness = initial_stiffness = HIL.t() * L - TT.t() * (trial_qtitt = current_qtitt = initial_qtitt = solve(mat(LI.t() * HILI), TT));

    pre_disp.zeros(m_size);

    current_qtifi.zeros(2);
    current_lambda.zeros(2);
    current_alpha.zeros(7);
    current_beta.zeros(7);

    trial_qtifi.zeros(2);
    trial_lambda.zeros(2);
    trial_alpha.zeros(7);
    trial_beta.zeros(7);
}

int QE2::update_status() {
    auto idx = 0;
    vec new_disp(m_size);
    for(const auto& t_ptr : node_ptr) {
        auto& t_disp = t_ptr.lock()->get_incre_displacement();
        for(unsigned pos = 0; pos < m_dof; ++pos) new_disp(idx++) = t_disp(pos);
    }

    const vec incre_disp = new_disp - pre_disp;
    const vec incre_lambda = -trial_qtitt * incre_disp - trial_qtifi; // eq. 65
    const vec incre_alpha = HIL * incre_disp + HILI * incre_lambda;   // eq. 57

    trial_lambda += incre_lambda; // eq. 66
    trial_alpha += incre_alpha;   // eq. 46

    vec local_stress(7, fill::zeros);
    mat local_stiffness(7, 7, fill::zeros);
    for(const auto& t_pt : int_pt) {
        if(t_pt.m_material->update_trial_status(t_pt.A * trial_alpha) != SUANPAN_SUCCESS) return -1;
        local_stiffness += t_pt.A.t() * (t_pt.m_material->get_trial_stiffness() * t_pt.factor) * t_pt.A; // eq. 56
        local_stress += t_pt.A.t() * (t_pt.m_material->get_trial_stress() * t_pt.factor);
    }

    const mat HILIHT = HILI.t() * local_stiffness, QT = HILIHT * HILI, TT = HILIHT * HIL; // eq. 60

    solve(trial_beta, HT, local_stress);
    solve(trial_qtitt, QT, TT);                  // eq. 65
    solve(trial_qtifi, QT, LI.t() * trial_beta); // eq. 65

    trial_resistance = L.t() * trial_beta - TT.t() * trial_qtifi;             // eq. 64
    trial_stiffness = HIL.t() * local_stiffness * HIL - TT.t() * trial_qtitt; // eq. 61

    pre_disp = new_disp;

    return 0;
}

int QE2::clear_status() {
    current_lambda.zeros();
    trial_lambda.zeros();
    current_alpha.zeros();
    trial_alpha.zeros();
    current_beta.zeros();
    trial_beta.zeros();
    current_qtifi.zeros();
    trial_qtifi.zeros();
    pre_disp.zeros();

    current_qtitt = trial_qtitt = initial_qtitt;

    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int QE2::commit_status() {
    current_lambda = trial_lambda;
    current_alpha = trial_alpha;
    current_beta = trial_beta;
    current_qtifi = trial_qtifi;
    current_qtitt = trial_qtitt;
    pre_disp.zeros();

    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int QE2::reset_status() {
    trial_lambda = current_lambda;
    trial_alpha = current_alpha;
    trial_beta = current_beta;
    trial_qtifi = current_qtifi;
    trial_qtitt = current_qtitt;
    pre_disp.zeros();

    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

vector<vec> QE2::record(const OutputType& T) {
    vector<vec> data;
    switch(T) {
    case OutputType::E:
        for(const auto& I : int_pt) data.emplace_back(I.A * current_alpha);
        break;
    case OutputType::S:
        for(const auto& I : int_pt) data.emplace_back(I.P * current_beta);
        break;
    default:
        break;
    }
    return data;
}

void QE2::print() {
    suanpan_info("Piltner's mixed quad element %u connects nodes:\n", get_tag());
    node_encoding.t().print();
    suanpan_info("Material model response:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("Integration Point %u:\n", I + 1);
        int_pt[I].m_material->print();
    }
    suanpan_info("Element model response:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("Integration Point %u:\n", I + 1);
        suanpan_info("Strain:\n");
        (int_pt[I].A * current_alpha).t().print();
        suanpan_info("Stress:\n");
        (int_pt[I].P * current_beta).t().print();
    }
}
