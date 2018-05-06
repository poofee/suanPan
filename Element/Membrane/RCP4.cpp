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

#include "RCP4.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.hpp>
#include <Toolbox/utility.h>
#include <utility>

const unsigned RCP4::m_node = 4;
const unsigned RCP4::m_dof = 2;
const unsigned RCP4::m_size = m_dof * m_node;

RCP4::IntegrationPoint::IntegrationPoint(vec C, const double W, const double J, unique_ptr<Material>&& M, mat PNPXY)
    : coor(std::move(C))
    , weight(W)
    , jacob_det(J)
    , m_material(std::move(M))
    , pn_pxy(std::move(PNPXY))
    , B(3, m_size, fill::zeros) {}

RCP4::RCP4(const unsigned T, const uvec& N, const unsigned M, const double TH, const double RX, const unsigned MX, const double RY, const unsigned MY, const bool R, const bool F)
    : MaterialElement(T, ET_RCP4, m_node, m_dof, N, uvec{ M }, F)
    , thickness(TH)
    , reduced_scheme(R)
    , rho_x(RX)
    , rho_y(RY)
    , mat_x(MX)
    , mat_y(MY) {}

void RCP4::initialize(const shared_ptr<DomainBase>& D) {
    mat ele_coor(m_node, m_dof);
    for(unsigned I = 0; I < m_node; ++I) {
        auto& t_coor = node_ptr[I].lock()->get_coordinate();
        for(unsigned J = 0; J < m_dof; ++J) ele_coor(I, J) = t_coor(J);
    }

    if(reduced_scheme) {
        hourglassing.zeros(m_size, m_size);
        const auto area = .5 * ((ele_coor(2, 0) - ele_coor(0, 0)) * (ele_coor(3, 1) - ele_coor(1, 1)) + (ele_coor(1, 0) - ele_coor(3, 0)) * (ele_coor(2, 1) - ele_coor(0, 1)));
        vec b1(4), b2(4);
        b1(0) = ele_coor(1, 1) - ele_coor(3, 1);
        b1(1) = ele_coor(2, 1) - ele_coor(0, 1);
        b1(2) = ele_coor(3, 1) - ele_coor(1, 1);
        b1(3) = ele_coor(0, 1) - ele_coor(2, 1);
        b2(0) = ele_coor(3, 0) - ele_coor(1, 0);
        b2(1) = ele_coor(0, 0) - ele_coor(2, 0);
        b2(2) = ele_coor(1, 0) - ele_coor(3, 0);
        b2(3) = ele_coor(2, 0) - ele_coor(0, 0);
        const vec h{ std::initializer_list<double>{ 1., -1., 1., -1. } };
        vec gamma = 2. * area * h - dot(h, ele_coor.col(0)) * b1 - dot(h, ele_coor.col(1)) * b2;
        mat t_hourglassing = gamma * gamma.t();
        for(unsigned I = 0; I < m_node; ++I)
            for(unsigned J = 0; J < m_node; ++J) hourglassing(m_dof * I + 1, m_dof * J + 1) = hourglassing(m_dof * I, m_dof * J) = t_hourglassing(I, J);
    }

    auto& material_proto = D->get_material(unsigned(material_tag(0)));

    auto ini_stiffness = material_proto->get_initial_stiffness();

    if(material_proto->material_type == MaterialType::D2 && std::dynamic_pointer_cast<Material2D>(material_proto)->plane_type == PlaneType::E) suanpan::hacker(thickness) = 1.;

    if(reinforced_x) {
        if(D->find_material(mat_x))
            ini_stiffness(0, 0) += D->get_material(mat_x)->get_initial_stiffness().at(0) * rho_x;
        else {
            suanpan_error("cannot find the reinforcement material.\n");
            D->disable_element(get_tag());
            return;
        }
    }
    if(reinforced_y) {
        if(D->find_material(mat_y))
            ini_stiffness(1, 1) += D->get_material(mat_y)->get_initial_stiffness().at(0) * rho_y;
        else {
            suanpan_error("cannot find the reinforcement material.\n");
            D->disable_element(get_tag());
            return;
        }
    }

    const IntegrationPlan plan(2, reduced_scheme ? 1 : 2, IntegrationType::GAUSS);

    initial_stiffness.zeros(m_size, m_size);

    int_pt.clear(), int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        const vec t_vec{ plan(I, 0), plan(I, 1) };
        const auto pn = shape::quad(t_vec, 1);
        const mat jacob = pn * ele_coor;
        const mat pn_pxy = solve(jacob, pn);
        int_pt.emplace_back(t_vec, plan(I, 2), det(jacob), material_proto->get_copy(), pn_pxy);

        if(reinforced_x) int_pt[I].rebar_x = D->get_material(mat_x)->get_copy();
        if(reinforced_y) int_pt[I].rebar_y = D->get_material(mat_y)->get_copy();

        for(unsigned J = 0; J < m_node; ++J) {
            const auto K = m_dof * J;
            int_pt[I].B(0, K) = int_pt[I].B(2, K + 1) = pn_pxy(0, J);
            int_pt[I].B(2, K) = int_pt[I].B(1, K + 1) = pn_pxy(1, J);
        }
        initial_stiffness += int_pt[I].jacob_det * int_pt[I].weight * thickness * int_pt[I].B.t() * ini_stiffness * int_pt[I].B;
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    if(nlgeom)
        for(auto&& I : int_pt) {
            I.BN.zeros(3, m_size);
            I.BG.zeros(m_dof * m_dof, m_size);
            for(unsigned J = 0; J < m_node; ++J) {
                I.BG(0, m_dof * J) = I.BG(2, m_dof * J + 1) = I.pn_pxy(0, J);
                I.BG(1, m_dof * J) = I.BG(3, m_dof * J + 1) = I.pn_pxy(1, J);
            }
        }

    initial_mass.zeros(m_size, m_size);
    const auto tmp_density = material_proto->get_parameter();
    if(tmp_density != 0.) {
        for(const auto& I : int_pt) {
            const auto n_int = shape::quad(I.coor, 0);
            const auto tmp_a = tmp_density * I.jacob_det * I.weight * thickness;
            for(unsigned J = 0; J < m_node; ++J)
                for(auto K = J; K < m_node; ++K) initial_mass(m_dof * J, m_dof * K) += tmp_a * n_int(J) * n_int(K);
        }
        for(unsigned I = 0; I < m_node * m_dof; I += m_dof) {
            initial_mass(I + 1, I + 1) = initial_mass(I, I);
            for(auto J = I + m_dof; J < m_size; J += m_dof) initial_mass(J, I) = initial_mass(I + 1, J + 1) = initial_mass(J + 1, I + 1) = initial_mass(I, J);
        }
    }
    trial_mass = current_mass = initial_mass;
}

int RCP4::update_status() {
    auto code = 0;

    vec t_strain(3);
    mat ele_disp(m_node, m_dof);
    for(unsigned I = 0; I < m_node; ++I) {
        auto& t_disp = node_ptr[I].lock()->get_trial_displacement();
        for(unsigned J = 0; J < m_dof; ++J) ele_disp(I, J) = t_disp(J);
    }

    if(nlgeom) trial_geometry.zeros(m_size, m_size);

    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);
    for(auto& I : int_pt) {
        if(nlgeom) {
            mat gradient = I.pn_pxy * ele_disp;
            gradient(0, 0) += 1., gradient(1, 1) += 1.;
            mat t_mat = .5 * gradient * gradient.t();
            t_strain(0) = t_mat(0, 0) - .5;
            t_strain(1) = t_mat(1, 1) - .5;
            t_strain(2) = t_mat(0, 1) + t_mat(1, 0);
            for(unsigned J = 0; J < m_node; ++J) {
                const auto tmp_a = m_dof * J;
                const auto tmp_b = tmp_a + 1;
                I.BN(0, tmp_a) = I.pn_pxy(0, J) * gradient(0, 0);
                I.BN(1, tmp_a) = I.pn_pxy(1, J) * gradient(1, 0);
                I.BN(2, tmp_a) = I.pn_pxy(0, J) * gradient(1, 0) + I.pn_pxy(1, J) * gradient(0, 0);
                I.BN(0, tmp_b) = I.pn_pxy(0, J) * gradient(0, 1);
                I.BN(1, tmp_b) = I.pn_pxy(1, J) * gradient(1, 1);
                I.BN(2, tmp_b) = I.pn_pxy(0, J) * gradient(1, 1) + I.pn_pxy(1, J) * gradient(0, 1);
            }
        } else {
            t_strain.zeros();
            for(unsigned J = 0; J < m_node; ++J) {
                const auto& t_disp = node_ptr[J].lock()->get_trial_displacement();
                t_strain(0) += t_disp(0) * I.pn_pxy(0, J);
                t_strain(1) += t_disp(1) * I.pn_pxy(1, J);
                t_strain(2) += t_disp(0) * I.pn_pxy(1, J) + t_disp(1) * I.pn_pxy(0, J);
            }
        }

        code += I.m_material->update_trial_status(t_strain);

        const auto t_factor = I.jacob_det * I.weight * thickness;

        auto t_stiff = I.m_material->get_trial_stiffness();
        auto t_stress = I.m_material->get_trial_stress();

        if(reinforced_x) {
            I.rebar_x->update_trial_status(t_strain(0));
            t_stress(0) += I.rebar_x->get_trial_stress().at(0) * rho_x;
            t_stiff(0, 0) += I.rebar_x->get_trial_stiffness().at(0) * rho_x;
        }
        if(reinforced_y) {
            I.rebar_y->update_trial_status(t_strain(1));
            t_stress(1) += I.rebar_y->get_trial_stress().at(0) * rho_y;
            t_stiff(1, 1) += I.rebar_y->get_trial_stiffness().at(0) * rho_y;
        }

        if(nlgeom) {
            mat sigma(4, 4, fill::zeros);
            sigma(0, 0) = sigma(2, 2) = t_stress(0);
            sigma(1, 1) = sigma(3, 3) = t_stress(1);
            sigma(0, 1) = sigma(1, 0) = sigma(2, 3) = sigma(3, 2) = t_stress(2);
            trial_geometry += t_factor * I.BG.t() * sigma * I.BG;
            trial_stiffness += t_factor * I.BN.t() * t_stiff * I.BN;
            trial_resistance += t_factor * I.BN.t() * t_stress;
        } else {
            trial_stiffness += t_factor * I.B.t() * t_stiff * I.B;
            trial_resistance += t_factor * I.B.t() * t_stress;
        }
    }

    if(nlgeom) trial_stiffness += trial_geometry;

    if(reduced_scheme) {
        trial_stiffness += hourglassing;
        trial_resistance += hourglassing * vectorise(ele_disp.t());
    }

    return code;
}

int RCP4::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    if(reinforced_x)
        for(const auto& I : int_pt) code += I.rebar_x->commit_status();
    if(reinforced_y)
        for(const auto& I : int_pt) code += I.rebar_y->commit_status();
    return code;
}

int RCP4::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    if(reinforced_x)
        for(const auto& I : int_pt) code += I.rebar_x->clear_status();
    if(reinforced_y)
        for(const auto& I : int_pt) code += I.rebar_y->clear_status();
    return code;
}

int RCP4::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    if(reinforced_x)
        for(const auto& I : int_pt) code += I.rebar_x->reset_status();
    if(reinforced_y)
        for(const auto& I : int_pt) code += I.rebar_y->reset_status();
    return code;
}

vector<vec> RCP4::record(const OutputType& P) {
    vector<vec> output;
    output.reserve(int_pt.size());

    for(const auto& I : int_pt)
        for(const auto& J : I.m_material->record(P)) output.emplace_back(J);

    return output;
}

void RCP4::print() {
    suanpan_info("Element %u is a four-node membrane element (RCP4)%s.\n", get_tag(), nlgeom ? " with nonlinear geomotry (TL formulation)" : "");
    suanpan_info("The nodes connected are:\n");
    node_encoding.t().print();
    suanpan_info("\nMaterial model response:");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("\nIntegration Point %u:\t", I + 1);
        int_pt[I].coor.t().print();
        int_pt[I].m_material->print();
    }
}
