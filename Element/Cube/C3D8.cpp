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

#include "C3D8.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material3D/Material3D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.hpp>
#include <utility>

const unsigned C3D8::c_node = 8;
const unsigned C3D8::c_dof = 3;
const unsigned C3D8::c_size = c_dof * c_node;

C3D8::IntegrationPoint::IntegrationPoint(vec C, const double W, const double J, unique_ptr<Material>&& M)
    : coor(std::move(C))
    , weight(W)
    , jacob_det(J)
    , c_material(move(M))
    , strain_mat(6, c_size, fill::zeros) {}

C3D8::C3D8(const unsigned T, const uvec& N, const unsigned M, const char R, const bool F)
    : MaterialElement(T, ET_C3D8, c_node, c_dof, N, uvec{ M }, F)
    , int_scheme(R) {}

void C3D8::initialize(const shared_ptr<DomainBase>& D) {
    mat ele_coor(c_node, c_dof);
    for(unsigned I = 0; I < c_node; ++I) {
        auto& t_coor = node_ptr[I].lock()->get_coordinate();
        for(unsigned J = 0; J < c_dof; ++J) ele_coor(I, J) = t_coor(J);
    }

    const auto& material_proto = D->get_material(unsigned(material_tag(0)));

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    const IntegrationPlan plan(3, int_scheme == 'R' ? 1 : 2, int_scheme == 'I' ? IntegrationType::IRONS : IntegrationType::GAUSS);

    initial_stiffness.zeros(c_size, c_size);

    int_pt.clear(), int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        const vec t_vec(std::initializer_list<double>{ plan(I, 0), plan(I, 1), plan(I, 2) });
        const auto pn = shape::cube(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(t_vec, plan(I, c_dof), det(jacob), material_proto->get_copy());

        const mat pn_pxy = solve(jacob, pn);
        for(unsigned J = 0; J < c_node; ++J) {
            const auto K = c_dof * J;
            int_pt[I].strain_mat(0, K) = int_pt[I].strain_mat(3, K + 1) = int_pt[I].strain_mat(5, K + 2) = pn_pxy(0, J);
            int_pt[I].strain_mat(3, K) = int_pt[I].strain_mat(1, K + 1) = int_pt[I].strain_mat(4, K + 2) = pn_pxy(1, J);
            int_pt[I].strain_mat(5, K) = int_pt[I].strain_mat(4, K + 1) = int_pt[I].strain_mat(2, K + 2) = pn_pxy(2, J);
        }
        initial_stiffness += int_pt[I].weight * int_pt[I].jacob_det * int_pt[I].strain_mat.t() * ini_stiffness * int_pt[I].strain_mat;
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    initial_mass.zeros(c_size, c_size);
    const auto t_density = material_proto->get_parameter();
    if(t_density != 0.) {
        for(const auto& I : int_pt) {
            const auto n_int = shape::cube(I.coor, 0);
            const auto tmp_a = t_density * I.jacob_det * I.weight;
            for(unsigned J = 0; J < c_node; ++J)
                for(auto K = J; K < c_node; ++K) initial_mass(c_dof * J, c_dof * K) += tmp_a * n_int(J) * n_int(K);
        }
        for(unsigned I = 0; I < c_size; I += c_dof) {
            initial_mass(I + 1, I + 1) = initial_mass(I + 2, I + 2) = initial_mass(I, I);
            for(auto J = I + c_dof; J < c_size; J += c_dof) initial_mass(J, I) = initial_mass(I + 1, J + 1) = initial_mass(I + 2, J + 2) = initial_mass(J + 1, I + 1) = initial_mass(J + 2, I + 2) = initial_mass(I, J);
        }
    }
    trial_mass = current_mass = initial_mass;
}

int C3D8::update_status() {
    auto code = 0, idx = 0;

    vec t_disp(c_size);
    for(const auto& I : node_ptr) {
        const auto& tmp_disp = I.lock()->get_trial_displacement();
        for(unsigned J = 0; J < c_dof; ++J) t_disp(idx++) = tmp_disp(J);
    }

    trial_stiffness.zeros(c_size, c_size);
    trial_resistance.zeros(c_size);
    for(const auto& I : int_pt) {
        code += I.c_material->update_trial_status(I.strain_mat * t_disp);
        trial_stiffness += I.jacob_det * I.weight * I.strain_mat.t() * I.c_material->get_trial_stiffness() * I.strain_mat;
        trial_resistance += I.jacob_det * I.weight * I.strain_mat.t() * I.c_material->get_trial_stress();
    }

    return code;
}

int C3D8::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->commit_status();
    return code;
}

int C3D8::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->clear_status();
    return code;
}

int C3D8::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->reset_status();
    return code;
}

vector<vec> C3D8::record(const OutputType& T) {
    vector<vec> data;
    switch(T) {
    case OutputType::E:
        for(const auto& I : int_pt) data.emplace_back(I.c_material->get_trial_strain());
        break;
    case OutputType::S:
        for(const auto& I : int_pt) data.emplace_back(I.c_material->get_trial_stress());
        break;
    default:
        for(const auto& I : int_pt)
            for(auto J : I.c_material->record(T)) data.emplace_back(J);
        break;
    }
    return data;
}

void C3D8::print() {
    suanpan_info("C3D8(R) element.\n");
    for(const auto& t_pt : int_pt) {
        suanpan_info("Strain:\t");
        t_pt.c_material->get_trial_strain().t().print();
        suanpan_info("Stress:\t");
        t_pt.c_material->get_trial_stress().t().print();
    }
}
