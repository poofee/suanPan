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

#include "F21H.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Recorder/OutputType.h>
#include <Section/Section.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.hpp>
#include <Toolbox/tensorToolbox.h>

const unsigned F21H::b_node = 2;
const unsigned F21H::b_dof = 3;

F21H::IntegrationPoint::IntegrationPoint(const double C, const double W, unique_ptr<Section>&& M)
    : coor(C)
    , weight(W)
    , b_section(move(M))
    , B(2, 3, fill::zeros) {
    B(0, 0) = 1.;
}

mat F21H::quick_inverse(const mat& stiffness) {
    mat flexibility(2, 2, fill::zeros);

    if(stiffness(0, 0) != 0.) flexibility(0, 0) = 1. / stiffness(0, 0);
    if(stiffness(1, 1) != 0.) flexibility(1, 1) = 1. / stiffness(1, 1);

    return flexibility;
}

F21H::F21H(const unsigned T, const uvec& N, const unsigned S, const double L, const bool F)
    : SectionElement(T, ET_F21H, b_node, b_dof, N, uvec{ S }, F)
    , hinge_length(L > .5 ? .5 : L) {}

void F21H::initialize(const shared_ptr<DomainBase>& D) {
    auto& coord_i = node_ptr.at(0).lock()->get_coordinate();
    auto& coord_j = node_ptr.at(1).lock()->get_coordinate();

    // chord vector
    const vec pos_diff = coord_j - coord_i;
    length = norm(pos_diff);
    direction_cosine = pos_diff / length;
    trans_mat = transform::beam::global_to_local(direction_cosine, length);
    inclination = transform::atan2(direction_cosine);

    auto& section_proto = D->get_section(unsigned(section_tag(0)));

    if(section_proto->section_type != SectionType::D2) {
        suanpan_warning("initialize() needs a 2D section.\n");
        D->disable_element(get_tag());
        return;
    }

    // quick computation of section flexibility
    elastic_section_flexibility = quick_inverse(section_proto->get_initial_stiffness());

    // perform integration of elastic region
    const IntegrationPlan plan(1, 2, IntegrationType::GAUSS);
    // add two inner points of Radau quadrature
    const auto int_pt_num = plan.n_rows + 2;
    const auto elastic_length = 1. - 8. * hinge_length;
    // elastic part will be reused in computation
    elastic_local_flexibility.zeros(3, 3);
    // build up the elastic interior
    elastic_int_pt.clear(), elastic_int_pt.reserve(int_pt_num);
    for(unsigned I = 0; I < int_pt_num; ++I) {
        double coor, weight;
        if(I == 0) {
            // left inner Radau point
            coor = 16. / 3. * hinge_length - 1.;
            weight = 3. * hinge_length;
        } else if(I == int_pt_num - 1) {
            // right inner Radau point
            coor = 1. - 16. / 3. * hinge_length;
            weight = 3. * hinge_length;
        } else {
            // Gauss points
            coor = plan(I - 1, 0) * elastic_length;
            weight = .5 * plan(I - 1, 1) * elastic_length;
        }
        elastic_int_pt.emplace_back(coor, weight, section_proto->get_copy());
        // first element moved to ctor
        elastic_int_pt[I].B(1, 1) = (coor - 1.) / 2.;
        elastic_int_pt[I].B(1, 2) = (coor + 1.) / 2.;
        elastic_local_flexibility += elastic_int_pt[I].B.t() * elastic_section_flexibility * elastic_int_pt[I].B * weight * length;
    }

    // perform integration of hinge part
    initial_local_flexibility = elastic_local_flexibility;
    int_pt.clear(), int_pt.reserve(2);
    int_pt.emplace_back(-1., hinge_length, section_proto->get_copy());
    int_pt.emplace_back(1., hinge_length, section_proto->get_copy());
    for(auto& I : int_pt) {
        I.B(1, 1) = .5 * (I.coor - 1.);
        I.B(1, 2) = .5 * (I.coor + 1.);
        initial_local_flexibility += I.B.t() * elastic_section_flexibility * I.B * I.weight * length;
    }

    initial_stiffness = trans_mat.t() * solve(initial_local_flexibility, trans_mat);

    trial_local_flexibility = current_local_flexibility = initial_local_flexibility;

    current_local_deformation.zeros(3);
    trial_local_deformation.zeros(3);
    current_local_resistance.zeros(3);
    trial_local_resistance.zeros(3);
}

int F21H::update_status() {
    auto& disp_i = node_ptr.at(0).lock()->get_trial_displacement();
    auto& disp_j = node_ptr.at(1).lock()->get_trial_displacement();

    vec t_disp(6);
    for(auto I = 0; I < 3; ++I) t_disp(I) = disp_i(I), t_disp(I + 3) = disp_j(I);

    vec residual_deformation = -trial_local_deformation;
    // transform global deformation to local one (remove rigid body motion)
    trial_local_deformation = trans_mat * t_disp;
    // initial residual be aware of how to compute it
    residual_deformation += trial_local_deformation;

    const auto new_length = length;

    auto counter = 0;
    while(true) {
        trial_local_resistance += solve(trial_local_flexibility, residual_deformation);
        residual_deformation.zeros();
        // consider elastic interior no residual generated
        trial_local_flexibility = elastic_local_flexibility;
        // computation of hinge part
        for(const auto& I : int_pt) {
            const vec target_section_resistance = I.B * trial_local_resistance;
            // compute unbalanced deformation use section stiffness
            const vec incre_deformation = (target_section_resistance - I.b_section->get_resistance()) / I.b_section->get_stiffness().diag();
            // update status
            I.b_section->update_trial_status(I.b_section->get_deformation() + incre_deformation);
            // collect new flexibility and deformation
            const mat t_flexibility = I.B.t() * quick_inverse(I.b_section->get_stiffness()) * I.weight * new_length;
            trial_local_flexibility += t_flexibility * I.B;
            residual_deformation += t_flexibility * (I.b_section->get_resistance() - target_section_resistance);
        }
        // computation of elastic part
        // only for deformation tracking no necessary for response
        for(const auto& I : elastic_int_pt) {
            // compute unbalanced deformation use initial section stiffness
            const vec incre_deformation = (I.B * trial_local_resistance - I.b_section->get_resistance()) % elastic_section_flexibility.diag();
            // update status
            I.b_section->update_trial_status(I.b_section->get_deformation() + incre_deformation);
        }
        // quit if converged
        if(norm(residual_deformation) < 1E-12) break;
        // impose a relatively more strict rule
        if(++counter > 5) {
            suanpan_extra_debug("iteration fails to converge at element level.\n");
            return -1;
        }
    }

    trial_stiffness = trans_mat.t() * solve(trial_local_flexibility, trans_mat);
    trial_resistance = trans_mat.t() * trial_local_resistance;

    return 0;
}

int F21H::clear_status() {
    trial_local_flexibility = current_local_flexibility = initial_local_flexibility;
    current_local_deformation.zeros();
    trial_local_deformation.zeros();
    current_local_resistance.zeros();
    trial_local_resistance.zeros();
    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->clear_status();
    return code;
}

int F21H::commit_status() {
    current_local_flexibility = trial_local_flexibility;
    current_local_deformation = trial_local_deformation;
    current_local_resistance = trial_local_resistance;
    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->commit_status();
    return code;
}

int F21H::reset_status() {
    trial_local_flexibility = current_local_flexibility;
    trial_local_deformation = current_local_deformation;
    trial_local_resistance = current_local_resistance;
    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->reset_status();
    return code;
}

vector<vec> F21H::record(const OutputType& P) {
    vector<vec> output;
    output.reserve(int_pt.size() + elastic_int_pt.size());

    if(P == OutputType::E) {
        output.emplace_back(int_pt[0].b_section->get_deformation());
        for(const auto& I : elastic_int_pt) output.emplace_back(I.b_section->get_deformation());
        output.emplace_back(int_pt[1].b_section->get_deformation());
    } else if(P == OutputType::S) {
        output.emplace_back(int_pt[0].b_section->get_resistance());
        for(const auto& I : int_pt) output.emplace_back(I.b_section->get_resistance());
        output.emplace_back(int_pt[1].b_section->get_resistance());
    } else if(P == OutputType::PE) {
        output.emplace_back(int_pt[0].b_section->get_deformation() - int_pt[0].b_section->get_resistance() / int_pt[0].b_section->get_initial_stiffness().diag());
        for(const auto& I : int_pt) output.emplace_back(I.b_section->get_deformation() - I.b_section->get_resistance() / I.b_section->get_initial_stiffness().diag());
        output.emplace_back(int_pt[1].b_section->get_deformation() - int_pt[1].b_section->get_resistance() / int_pt[1].b_section->get_initial_stiffness().diag());
    }

    return output;
}

void F21H::print() { suanpan_info("A Forced-Based Beam Element.\nhttps://doi.org/10.1016/0045-7949(95)00103-N\n"); }
