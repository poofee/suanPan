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

#include "MVLEM.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material.h>

const unsigned MVLEM::b_node = 2;
const unsigned MVLEM::b_dof = 3;
const unsigned MVLEM::b_size = b_node * b_dof;

MVLEM::Fibre::Fibre(const double B, const double H, const double R)
    : width(B)
    , height(H)
    , c_area(B * H * (1. - R))
    , s_area(B * H * R) {}

MVLEM::MVLEM(const unsigned T, const uvec& NT, const vector<double>& B, const vector<double>& H, const vector<double>& R, const uvec& CRT, const uvec& STT, const unsigned SST, const double CH)
    : MaterialElement(T, ET_MVLEM, b_node, b_dof, NT, join_cols(CRT, STT))
    , shear_height(CH)
    , shear_spring_tag(SST) {
    axial_spring.clear(), axial_spring.reserve(B.size());
    auto width_indicator = 0.;
    for(size_t I = 0; I < B.size(); ++I) {
        axial_spring.emplace_back(B[I], H[I], R[I]);
        width_indicator += B[I];
        total_area += B[I] * H[I];
    }
    width_indicator *= -.5;
    for(size_t I = 0; I < B.size(); ++I) {
        axial_spring[I].eccentricity = width_indicator + .5 * axial_spring[I].width;
        width_indicator += axial_spring[I].width;
    }
}

void MVLEM::initialize(const shared_ptr<DomainBase>& D) {
    auto& coord_i = node_ptr.at(0).lock()->get_coordinate();
    auto& coord_j = node_ptr.at(1).lock()->get_coordinate();

    // chord vector
    const vec pos_diff = coord_j(span(0, 1)) - coord_i(span(0, 1));
    length = norm(pos_diff);
    shear_height_a = shear_height * length;
    shear_height_b = shear_height_a - length;

    trans_mat.zeros(6, 6);
    trans_mat(2, 2) = trans_mat(5, 5) = 1.;
    trans_mat(0, 0) = trans_mat(1, 1) = trans_mat(3, 3) = trans_mat(4, 4) = pos_diff(1) / length;
    trans_mat(0, 1) = trans_mat(3, 4) = -(trans_mat(1, 0) = trans_mat(4, 3) = pos_diff(0)) / length;

    const auto& total_fibre_num = axial_spring.size();
    for(size_t I = 0; I < total_fibre_num; ++I) {
        axial_spring[I].c_material = suanpan::make_copy(D->get_material(unsigned(material_tag(I))));
        axial_spring[I].s_material = suanpan::make_copy(D->get_material(unsigned(material_tag(I + total_fibre_num))));
    }

    shear_spring = suanpan::make_copy(D->get_material(shear_spring_tag));

    // form initial stiffness
    auto t_a = 0., t_b = 0., t_c = 0.;
    for(const auto& I : axial_spring) {
        auto t_stiff = I.c_material->get_initial_stiffness().at(0) * I.c_area + I.s_material->get_initial_stiffness().at(0) * I.s_area;
        t_a += t_stiff;
        t_b += t_stiff *= I.eccentricity;
        t_c += t_stiff *= I.eccentricity;
    }

    initial_stiffness.zeros(6, 6);

    t_a /= length, t_b /= length, t_c /= length;
    initial_stiffness(1, 4) = -(initial_stiffness(1, 1) = initial_stiffness(4, 4) = t_a);
    initial_stiffness(2, 5) = -(initial_stiffness(2, 2) = initial_stiffness(5, 5) = t_c);
    initial_stiffness(1, 5) = initial_stiffness(2, 4) = -(initial_stiffness(1, 2) = initial_stiffness(4, 5) = t_b);

    t_c = total_area / length * shear_spring->get_initial_stiffness().at(0);
    t_a = shear_height_a * t_c;
    t_b = shear_height_b * t_c;

    initial_stiffness(0, 2) -= t_a;
    initial_stiffness(2, 3) += t_a;
    initial_stiffness(0, 3) -= t_c;
    initial_stiffness(0, 5) += t_b;
    initial_stiffness(3, 5) -= t_b;
    initial_stiffness(2, 5) -= shear_height_b * t_a;

    initial_stiffness(0, 0) += t_c;
    initial_stiffness(3, 3) += t_c;
    initial_stiffness(2, 2) += shear_height_a * t_a;
    initial_stiffness(5, 5) += shear_height_b * t_b;

    for(auto I = 0; I < 5; ++I)
        for(auto J = I + 1; J < 6; ++J) initial_stiffness(J, I) = initial_stiffness(I, J);

    initial_stiffness = trans_mat.t() * initial_stiffness * trans_mat;

    trial_stiffness = current_stiffness = initial_stiffness;
}

int MVLEM::update_status() {
    const auto& node_i = node_ptr.at(0).lock();
    const auto& node_j = node_ptr.at(1).lock();

    auto& disp_i = node_i->get_trial_displacement();
    auto& disp_j = node_j->get_trial_displacement();

    vec trial_disp(6);
    for(unsigned I = 0; I < b_dof; ++I) {
        trial_disp(I) = disp_i(I);
        trial_disp(I + 3) = disp_j(I);
    }
    // local displacement
    trial_disp = trans_mat * trial_disp;

    vec converter(6, fill::zeros);
    converter(1) = -(converter(4) = 1.);

    auto t_a = 0., t_b = 0., t_c = 0., t_d = 0., t_e = 0.;
    for(const auto& I : axial_spring) {
        converter(2) = -(converter(5) = I.eccentricity);
        const auto trial_strain = dot(converter, trial_disp) / length;
        I.c_material->update_trial_status(trial_strain);
        I.s_material->update_trial_status(trial_strain);
        auto t_stiff = I.c_material->get_trial_stiffness().at(0) * I.c_area + I.s_material->get_trial_stiffness().at(0) * I.s_area;
        t_a += t_stiff;
        t_b += t_stiff *= I.eccentricity;
        t_c += t_stiff *= I.eccentricity;
        const auto t_stress = I.c_material->get_trial_stress().at(0) * I.c_area + I.s_material->get_trial_stress().at(0) * I.s_area;
        t_d += t_stress;
        t_e += t_stress * I.eccentricity;
    }

    trial_stiffness.zeros(6, 6);

    t_a /= length, t_b /= length, t_c /= length;
    trial_stiffness(1, 4) = -(trial_stiffness(1, 1) = trial_stiffness(4, 4) = t_a);
    trial_stiffness(2, 5) = -(trial_stiffness(2, 2) = trial_stiffness(5, 5) = t_c);
    trial_stiffness(1, 5) = trial_stiffness(2, 4) = -(trial_stiffness(1, 2) = trial_stiffness(4, 5) = t_b);

    converter.zeros();
    converter(3) = -(converter(0) = 1.);
    converter(2) = -shear_height_a;
    converter(5) = shear_height_b;

    shear_spring->update_trial_status(dot(converter, trial_disp) / length);

    trial_resistance.zeros(6);
    trial_resistance(3) = -(trial_resistance(0) = shear_spring->get_trial_stress().at(0) * total_area);
    trial_resistance(1) = -(trial_resistance(4) = t_d);
    trial_resistance(2) = -shear_height_a * trial_resistance(0) - t_e;
    trial_resistance(5) = t_e + shear_height_b * trial_resistance(0);

    t_c = total_area / length * shear_spring->get_trial_stiffness().at(0);
    t_d = shear_height_a * t_c;
    t_e = shear_height_b * t_c;

    trial_stiffness(0, 2) -= t_d;
    trial_stiffness(2, 3) += t_d;
    trial_stiffness(0, 3) -= t_c;
    trial_stiffness(0, 5) += t_e;
    trial_stiffness(3, 5) -= t_e;
    trial_stiffness(2, 5) -= shear_height_b * t_d;

    trial_stiffness(0, 0) += t_c;
    trial_stiffness(3, 3) += t_c;
    trial_stiffness(2, 2) += shear_height_a * t_d;
    trial_stiffness(5, 5) += shear_height_b * t_e;

    for(auto I = 0; I < 5; ++I)
        for(auto J = I + 1; J < 6; ++J) trial_stiffness(J, I) = trial_stiffness(I, J);

    // transform back to the gloabl coordinate system
    trial_stiffness = trans_mat.t() * trial_stiffness * trans_mat;
    trial_resistance = trans_mat.t() * trial_resistance;

    return SUANPAN_SUCCESS;
}

int MVLEM::commit_status() {
    auto code = 0;

    code += shear_spring->commit_status();

    for(const auto& I : axial_spring) code += I.c_material->commit_status() + I.s_material->commit_status();

    return code;
}

int MVLEM::clear_status() {
    auto code = 0;

    code += shear_spring->clear_status();

    for(const auto& I : axial_spring) code += I.c_material->clear_status() + I.s_material->clear_status();

    return code;
}

int MVLEM::reset_status() {
    auto code = 0;

    code += shear_spring->reset_status();

    for(const auto& I : axial_spring) code += I.c_material->reset_status() + I.s_material->reset_status();

    return code;
}
