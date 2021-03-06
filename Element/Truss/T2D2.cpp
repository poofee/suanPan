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

#include "T2D2.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material1D/Material1D.h>

const unsigned T2D2::t_node = 2;
const unsigned T2D2::t_dof = 2;
const unsigned T2D2::t_size = t_dof * t_node;

T2D2::T2D2(const unsigned T, const uvec& N, const unsigned M, const double A, const bool F, const bool UA, const bool LS)
    : MaterialElement(T, ET_T2D2, t_node, t_dof, N, uvec{ M }, F)
    , area(A)
    , update_area(UA)
    , log_strain(LS) {}

void T2D2::initialize(const shared_ptr<DomainBase>& D) {
    auto& coord_i = node_ptr.at(0).lock()->get_coordinate();
    auto& coord_j = node_ptr.at(1).lock()->get_coordinate();

    if(coord_i.size() < 2 || coord_j.size() < 2) {
        suanpan_error("initialize() finds incompatible nodes.\n");
        return;
    }

    const vec pos_diff = coord_j(span(0, 1)) - coord_i(span(0, 1));

    length = norm(pos_diff);

    t_material = D->get_material(unsigned(material_tag(0)))->get_copy();

    direction_cosine = pos_diff / length;

    const auto tmp_d = area / length * as_scalar(t_material->get_initial_stiffness());

    initial_stiffness.zeros(t_size, t_size);
    initial_stiffness(0, 2) = initial_stiffness(2, 0) = -(initial_stiffness(0, 0) = initial_stiffness(2, 2) = tmp_d * direction_cosine(0) * direction_cosine(0));
    initial_stiffness(1, 3) = initial_stiffness(3, 1) = -(initial_stiffness(1, 1) = initial_stiffness(3, 3) = tmp_d * direction_cosine(1) * direction_cosine(1));
    initial_stiffness(0, 3) = initial_stiffness(1, 2) = initial_stiffness(2, 1) = initial_stiffness(3, 0) = -(initial_stiffness(0, 1) = initial_stiffness(1, 0) = initial_stiffness(2, 3) = initial_stiffness(3, 2) = tmp_d * direction_cosine(0) * direction_cosine(1));
}

int T2D2::update_status() {
    const auto& node_i = node_ptr.at(0).lock();
    const auto& node_j = node_ptr.at(1).lock();

    // in a truss-beam system a node may have either 2 or 3 dofs depends on the type of elements connected
    // resize the displacement vectors to make sure they are compatiable with the truss formulation
    auto& disp_i = node_i->get_trial_displacement();
    auto& disp_j = node_j->get_trial_displacement();
    vec disp_diff(2);
    disp_diff(0) = disp_j(0) - disp_i(0);
    disp_diff(1) = disp_j(1) - disp_i(1);

    double trial_strain;

    auto new_area = area;
    auto new_length = length;

    if(nlgeom) {
        disp_diff += node_j->get_coordinate() - node_i->get_coordinate();

        new_length = norm(disp_diff);

        direction_cosine = disp_diff / new_length;

        if(update_area) new_area *= length / new_length;

        trial_strain = log_strain ? log(new_length / length) : new_length / length - 1.;
    } else
        trial_strain = dot(disp_diff, direction_cosine) / new_length;

    t_material->update_trial_status(trial_strain);

    const auto t_factor = new_area / new_length * as_scalar(t_material->get_trial_stiffness());

    trial_stiffness.zeros(t_size, t_size);
    trial_stiffness(0, 2) = -(trial_stiffness(0, 0) = trial_stiffness(2, 2) = t_factor * direction_cosine(0) * direction_cosine(0));
    trial_stiffness(1, 3) = -(trial_stiffness(1, 1) = trial_stiffness(3, 3) = t_factor * direction_cosine(1) * direction_cosine(1));
    trial_stiffness(0, 3) = trial_stiffness(1, 2) = -(trial_stiffness(0, 1) = trial_stiffness(2, 3) = t_factor * direction_cosine(0) * direction_cosine(1));

    if(nlgeom) {
        trial_geometry.zeros(t_size, t_size);
        trial_geometry(0, 2) = trial_geometry(1, 3) = -(trial_geometry(0, 0) = trial_geometry(1, 1) = trial_geometry(2, 2) = trial_geometry(3, 3) = new_area / new_length * as_scalar(t_material->get_trial_stress()));
        trial_stiffness += trial_geometry;
    }

    for(auto I = 0; I < 3; ++I)
        for(auto J = I + 1; J < 4; ++J) trial_stiffness(J, I) = trial_stiffness(I, J);

    trial_resistance.zeros(t_size);
    const auto tmp_f = new_area * as_scalar(t_material->get_trial_stress());
    trial_resistance(0) = -(trial_resistance(2) = tmp_f * direction_cosine(0));
    trial_resistance(1) = -(trial_resistance(3) = tmp_f * direction_cosine(1));

    return 0;
}

int T2D2::commit_status() { return t_material->commit_status(); }

int T2D2::clear_status() { return t_material->clear_status(); }

int T2D2::reset_status() { return t_material->reset_status(); }

void T2D2::print() {
    suanpan_info("2-D truss element with ");
    if(nlgeom)
        suanpan_info("corotational formulation, assuming constant %s and %s strain. ", update_area ? "volume" : "area", log_strain ? "logarithmic" : "engineering");
    else
        suanpan_info("linear formulation. ");
    suanpan_info("The nodes connected are\n");
    node_encoding.t().print();
    suanpan_info("The area is %.4E. The initial element length is %.4E.\n", area, length);
    suanpan_info("Material Model: ");
    t_material->print();
}
