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

#include "T3D2.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material1D/Material1D.h>

const unsigned T3D2::t_node = 2;
const unsigned T3D2::t_dof = 3;
const unsigned T3D2::t_size = t_dof * t_node;

T3D2::T3D2(const unsigned T, const uvec& N, const unsigned M, const double A, const bool F, const bool UA, const bool LS)
    : MaterialElement(T, ET_T3D2, t_node, t_dof, N, uvec{ M }, F)
    , area(A)
    , update_area(UA)
    , log_strain(LS) {}

void T3D2::initialize(const shared_ptr<DomainBase>& D) {
    auto& coord_i = node_ptr.at(0).lock()->get_coordinate();
    auto& coord_j = node_ptr.at(1).lock()->get_coordinate();

    if(coord_i.size() < 3 || coord_j.size() < 3) {
        suanpan_error("initialize() finds incompatible nodes.\n");
        D->disable_element(get_tag());
        return;
    }

    const vec pos_diff = coord_j(span(0, 2)) - coord_i(span(0, 2));

    length = norm(pos_diff);

    direction_cosine = pos_diff / length;

    t_material = D->get_material(unsigned(material_tag(0)))->get_copy();
}

int T3D2::update_status() {
    const auto& node_i = node_ptr.at(0).lock();
    const auto& node_j = node_ptr.at(1).lock();

    const auto& disp_i = node_i->get_trial_displacement();
    const auto& disp_j = node_j->get_trial_displacement();

    // in a truss-beam system a node may have either 2 or 3 dofs depends on the type of elements connected
    // resize the displacement vectors to make sure they are compatiable with the truss formulation
    vec disp_diff = disp_j(span(0, 2)) - disp_i(span(0, 2));

    double trial_strain;

    auto new_area = area;
    auto new_length = length;

    if(nlgeom) {
        disp_diff += node_j->get_coordinate()(span(0, 2)) - node_i->get_coordinate()(span(0, 2));

        new_length = norm(disp_diff);

        direction_cosine = disp_diff / new_length;

        if(update_area) new_area *= length / new_length;

        trial_strain = log_strain ? log(new_length / length) : new_length / length - 1.;
    } else
        trial_strain = dot(disp_diff, direction_cosine) / new_length;

    t_material->update_trial_status(trial_strain);

    const auto t_factor = new_area / new_length * as_scalar(t_material->get_trial_stiffness());

    vec t_vec(t_size);
    t_vec(span(0, 2)) = -direction_cosine;
    t_vec(span(3, 5)) = direction_cosine;

    trial_stiffness = t_factor * t_vec * t_vec.t();

    if(nlgeom) {
        trial_geometry.zeros(t_size, t_size);
        trial_geometry(3, 0) = trial_geometry(4, 1) = trial_geometry(5, 2) = trial_geometry(0, 3) = trial_geometry(1, 4) = trial_geometry(2, 5) = -(trial_geometry(0, 0) = trial_geometry(1, 1) = trial_geometry(2, 2) = trial_geometry(3, 3) = trial_geometry(4, 4) = trial_geometry(5, 5) = new_area / new_length * as_scalar(t_material->get_trial_stress()));
        trial_stiffness += trial_geometry;
    }

    trial_resistance = new_area * as_scalar(t_material->get_trial_stress()) * t_vec;

    return SUANPAN_SUCCESS;
}

int T3D2::commit_status() { return t_material->commit_status(); }

int T3D2::clear_status() { return t_material->clear_status(); }

int T3D2::reset_status() { return t_material->reset_status(); }

void T3D2::print() {
    suanpan_info("3-D truss element with ");
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
