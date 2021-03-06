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

#include "Damper01.h"
#include <Domain/Node.h>
#include <Toolbox/utility.h>

unsigned Damper01::d_node = 2;
unsigned Damper01::d_dof = 2;
unsigned Damper01::d_size = d_dof * d_node;

Damper01::Damper01(const unsigned T, const uvec& NT, const double C, const double A)
    : Element(T, ET_DAMPER01, d_node, d_dof, NT)
    , damper(0, A, C, C, C, C) {}

void Damper01::initialize(const shared_ptr<DomainBase>& D) {
    damper.initialize(D);

    const auto& node_i = node_ptr[0].lock();
    const auto& node_j = node_ptr[1].lock();

    auto& coor_i = node_i->get_coordinate();
    auto& coor_j = node_j->get_coordinate();

    vec pos_diff(2);
    pos_diff(0) = coor_j(0) - coor_i(0);
    pos_diff(1) = coor_j(1) - coor_i(1);

    length = norm(pos_diff);

    if(length < datum::eps) suanpan_warning("zero length detected, check the definition of damper location.\n");
}

int Damper01::update_status() {
    const auto& node_i = node_ptr[0].lock();
    const auto& node_j = node_ptr[1].lock();

    auto& coor_i = node_i->get_coordinate();
    auto& coor_j = node_j->get_coordinate();

    auto& disp_i = node_i->get_trial_displacement();
    auto& disp_j = node_j->get_trial_displacement();

    auto& velocity_i = node_i->get_trial_velocity();
    auto& velocity_j = node_j->get_trial_velocity();

    vec pos_diff(2);
    pos_diff(0) = coor_j(0) - coor_i(0) + disp_j(0) - disp_i(0);
    pos_diff(1) = coor_j(1) - coor_i(1) + disp_j(1) - disp_i(1);

    new_length = norm(pos_diff);

    vec vel_diff(2);
    vel_diff(0) = velocity_j(0) - velocity_i(0);
    vel_diff(1) = velocity_j(1) - velocity_i(1);

    if(new_length == 0.) {
        if(length == 0.) return 0;
        direction_cosine = pos_diff / length;
    } else
        direction_cosine = pos_diff / new_length;

    damper.update_trial_status(vec{ (new_length - length) / length }, vec{ dot(direction_cosine, vel_diff) });

    const auto& t_resistance = damper.get_trial_stress().at(0);

    trial_resistance.zeros(d_size);
    trial_resistance(0) = -(trial_resistance(2) = direction_cosine(0) * t_resistance);
    trial_resistance(1) = -(trial_resistance(3) = direction_cosine(1) * t_resistance);

    return 0;
}

int Damper01::commit_status() { return damper.commit_status(); }

int Damper01::clear_status() { return damper.clear_status(); }

int Damper01::reset_status() { return damper.reset_status(); }
