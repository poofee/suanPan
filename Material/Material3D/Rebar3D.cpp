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

#include "Rebar3D.h"
#include <Domain/DomainBase.h>
#include <Toolbox/tensorToolbox.h>

Rebar3D::Rebar3D(const unsigned T, const unsigned XT, const unsigned YT, const unsigned ZT, const double XR, const double YR, const double ZR, const double R)
    : Material3D(T, MT_REBAR3D, R)
    , tag_x(XT)
    , tag_y(YT)
    , tag_z(ZT)
    , ratio_x(XR)
    , ratio_y(YR)
    , ratio_z(ZR) {}

Rebar3D::Rebar3D(const Rebar3D& P)
    : Material3D(P)
    , tag_x(P.tag_x)
    , tag_y(P.tag_y)
    , tag_z(P.tag_z)
    , ratio_x(P.ratio_x)
    , ratio_y(P.ratio_y)
    , ratio_z(P.ratio_z)
    , rebar_x(P.rebar_x->get_copy())
    , rebar_y(P.rebar_y->get_copy())
    , rebar_z(P.rebar_z->get_copy()) {
    Material::initialize();
    Rebar3D::initialize();
}

void Rebar3D::initialize(const shared_ptr<DomainBase>&) {}

unique_ptr<Material> Rebar3D::get_copy() { return make_unique<Rebar3D>(*this); }

int Rebar3D::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;

    rebar_x->update_trial_status(trial_strain(0));
    rebar_y->update_trial_status(trial_strain(1));
    rebar_z->update_trial_status(trial_strain(2));

    trial_stress.zeros(6);
    trial_stress(0) = rebar_x->get_trial_stress().at(0);
    trial_stress(1) = rebar_y->get_trial_stress().at(0);
    trial_stress(2) = rebar_z->get_trial_stress().at(0);

    trial_stiffness.zeros(6, 6);
    trial_stiffness(0, 0) = rebar_x->get_trial_stiffness().at(0);
    trial_stiffness(1, 1) = rebar_y->get_trial_stiffness().at(0);
    trial_stiffness(2, 2) = rebar_z->get_trial_stiffness().at(0);

    return 0;
}

int Rebar3D::clear_status() { return rebar_x->clear_status() + rebar_y->clear_status() + rebar_z->clear_status(); }

int Rebar3D::commit_status() { return rebar_x->commit_status() + rebar_y->commit_status() + rebar_z->commit_status(); }

int Rebar3D::reset_status() { return rebar_x->reset_status() + rebar_y->reset_status() + rebar_z->reset_status(); }
