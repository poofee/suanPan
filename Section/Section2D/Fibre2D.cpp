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

#include "Fibre2D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

Fibre2D::Fibre2D(const unsigned T)
    : Section2D(T, ST_FIBRE2D) {}

Fibre2D::Fibre2D(const Fibre2D& old_obj)
    : Section2D(old_obj.get_tag(), ST_FIBRE2D)
    , fibre_tag(old_obj.fibre_tag) {
    fibre.clear();
    fibre.reserve(old_obj.fibre.size());
    for(const auto& I : old_obj.fibre) fibre.emplace_back(I->get_copy());
}

void Fibre2D::initialize(const shared_ptr<DomainBase>& D) {
    fibre.clear();
    fibre.reserve(fibre_tag.n_elem);
    for(uword I = 0; I < fibre_tag.n_elem; ++I)
        if(D->find_section(fibre_tag(I))) fibre.emplace_back(D->get_section(fibre_tag(I))->get_copy());
}

unique_ptr<Section> Fibre2D::get_copy() { return make_unique<Fibre2D>(*this); }

double Fibre2D::get_parameter(const ParameterType&) { return 0.; }

int Fibre2D::update_trial_status(const vec&) { return 0; }

int Fibre2D::clear_status() { return 0; }

int Fibre2D::commit_status() { return 0; }

int Fibre2D::reset_status() { return 0; }

void Fibre2D::print() {}
