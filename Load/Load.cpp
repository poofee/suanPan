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

#include "Load.h"
#include <Domain/DomainBase.h>
#include <Load/Amplitude/Ramp.h>
#include <Step/Step.h>
#include <utility>

const double Load::multiplier = 1E8;

Load::Load(const unsigned T, const unsigned CT, const unsigned ST, const unsigned AT, uvec NT, uvec DT, const double PT)
    : Tag(T, CT)
    , start_step(ST)
    , amplitude_tag(AT)
    , nodes(std::move(NT))
    , dofs(std::move(DT))
    , pattern(PT) {}

Load::~Load() { suanpan_debug("Load %u dtor() called.\n", get_tag()); }

int Load::initialize(const shared_ptr<DomainBase>& D) {
    if(magnitude != nullptr) {
        if(!magnitude->is_active()) magnitude = nullptr;
    } else if(amplitude_tag != 0 && D->find_amplitude(amplitude_tag))
        magnitude = D->get_amplitude(amplitude_tag);

    if(amplitude_tag == 0 || magnitude == nullptr) {
        auto t_tag = unsigned(D->get_amplitude()) + 1;
        while(D->find_amplitude(t_tag)) t_tag++;
        magnitude = make_shared<Ramp>(t_tag);
        // D->insert(magnitude);
    }

    auto start_time = 0.;
    for(const auto& I : D->get_step_pool()) {
        if(I.second->get_tag() >= start_step) break;
        start_time += I.second->get_time_period();
    }

    magnitude->set_start_step(start_step);
    magnitude->set_start_time(start_time);

    access::rw(initialized) = true;

    return 0;
}

int Load::process(const shared_ptr<DomainBase>&) { throw; }

void Load::set_start_step(const unsigned T) { start_step = T; }

unsigned Load::get_start_step() const { return start_step; }

void Load::set_end_step(const unsigned T) { end_step = T; }

unsigned Load::get_end_step() const { return end_step; }
