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

#include "Displacement.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Load/Amplitude/Amplitude.h>

Displacement::Displacement(const unsigned T, const unsigned ST, const double L, const uvec& N, const unsigned D, const unsigned AT)
    : Load(T, CT_DISPLACEMENT, ST == 0 ? 1 : ST, AT, N, uvec{ D }, L) {}

Displacement::Displacement(const unsigned T, const unsigned ST, const double L, const uvec& N, const uvec& D, const unsigned AT)
    : Load(T, CT_DISPLACEMENT, ST == 0 ? 1 : ST, AT, N, D, L) {}

int Displacement::initialize(const shared_ptr<DomainBase>& D) {
    if(initialized) return 0;

    const auto code = Load::initialize(D);

    if(code != 0) return code;

    set_end_step(start_step + 1);

    const auto& t_factory = D->get_factory();

    vector<uword> r_dof;

    for(auto I : t_factory->get_reference_dof()) r_dof.emplace_back(I);

    for(auto I : nodes) {
        auto& t_node = D->get_node(unsigned(I));
        if(t_node != nullptr && t_node->is_active()) {
            auto& t_dof = t_node->get_reordered_dof();
            for(auto J : dofs)
                if(J <= t_dof.n_elem) r_dof.emplace_back(t_dof(J - 1));
        }
    }

    t_factory->set_reference_dof(uvec(r_dof));
    t_factory->set_reference_size(unsigned(r_dof.size()));

    return 0;
}

int Displacement::process(const shared_ptr<DomainBase>& D) {
    const auto& t_step = D->get_current_step_tag();
    if(t_step < start_step || t_step >= end_step) return 0;

    const auto& t_factory = D->get_factory();

    const auto final_settlement = pattern * magnitude->get_amplitude(t_factory->get_trial_time());

    auto& t_settlement = get_trial_settlement(t_factory);

    for(const auto& I : nodes)
        if(D->find_node(unsigned(I))) {
            auto& t_node = D->get_node(unsigned(I));
            if(t_node->is_active()) {
                auto& t_dof = t_node->get_reordered_dof();
                for(const auto& J : dofs)
                    if(J <= t_dof.n_elem) t_settlement(t_dof(J - 1)) = final_settlement;
            }
        }

    t_factory->set_incre_settlement(t_settlement - t_factory->get_current_settlement());

    // const auto final_load = pattern * magnitude->get_amplitude(t_factory->get_trial_time());

    // auto& t_stiff = get_stiffness(t_factory);
    // auto& t_load = get_trial_load(t_factory);

    // for(const auto& I : nodes) {
    //    auto& t_node = D->get_node(static_cast<unsigned>(I));
    //    if(t_node != nullptr && t_node->is_active()) {
    //        auto& t_dof = t_node->get_reordered_dof();
    //        for(const auto& J : dofs)
    //            if(J <= t_dof.n_elem) {
    //                auto& t_idx = t_dof(J - 1);
    //                if(D->insert_constrained_dof(unsigned(t_idx))) {
    //                    if(t_stiff(t_idx, t_idx) == 0) {
    //                        auto& t_set = D->get_constrained_dof();
    //                        t_stiff.at(t_idx, t_idx) = t_set.size() == 1 ? multiplier * t_stiff.max() : *t_set.cbegin() == t_idx ? t_stiff(*++t_set.cbegin(), *++t_set.cbegin()) : t_stiff(*t_set.cbegin(), *t_set.cbegin());
    //                    } else
    //                        t_stiff.at(t_idx, t_idx) *= multiplier;
    //                }
    //                t_load(t_idx) = t_stiff(t_idx, t_idx) * (final_load - t_factory->get_trial_displacement().at(t_idx)); // add unbalanced eqv. load
    //            }
    //    }
    //}

    // t_factory->set_incre_load(t_load - t_factory->get_current_load());

    return 0;
}
