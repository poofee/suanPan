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

#include "NodeRecorder.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

NodeRecorder::NodeRecorder(const unsigned T, const uvec& B, const OutputType& L, const bool R, const bool H)
    : Recorder(T, CT_NODERECORDER, B, L, R, H) {}

void NodeRecorder::record(const shared_ptr<DomainBase>& D) {
    auto& obj_tag = get_object_tag();

    for(const auto I : obj_tag)
        if(!D->find_node(unsigned(I))) {
            D->disable_recorder(get_tag());
            return;
        }

    for(auto I = 0; I < int(obj_tag.n_elem); ++I) insert(D->get_node(unsigned(obj_tag(I)))->record(get_variable_type()), I);

    if(if_record_time()) insert(D->get_factory()->get_current_time());
}

void NodeRecorder::print() { suanpan_info("A Node Recorder.\n"); }
