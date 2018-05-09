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

#include <Domain/DomainBase.h>
#include <Domain/ExternalModule.h>
#include <Element/Element>
#include <Toolbox/utility.h>

using std::vector;

int create_new_element(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string element_id;
    if(!get_input(command, element_id)) {
        suanpan_info("create_new_element() needs element type.\n");
        return 0;
    }

    unique_ptr<Element> new_element = nullptr;

    if(is_equal(element_id, "B21"))
        new_b21(new_element, command);
    else if(is_equal(element_id, "B21H"))
        new_b21h(new_element, command);
    else if(is_equal(element_id, "C3D20"))
        new_c3d20(new_element, command);
    else if(is_equal(element_id, "C3D8"))
        new_c3d8(new_element, command);
    else if(is_equal(element_id, "CP3"))
        new_cp3(new_element, command);
    else if(is_equal(element_id, "CP4"))
        new_cp4(new_element, command);
    else if(is_equal(element_id, "CP4R"))
        new_cp4r(new_element, command);
    else if(is_equal(element_id, "CP6"))
        new_cp6(new_element, command);
    else if(is_equal(element_id, "CP8"))
        new_cp8(new_element, command);
    else if(is_equal(element_id, "Damper01"))
        new_damper01(new_element, command);
    else if(is_equal(element_id, "EB21"))
        new_eb21(new_element, command);
    else if(is_equal(element_id, "F21"))
        new_f21(new_element, command);
    else if(is_equal(element_id, "F21H"))
        new_f21h(new_element, command);
    else if(is_equal(element_id, "GCMQ"))
        new_gcmq(new_element, command);
    else if(is_equal(element_id, "GCMT"))
        new_gcmt(new_element, command);
    else if(is_equal(element_id, "GQ12"))
        new_gq12(new_element, command);
    else if(is_equal(element_id, "Mass"))
        new_mass(new_element, command);
    else if(is_equal(element_id, "MVLEM"))
        new_mvlem(new_element, command);
    else if(is_equal(element_id, "PS"))
        new_ps(new_element, command);
    else if(is_equal(element_id, "QE2"))
        new_qe2(new_element, command);
    else if(is_equal(element_id, "RCP4"))
        new_rcp4(new_element, command);
    else if(is_equal(element_id, "RebarLayer"))
        new_rebarlayer(new_element, command);
    else if(is_equal(element_id, "RGCMQ"))
        new_rgcmq(new_element, command);
    else if(is_equal(element_id, "S3"))
        new_s3(new_element, command);
    else if(is_equal(element_id, "S4"))
        new_s4(new_element, command);
    else if(is_equal(element_id, "SingleSection"))
        new_singlesection(new_element, command);
    else if(is_equal(element_id, "T2D2"))
        new_t2d2(new_element, command);
    else if(is_equal(element_id, "T3D2"))
        new_t3d2(new_element, command);
    else {
        // check if the library is already loaded
        auto code = 0;
        for(const auto& I : domain->get_external_module_pool())
            if(I->library_name == element_id) {
                code = 1;
                break;
            }

        // not loaded then try load it
        if(code == 0 && domain->insert(make_shared<ExternalModule>(element_id))) code = 1;

        // if loaded find corresponding function
        if(code == 1)
            for(const auto& I : domain->get_external_module_pool()) {
                if(I->locate_module(element_id)) I->new_object(new_element, command);
                if(new_element != nullptr) break;
            }
    }

    if(new_element == nullptr || !domain->insert(move(new_element))) suanpan_error("create_new_element() fails to create new element.\n");

    return 0;
}

void new_cp3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_cp3() needs a tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 3; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_cp3() needs three valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_debug("new_cp3() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) suanpan_debug("new_cp3() assumes thickness to be unit.\n");

    return_obj = make_unique<CP3>(tag, uvec(node_tag), material_tag, thickness);
}

void new_cp4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_debug("new_cp4() needs a tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 4; ++I) {
        if((command >> node).fail()) {
            suanpan_debug("new_cp4() needs four valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if((command >> material_tag).fail()) {
        suanpan_debug("new_cp4() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof() && !get_input(command, thickness)) {
        suanpan_debug("new_cp4() needs a valid thickness.\n");
        return;
    }

    string reduced_scheme = "N";
    if(!command.eof() && !get_input(command, reduced_scheme)) {
        suanpan_debug("new_cp4() needs a valid reduced scheme switch.\n");
        return;
    }

    string nonlinear = "N";
    if(!command.eof() && !get_input(command, nonlinear)) {
        suanpan_debug("new_cp4() needs a valid nonlinear geometry switch.\n");
        return;
    }

    return_obj = make_unique<CP4>(tag, uvec(node_tag), material_tag, thickness, is_true(reduced_scheme), is_true(nonlinear));
}

void new_cp4r(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_cp4r() needs a tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 4; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_cp4r() needs four valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_debug("new_cp4r() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof()) {
        if(!get_input(command, thickness)) suanpan_debug("new_cp4r() needs a valid thickness.\n");
    } else
        suanpan_debug("new_cp4r() assumes thickness to be unit.\n");

    unsigned nonlinear = 0;
    if(!command.eof()) {
        if(!get_input(command, nonlinear)) suanpan_debug("new_cp4r() needs a valid nonlinear geometry switch (0,1).\n");
    } else
        suanpan_debug("new_cp4r() assumes linear geometry.\n");

    return_obj = make_unique<CP4>(tag, uvec(node_tag), material_tag, thickness, true, !!nonlinear);
}

void new_cp6(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_cp6() needs a tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 6; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_cp6() needs six valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_debug("new_cp6() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof()) {
        if(!get_input(command, thickness)) suanpan_debug("new_cp6() needs a valid thickness.\n");
    } else
        suanpan_extra_debug("new_cp6() assumes thickness to be unit.\n");

    return_obj = make_unique<CP6>(tag, uvec(node_tag), material_tag, thickness);
}

void new_cp8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_debug("new_cp8() needs a tag.\n");
        return;
    }

    uword node;
    vector<uword> node_tag;
    for(auto I = 0; I < 8; ++I) {
        if((command >> node).fail()) {
            suanpan_debug("new_cp8() needs eight valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if((command >> material_tag).fail()) {
        suanpan_debug("new_cp8() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof() && (command >> thickness).fail()) {
        suanpan_debug("new_cp8() needs a valid thickness.\n");
        return;
    }

    unsigned reduced_scheme = 0;
    if(!command.eof() && (command >> reduced_scheme).fail()) {
        suanpan_debug("new_cp8() needs a valid reduced integration switch (0,1).\n");
        return;
    }

    unsigned nonlinear = 0;
    if(!command.eof() && (command >> nonlinear).fail()) {
        suanpan_debug("new_cp8() needs a valid nonlinear geometry switch (0,1).\n");
        return;
    }

    return_obj = make_unique<CP8>(tag, uvec(node_tag), material_tag, thickness, !!reduced_scheme, !!nonlinear);
}

void new_gq12(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_debug("new_gq12() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 4; ++I) {
        if((command >> node).fail()) {
            suanpan_debug("new_gq12() needs four valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if((command >> material_tag).fail()) {
        suanpan_debug("new_gq12() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof()) {
        if((command >> thickness).fail()) suanpan_debug("new_gq12() needs a valid thickness.\n");
    } else
        suanpan_debug("new_gq12() assumes thickness to be unit.\n");

    return_obj = make_unique<GQ12>(tag, uvec(node_tag), material_tag, thickness);
}

void new_ps(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_debug("new_ps() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 4; ++I) {
        if((command >> node).fail()) {
            suanpan_debug("new_ps() needs four valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if((command >> material_tag).fail()) {
        suanpan_debug("new_ps() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof()) {
        if((command >> thickness).fail()) suanpan_debug("new_ps() needs a valid thickness.\n");
    } else
        suanpan_debug("new_ps() assumes thickness to be unit.\n");

    return_obj = make_unique<PS>(tag, uvec(node_tag), material_tag, thickness);
}

void new_qe2(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_debug("new_qe2() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 4; ++I) {
        if((command >> node).fail()) {
            suanpan_debug("new_qe2() needs four valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if((command >> material_tag).fail()) {
        suanpan_debug("new_qe2() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof()) {
        if((command >> thickness).fail()) suanpan_debug("new_qe2() needs a valid thickness.\n");
    } else
        suanpan_debug("new_qe2() assumes thickness to be unit.\n");

    return_obj = make_unique<QE2>(tag, uvec(node_tag), material_tag, thickness);
}

void new_rcp4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_debug("new_rcp4() needs a tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 4; ++I) {
        if((command >> node).fail()) {
            suanpan_debug("new_rcp4() needs four valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if((command >> material_tag).fail()) {
        suanpan_debug("new_rcp4() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof() && !get_input(command, thickness)) {
        suanpan_debug("new_rcp4() needs a valid thickness.\n");
        return;
    }

    auto rho_x = 0., rho_y = 0.;
    auto mat_tag_x = 0, mat_tag_y = 0;
    if(!command.eof() && !get_input(command, rho_x)) suanpan_debug("new_rcp4() needs a valid rho_x.\n");
    if(!command.eof() && !get_input(command, mat_tag_x)) suanpan_debug("new_rcp4() needs a valid mat_x.\n");
    if(!command.eof() && !get_input(command, rho_y)) suanpan_debug("new_rcp4() needs a valid rho_y.\n");
    if(!command.eof() && !get_input(command, mat_tag_y)) suanpan_debug("new_rcp4() needs a valid mat_y.\n");

    string reduced_scheme = "N";
    if(!command.eof() && !get_input(command, reduced_scheme)) {
        suanpan_debug("new_cp4() needs a valid reduced scheme switch.\n");
        return;
    }

    string nonlinear = "N";
    if(!command.eof() && !get_input(command, nonlinear)) {
        suanpan_debug("new_cp4() needs a valid nonlinear geometry switch.\n");
        return;
    }

    return_obj = make_unique<RCP4>(tag, uvec(node_tag), material_tag, thickness, rho_x, mat_tag_x, rho_y, mat_tag_y, is_true(reduced_scheme), is_true(nonlinear));
}

void new_t2d2(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_t2d2() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 2; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_t2d2() needs two valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_debug("new_t2d2() needs a valid material/section tag.\n");
        return;
    }

    auto area = 1.;
    if(!command.eof()) {
        if(!get_input(command, area)) suanpan_debug("new_t2d2() needs a valid area.\n");
    } else
        suanpan_extra_debug("new_t2d2() assumes a unit area.\n");

    string nonlinear = "N", update_area = "N", log_strain = "N";

    if(!command.eof()) {
        if(!get_input(command, nonlinear)) suanpan_debug("new_t2d2() needs a valid nonlinear geometry switch.\n");
    } else
        suanpan_extra_debug("new_t2d2() assumes linear geometry.\n");

    if(!command.eof()) {
        if(!get_input(command, update_area)) suanpan_debug("new_t2d2() needs a valid switch (0,1) to indicate if update area.\n");
    } else
        suanpan_extra_debug("new_truss2d() assumes constant area.\n");

    if(!command.eof()) {
        if(!get_input(command, log_strain)) suanpan_debug("new_t2d2() needs a valid switch (0,1) to indicate if to use engineering strain.\n");
    } else
        suanpan_extra_debug("new_t2d2() assumes engineering strain.\n");

    return_obj = make_unique<T2D2>(tag, uvec(node_tag), material_tag, area, is_true(nonlinear), is_true(update_area), is_true(log_strain));
}

void new_t3d2(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_t3d2() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 2; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_t3d2() needs two valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_debug("new_t3d2() needs a valid material tag.\n");
        return;
    }

    auto area = 1.;
    if(!command.eof()) {
        if(!get_input(command, area)) suanpan_debug("new_t3d2() needs a valid area.\n");
    } else
        suanpan_extra_debug("new_t3d2() assumes a unit area.\n");

    string nonlinear = "N", update_area = "N", log_strain = "N";

    if(!command.eof()) {
        if((command >> nonlinear).fail()) suanpan_debug("new_t3d2() needs a valid nonlinear geometry switch (0,1).\n");
    } else
        suanpan_extra_debug("new_t2d2() assumes linear geometry.\n");

    if(!command.eof()) {
        if((command >> update_area).fail()) suanpan_debug("new_t3d2() needs a valid switch (0,1) to indicate if update area.\n");
    } else
        suanpan_extra_debug("new_truss2d() assumes constant area.\n");

    if(!command.eof()) {
        if((command >> log_strain).fail()) suanpan_debug("new_t3d2() needs a valid switch (0,1) to indicate if to use engineering strain.\n");
    } else
        suanpan_extra_debug("new_t3d2() assumes engineering strain.\n");

    return_obj = make_unique<T3D2>(tag, uvec(node_tag), material_tag, area, is_true(nonlinear), is_true(update_area), is_true(log_strain));
}

void new_c3d8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_debug("new_c3d8() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 8; ++I) {
        if((command >> node).fail()) {
            suanpan_debug("new_c3d8() needs two valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if((command >> material_tag).fail()) {
        suanpan_debug("new_c3d8() needs a valid material tag.\n");
        return;
    }

    string reduced_scheme = "I";
    if(!command.eof()) {
        if((command >> reduced_scheme).fail()) suanpan_debug("new_c3d8() needs a valid reduced integration switch (0,1).\n");
    } else
        suanpan_debug("new_c3d8() assumes standard integration scheme (2*2).\n");

    unsigned nonlinear = 0;
    if(!command.eof()) {
        if((command >> nonlinear).fail()) suanpan_debug("new_c3d8() needs a valid nonlinear geomtery switch (0,1).\n");
    } else
        suanpan_debug("new_c3d8() assumes linear geometry.\n");

    return_obj = make_unique<C3D8>(tag, uvec(node_tag), material_tag, suanpan::to_upper(reduced_scheme[0]), !!nonlinear);
}

void new_c3d20(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_debug("new_c3d20() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 20; ++I) {
        if((command >> node).fail()) {
            suanpan_debug("new_c3d20() needs twenty valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if((command >> material_tag).fail()) {
        suanpan_debug("new_c3d20() needs a valid material tag.\n");
        return;
    }

    unsigned reduced_scheme = 0;
    if(!command.eof()) {
        if((command >> reduced_scheme).fail()) suanpan_debug("new_c3d20() needs a valid reduced integration switch (0,1).\n");
    } else
        suanpan_debug("new_c3d20() assumes standard integration scheme (3*3).\n");

    unsigned nonlinear = 0;
    if(!command.eof()) {
        if((command >> nonlinear).fail()) suanpan_debug("new_c3d20() needs a valid nonlinear geomtery switch (0,1).\n");
    } else
        suanpan_debug("new_c3d20() assumes linear geometry.\n");

    return_obj = make_unique<C3D20>(tag, uvec(node_tag), material_tag, !!reduced_scheme, !!nonlinear);
}

void new_eb21(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_debug("new_eb21() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 2; ++I) {
        if((command >> node).fail()) {
            suanpan_debug("new_eb21() needs two valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    double area;
    if((command >> area).fail()) {
        suanpan_debug("new_eb21() needs a valid area.\n");
        return;
    }

    double moment_inertia;
    if((command >> moment_inertia).fail()) {
        suanpan_debug("new_eb21() needs a valid moment of inertia.\n");
        return;
    }

    unsigned material_tag;
    if((command >> material_tag).fail()) {
        suanpan_debug("new_eb21() needs a valid material tag.\n");
        return;
    }

    unsigned nonlinear = 0;
    if(!command.eof()) {
        if((command >> nonlinear).fail()) suanpan_debug("new_eb21() needs a valid nonlinear geomtery switch (0,1).\n");
    } else
        suanpan_debug("new_eb21() assumes linear geometry.\n");

    return_obj = make_unique<EB21>(tag, uvec(node_tag), area, moment_inertia, material_tag, !!nonlinear);
}

void new_b21(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_b21() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 2; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_b21() needs two valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_debug("new_b21() needs a valid section tag.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!command.eof() && !get_input(command, int_pt)) {
        suanpan_debug("new_b21() needs a valid number of integration points.\n");
        return;
    }

    unsigned nonlinear = 0;
    if(!command.eof()) {
        if(!get_input(command, nonlinear)) suanpan_debug("new_b21() needs a valid nonlinear geomtery switch (0,1).\n");
    } else
        suanpan_debug("new_b21() assumes linear geometry.\n");

    return_obj = make_unique<B21>(tag, uvec(node_tag), section_id, int_pt, !!nonlinear);
}

void new_b21h(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_b21h() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 2; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_b21h() needs two valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_debug("new_b21h() needs a valid section tag.\n");
        return;
    }

    auto elastic_length = .2;
    if(!command.eof() && !get_input(command, elastic_length)) {
        suanpan_debug("new_b21h() needs a valid number of integration points.\n");
        return;
    }

    unsigned nonlinear = 0;
    if(!command.eof()) {
        if(!get_input(command, nonlinear)) suanpan_debug("new_b21h() needs a valid nonlinear geomtery switch (0,1).\n");
    } else
        suanpan_extra_debug("new_b21h() assumes linear geometry.\n");

    return_obj = make_unique<B21H>(tag, uvec(node_tag), section_id, elastic_length, !!nonlinear);
}

void new_f21(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_f21() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 2; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_f21() needs two valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_debug("new_f21() needs a valid section tag.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!command.eof() && !get_input(command, int_pt)) {
        suanpan_debug("new_f21() needs a valid number of integration points.\n");
        return;
    }

    unsigned nonlinear = 0;
    if(!command.eof()) {
        if(!get_input(command, nonlinear)) suanpan_debug("new_f21() needs a valid nonlinear geomtery switch (0,1).\n");
    } else
        suanpan_extra_debug("new_f21() assumes linear geometry.\n");

    return_obj = make_unique<F21>(tag, uvec(node_tag), section_id, int_pt, !!nonlinear);
}

void new_f21h(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_f21h() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 2; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_f21h() needs two valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_debug("new_f21h() needs a valid section tag.\n");
        return;
    }

    auto elastic_length = .2;
    if(!command.eof() && !get_input(command, elastic_length)) {
        suanpan_debug("new_f21h() needs a valid number of integration points.\n");
        return;
    }

    unsigned nonlinear = 0;
    if(!command.eof()) {
        if(!get_input(command, nonlinear)) suanpan_debug("new_f21h() needs a valid nonlinear geomtery switch (0,1).\n");
    } else
        suanpan_extra_debug("new_f21h() assumes linear geometry.\n");

    return_obj = make_unique<F21H>(tag, uvec(node_tag), section_id, elastic_length, !!nonlinear);
}

void new_mass(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_mass() needs a valid tag.\n");
        return;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_debug("new_mass() needs one valid node.\n");
        return;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_debug("new_mass() needs a valid magnitude.\n");
        return;
    }

    unsigned dof;
    vector<uword> dof_tag;
    while(get_input(command, dof)) dof_tag.push_back(dof);

    return_obj = make_unique<Mass>(tag, node, magnitude, uvec(dof_tag));
}

void new_mvlem(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_mvlem() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 2; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_mvlem() needs two valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned shear_tag;
    if(!get_input(command, shear_tag)) {
        suanpan_debug("new_mvlem() needs a valid shear material.\n");
        return;
    }

    double c_height;
    if(!get_input(command, c_height)) {
        suanpan_debug("new_mvlem() needs a valid c.\n");
        return;
    }

    vector<double> B, H, R;
    vector<uword> CT, ST;
    uword t_tag;
    double t_value;
    while(!command.eof()) {
        if(!get_input(command, t_value)) {
            suanpan_debug("new_mvlem() needs a valid fibre width.\n");
            return;
        }
        B.emplace_back(t_value);
        if(!get_input(command, t_value)) {
            suanpan_debug("new_mvlem() needs a valid fibre thickness.\n");
            return;
        }
        H.emplace_back(t_value);
        if(!get_input(command, t_value)) {
            suanpan_debug("new_mvlem() needs a valid fibre reinforcement ratio.\n");
            return;
        }
        R.emplace_back(t_value);
        if(!get_input(command, t_tag)) {
            suanpan_debug("new_mvlem() needs a valid material tag.\n");
            return;
        }
        CT.emplace_back(t_tag);
        if(!get_input(command, t_tag)) {
            suanpan_debug("new_mvlem() needs a valid material tag.\n");
            return;
        }
        ST.emplace_back(t_tag);
    }

    return_obj = make_unique<MVLEM>(tag, uvec(node_tag), B, H, R, uvec(CT), uvec(ST), shear_tag, c_height);
}

void new_damper01(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_damper01() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 2; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_damper01() needs two valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    double damping;
    if(!get_input(command, damping)) {
        suanpan_debug("new_damper01() needs a valid damping coefficient.\n");
        return;
    }

    double alpha;
    if(!get_input(command, alpha)) {
        suanpan_debug("new_damper01() needs a valid alpha.\n");
        return;
    }

    return_obj = make_unique<Damper01>(tag, node_tag, damping, alpha);
}

void new_singlesection(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_singlesection() needs a valid tag.\n");
        return;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_debug("new_singlesection() needs one valid node.\n");
        return;
    }

    unsigned section_tag;
    if(!get_input(command, section_tag)) {
        suanpan_debug("new_singlesection() needs a valid section tag.\n");
        return;
    }

    return_obj = make_unique<SingleSection>(tag, node, section_tag);
}

void new_rebarlayer(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_debug("new_rebarlayer() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 4; ++I) {
        if((command >> node).fail()) {
            suanpan_debug("new_rebarlayer() needs four valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    double config;
    vector<double> configuration;
    for(auto I = 0; I < 9; ++I) {
        if((command >> config).fail()) {
            suanpan_debug("new_rebarlayer() needs nine valid paramters.\n");
            return;
        }
        configuration.push_back(config);
    }

    auto thickness = 1.;
    if(!command.eof()) {
        if((command >> thickness).fail()) suanpan_debug("new_rebarlayer() needs a valid thickness.\n");
    } else
        suanpan_debug("new_rebarlayer() assumes thickness to be unit.\n");

    return_obj = make_unique<RebarLayer>(tag, uvec(node_tag), vec(configuration), thickness);
}

void new_s3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_debug("new_s3() needs a tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 3; ++I) {
        if((command >> node).fail()) {
            suanpan_debug("new_s3() needs four valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if((command >> material_tag).fail()) {
        suanpan_debug("new_s3() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof() && !get_input(command, thickness)) {
        suanpan_debug("new_s3() needs a valid thickness.\n");
        return;
    }

    return_obj = make_unique<S3>(tag, uvec(node_tag), material_tag, thickness);
}

void new_s4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_debug("new_s4() needs a tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 4; ++I) {
        if((command >> node).fail()) {
            suanpan_debug("new_s4() needs four valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if((command >> material_tag).fail()) {
        suanpan_debug("new_s4() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof() && !get_input(command, thickness)) {
        suanpan_debug("new_s4() needs a valid thickness.\n");
        return;
    }

    string reduced_scheme = "G";
    if(!command.eof() && !get_input(command, reduced_scheme)) {
        suanpan_debug("new_s4() needs a valid reduced scheme switch.\n");
        return;
    }

    return_obj = make_unique<S4>(tag, uvec(node_tag), material_tag, thickness, suanpan::to_upper(reduced_scheme[0]));
}

void new_gcmt(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_gcmt() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 3; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_gcmt() needs four valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_debug("new_gcmt() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof()) {
        if((command >> thickness).fail()) suanpan_debug("new_gcmt() needs a valid thickness.\n");
    } else
        suanpan_debug("new_gcmt() assumes thickness to be unit.\n");

    string int_scheme = "I";
    if(!command.eof() && (command >> int_scheme).fail()) suanpan_debug("new_gcmt() needs a valid integration scheme switch.\n");

    return_obj = make_unique<GCMT>(tag, uvec(node_tag), material_tag, thickness, suanpan::to_upper(int_scheme[0]));
}

void new_gcmq(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_gcmq() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 4; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_gcmq() needs four valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_debug("new_gcmq() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof()) {
        if((command >> thickness).fail()) suanpan_debug("new_gcmq() needs a valid thickness.\n");
    } else
        suanpan_debug("new_gcmq() assumes thickness to be unit.\n");

    string int_scheme = "I";
    if(!command.eof() && (command >> int_scheme).fail()) suanpan_debug("new_gcmq() needs a valid reduced scheme switch.\n");

    return_obj = make_unique<GCMQ>(tag, uvec(node_tag), material_tag, thickness, suanpan::to_upper(int_scheme[0]));
}

void new_rgcmq(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("new_rgcmq() needs a valid tag.\n");
        return;
    }

    unsigned node;
    vector<uword> node_tag;
    for(auto I = 0; I < 4; ++I) {
        if(!get_input(command, node)) {
            suanpan_debug("new_rgcmq() needs four valid nodes.\n");
            return;
        }
        node_tag.push_back(node);
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_debug("new_rgcmq() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof()) {
        if((command >> thickness).fail()) suanpan_debug("new_rgcmq() needs a valid thickness.\n");
    } else
        suanpan_debug("new_rgcmq() assumes thickness to be unit.\n");

    string reduced = "I";
    if(!command.eof() && (command >> reduced).fail()) suanpan_debug("new_rgcmq() needs a valid reduced scheme switch.\n");

    auto rho_x = 0., rho_y = 0.;
    auto mat_tag_x = 0, mat_tag_y = 0;
    if(!command.eof() && (command >> rho_x).fail()) suanpan_debug("new_rgcmq() needs a valid rho_x.\n");
    if(!command.eof() && (command >> mat_tag_x).fail()) suanpan_debug("new_rgcmq() needs a valid mat_x.\n");
    if(!command.eof() && (command >> rho_y).fail()) suanpan_debug("new_rgcmq() needs a valid rho_y.\n");
    if(!command.eof() && (command >> mat_tag_y).fail()) suanpan_debug("new_rgcmq() needs a valid mat_y.\n");

    return_obj = make_unique<RGCMQ>(tag, uvec(node_tag), material_tag, thickness, suanpan::to_upper(reduced[0]), rho_x, mat_tag_x, rho_y, mat_tag_y);
}
