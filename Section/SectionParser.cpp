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
#include <Section/Section>
#include <Toolbox/utility.h>

int create_new_section(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string section_id;
    if(!get_input(command, section_id)) {
        suanpan_info("create_new_section() needs a section type.\n");
        return 0;
    }

    unique_ptr<Section> new_section = nullptr;

    if(is_equal(section_id, "Circle1D"))
        new_circle2d(new_section, command);
    else if(is_equal(section_id, "Rectangle1D"))
        new_circle2d(new_section, command);
    else if(is_equal(section_id, "TrussSection"))
        new_circle2d(new_section, command);
    else if(is_equal(section_id, "Rectangle2D"))
        new_rectangle2d(new_section, command);
    else if(is_equal(section_id, "Circle2D"))
        new_circle2d(new_section, command);
    else if(is_equal(section_id, "ISection2D"))
        new_isection2d(new_section, command);
    else if(is_equal(section_id, "HSection2D"))
        new_hsection2d(new_section, command);
    else if(is_equal(section_id, "Fibre1D"))
        new_fibre1d(new_section, command);
    else if(is_equal(section_id, "Fibre2D"))
        new_fibre2d(new_section, command);
    else if(is_equal(section_id, "Fibre3D"))
        new_fibre3d(new_section, command);
    else {
        // check if the library is already loaded
        auto code = 0;
        for(const auto& I : domain->get_external_module_pool())
            if(I->library_name == section_id) {
                code = 1;
                break;
            }

        // not loaded then try load it
        if(code == 0 && domain->insert(make_shared<ExternalModule>(section_id))) code = 1;

        // if loaded find corresponding function
        if(code == 1)
            for(const auto& I : domain->get_external_module_pool()) {
                if(I->locate_module(section_id)) I->new_object(new_section, command);
                if(new_section != nullptr) break;
            }
    }

    if(new_section == nullptr || !domain->insert(move(new_section))) suanpan_debug("create_new_section() fails to insert new section.\n");

    return 0;
}

void new_rectangle1d(unique_ptr<Section>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_rectangle1d() requires a valid tag.\n");
        return;
    }

    double width;
    if(!get_input(command, width)) {
        suanpan_error("new_rectangle1d() requires a valid width.\n");
        return;
    }

    double height;
    if(!get_input(command, height)) {
        suanpan_error("new_rectangle1d() requires a valid height.\n");
        return;
    }

    unsigned material_id;
    if(!get_input(command, material_id)) {
        suanpan_error("new_rectangle1d() requires a material tag.\n");
        return;
    }

    return_obj = make_unique<Rectangle1D>(tag, width, height, material_id);
}

void new_circle1d(unique_ptr<Section>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_circle1d() requires a valid tag.\n");
        return;
    }

    double radius;
    if(!get_input(command, radius)) {
        suanpan_error("new_circle1d() requires a valid radius.\n");
        return;
    }

    unsigned material_id;
    if(!get_input(command, material_id)) {
        suanpan_error("new_circle1d() requires a material tag.\n");
        return;
    }

    return_obj = make_unique<Circle1D>(tag, radius, material_id);
}

void new_trusssection(unique_ptr<Section>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_trusssection() requires a valid tag.\n");
        return;
    }

    double area;
    if(!get_input(command, area)) {
        suanpan_error("new_trusssection() requires a valid area.\n");
        return;
    }

    unsigned material_id;
    if(!get_input(command, material_id)) {
        suanpan_error("new_trusssection() requires a material tag.\n");
        return;
    }

    return_obj = make_unique<TrussSection>(tag, area, material_id);
}

void new_fibre1d(unique_ptr<Section>&, istringstream&) {}

void new_rectangle2d(unique_ptr<Section>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_rectangle2D() requires a valid tag.\n");
        return;
    }

    double width;
    if(!get_input(command, width)) {
        suanpan_error("new_rectangle2D() requires a valid width.\n");
        return;
    }

    double height;
    if(!get_input(command, height)) {
        suanpan_error("new_rectangle2D() requires a valid height.\n");
        return;
    }

    unsigned material_id;
    if(!get_input(command, material_id)) {
        suanpan_error("new_rectangle2D() requires a material tag.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!get_optional_input(command, int_pt)) suanpan_extra_debug("new_rectangle2D() uses six integration points.\n");

    return_obj = make_unique<Rectangle2D>(tag, width, height, material_id, int_pt);
}

void new_circle2d(unique_ptr<Section>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_circle2D() requires a valid tag.\n");
        return;
    }

    double radius;
    if(!get_input(command, radius)) {
        suanpan_error("new_circle2D() requires a valid radius.\n");
        return;
    }

    unsigned material_id;
    if(!get_input(command, material_id)) {
        suanpan_error("new_circle2D() requires a material tag.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!command.eof() && !get_input(command, int_pt)) {
        suanpan_error("new_circle2D() requires a number of integration points.\n");
        return;
    }

    return_obj = make_unique<Circle2D>(tag, radius, material_id, int_pt);
}

void new_isection2d(unique_ptr<Section>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_isection2d() requires a valid tag.\n");
        return;
    }

    podarray<double> dim(6);
    double size;
    for(auto I = 0; I < 6; ++I) {
        if(!get_input(command, size)) {
            suanpan_error("new_isection2d() requires a valid dimension.\n");
            return;
        }
        dim(I) = size;
    }

    unsigned material_id;
    if(!get_input(command, material_id)) {
        suanpan_error("new_isection2d() requires a material tag.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!command.eof() && !get_input(command, int_pt)) {
        suanpan_error("new_isection2d() requires a number of integration points.\n");
        return;
    }

    return_obj = make_unique<ISection2D>(tag, dim(0), dim(1), dim(2), dim(3), dim(4), dim(5), material_id, int_pt);
}

void new_hsection2d(unique_ptr<Section>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_hsection2d() requires a valid tag.\n");
        return;
    }

    podarray<double> dim(6);
    double size;
    for(auto I = 0; I < 6; ++I) {
        if(!get_input(command, size)) {
            suanpan_error("new_hsection2d() requires a valid dimension.\n");
            return;
        }
        dim(I) = size;
    }

    unsigned material_id;
    if(!get_input(command, material_id)) {
        suanpan_error("new_hsection2d() requires a material tag.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!command.eof() && !get_input(command, int_pt)) {
        suanpan_error("new_hsection2d() requires a number of integration points.\n");
        return;
    }

    return_obj = make_unique<HSection2D>(tag, dim(0), dim(1), dim(2), dim(3), dim(4), dim(5), material_id, int_pt);
}

void new_fibre2d(unique_ptr<Section>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_fibre2d() requires a valid tag.\n");
        return;
    }

    return_obj = make_unique<Fibre2D>(tag);
}

void new_fibre3d(unique_ptr<Section>&, istringstream&) {}
