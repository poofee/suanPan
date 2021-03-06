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

#include <suanPan>

using std::ifstream;
using std::string;
using std::vector;

int process_command(const shared_ptr<Bead>& model, istringstream& command) {
    if(model == nullptr) return 0;

    const auto command_id = get_input<string>(command);
    if(command.fail()) return 0;

    if(is_equal(command_id, "exit")) return SUANPAN_EXIT;
    if(is_equal(command_id, "quit")) return SUANPAN_EXIT;

    if(is_equal(command_id, "file")) return process_file(model, command);

    if(is_equal(command_id, "domain")) return create_new_domain(model, command);

    if(is_equal(command_id, "enable")) return enable_object(model, command);
    if(is_equal(command_id, "disable")) return disable_object(model, command);
    if(is_equal(command_id, "mute")) return disable_object(model, command);
    if(is_equal(command_id, "erase")) return erase_object(model, command);
    if(is_equal(command_id, "delete")) return erase_object(model, command);
    if(is_equal(command_id, "remove")) return erase_object(model, command);
    if(is_equal(command_id, "save")) return save_object(model, command);

    const auto& domain = get_current_domain(model);

    if(is_equal(command_id, "acceleration")) return create_new_acceleration(domain, command);
    if(is_equal(command_id, "amplitude")) return create_new_amplitude(domain, command);
    if(is_equal(command_id, "cload")) return create_new_cload(domain, command);
    if(is_equal(command_id, "converger")) return create_new_converger(domain, command);
    if(is_equal(command_id, "criterion")) return create_new_criterion(domain, command);
    if(is_equal(command_id, "disp")) return create_new_displacement(domain, command);
    if(is_equal(command_id, "displacement")) return create_new_displacement(domain, command);
    if(is_equal(command_id, "dispload")) return create_new_displacement(domain, command);
    if(is_equal(command_id, "element")) return create_new_element(domain, command);
    if(is_equal(command_id, "fix")) return create_new_bc(domain, command);
    if(is_equal(command_id, "import")) return create_new_external_module(domain, command);
    if(is_equal(command_id, "integrator")) return create_new_integrator(domain, command);
    if(is_equal(command_id, "material")) return create_new_material(domain, command);
    if(is_equal(command_id, "mass")) return create_new_mass(domain, command);
    if(is_equal(command_id, "node")) return create_new_node(domain, command);
    if(is_equal(command_id, "recorder")) return create_new_recorder(domain, command);
    if(is_equal(command_id, "plainrecorder")) return create_new_plainrecorder(domain, command);
    if(is_equal(command_id, "hdf5recorder")) return create_new_hdf5recorder(domain, command);
    if(is_equal(command_id, "section")) return create_new_section(domain, command);
    if(is_equal(command_id, "solver")) return create_new_solver(domain, command);
    if(is_equal(command_id, "step")) return create_new_step(domain, command);

    if(is_equal(command_id, "set")) return set_property(domain, command);

    if(is_equal(command_id, "materialtest")) return test_material(domain, command);

    if(is_equal(command_id, "peek")) return print_info(domain, command);

    if(is_equal(command_id, "analyze")) return model->analyze();

    if(is_equal(command_id, "clear")) {
        domain->clear_status();
        return 0;
    }

    if(is_equal(command_id, "summary")) {
        domain->summary();
        return 0;
    }

    if(is_equal(command_id, "version"))
        print_version();
    else if(is_equal(command_id, "help"))
        print_command_usage(command);

    return 0;
}

int process_file(const shared_ptr<Bead>& model, const char* file_name) {
    ifstream input_file(file_name);

    if(!input_file.is_open()) {
        string new_name = file_name;
        new_name += ".supan";
        input_file.open(new_name);
        if(!input_file.is_open()) {
            suanpan_error("process_file() cannot open the input file.\n");
            return 0;
        }
    }

    string command_line;
    while(!getline(input_file, command_line).fail()) {
        if(!command_line.empty() && command_line[0] != '#') {
            istringstream tmp_str(command_line);
            if(process_command(model, tmp_str) == SUANPAN_EXIT) return SUANPAN_EXIT;
        }
    }
    return 0;
}

int process_file(const shared_ptr<Bead>& model, istringstream& command) {
    string file_name;
    if(!get_input(command, file_name)) {
        suanpan_info("process_file() needs a file name.\n");
        return 0;
    }

    return process_file(model, file_name.c_str());
}

int create_new_domain(const shared_ptr<Bead>& model, istringstream& command) {
    unsigned domain_id;
    if((command >> domain_id).fail()) {
        suanpan_info("create_new_domain() requires a tag.\n");
        return 0;
    }

    model->set_current_domain_tag(domain_id);

    auto& tmp_domain = get_domain(model, domain_id);

    if(tmp_domain == nullptr) {
        tmp_domain = make_shared<Domain>(domain_id);
        if(tmp_domain != nullptr) suanpan_info("create_new_domain() successfully creates Domain %u.\n", domain_id);
    } else
        suanpan_info("create_new_domain() switches to Domain %u.\n", domain_id);

    return 0;
}

int disable_object(const shared_ptr<Bead>& model, istringstream& command) {
    const auto& domain = get_current_domain(model);
    if(domain == nullptr) {
        suanpan_info("disable_object() needs a valid domain.\n");
        return 0;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_info("disable_object() needs object type.\n");
        return 0;
    }

    unsigned tag;
    if(is_equal(object_type, "domain"))
        while(get_input(command, tag)) model->disable_domain(tag);
    else if(is_equal(object_type, "step"))
        while(get_input(command, tag)) domain->disable_step(tag);
    else if(is_equal(object_type, "converger"))
        while(get_input(command, tag)) domain->disable_converger(tag);
    else if(is_equal(object_type, "bc"))
        while(get_input(command, tag)) domain->disable_constraint(tag);
    else if(is_equal(object_type, "constraint"))
        while(get_input(command, tag)) domain->disable_constraint(tag);
    else if(is_equal(object_type, "element"))
        while(get_input(command, tag)) domain->disable_element(tag);
    else if(is_equal(object_type, "load"))
        while(get_input(command, tag)) domain->disable_load(tag);
    else if(is_equal(object_type, "material"))
        while(get_input(command, tag)) domain->disable_material(tag);
    else if(is_equal(object_type, "node"))
        while(get_input(command, tag)) domain->disable_node(tag);
    else if(is_equal(object_type, "recorder"))
        while(get_input(command, tag)) domain->disable_recorder(tag);

    return 0;
}

int enable_object(const shared_ptr<Bead>& model, istringstream& command) {
    const auto& domain = get_current_domain(model);
    if(domain == nullptr) {
        suanpan_info("enable_object() needs a valid domain.\n");
        return 0;
    }

    string object_type;
    if((command >> object_type).fail()) {
        suanpan_info("enable_object() needs object type.\n");
        return 0;
    }

    unsigned tag;
    if(is_equal(object_type, "domain"))
        while(get_input(command, tag)) model->enable_domain(tag);
    else if(is_equal(object_type, "step"))
        while(get_input(command, tag)) domain->enable_step(tag);
    else if(is_equal(object_type, "converger"))
        while(get_input(command, tag)) domain->enable_converger(tag);
    else if(is_equal(object_type, "bc"))
        while(get_input(command, tag)) domain->enable_constraint(tag);
    else if(is_equal(object_type, "constraint"))
        while(get_input(command, tag)) domain->enable_constraint(tag);
    else if(is_equal(object_type, "element"))
        while(get_input(command, tag)) domain->enable_element(tag);
    else if(is_equal(object_type, "load"))
        while(get_input(command, tag)) domain->enable_load(tag);
    else if(is_equal(object_type, "material"))
        while(get_input(command, tag)) domain->enable_material(tag);
    else if(is_equal(object_type, "node"))
        while(get_input(command, tag)) domain->enable_node(tag);
    else if(is_equal(object_type, "recorder"))
        while(get_input(command, tag)) domain->enable_recorder(tag);

    return 0;
}

int erase_object(const shared_ptr<Bead>& model, istringstream& command) {
    const auto& domain = get_current_domain(model);
    if(domain == nullptr) {
        suanpan_info("erase_object() needs a valid domain.\n");
        return 0;
    }

    string object_type;
    if((command >> object_type).fail()) {
        suanpan_info("erase_object() needs object type.\n");
        return 0;
    }

    unsigned tag;
    if(is_equal(object_type, "domain"))
        while(get_input(command, tag)) model->erase_domain(tag);
    else if(is_equal(object_type, "step"))
        while(get_input(command, tag)) domain->erase_step(tag);
    else if(is_equal(object_type, "converger"))
        while(get_input(command, tag)) domain->erase_converger(tag);
    else if(is_equal(object_type, "bc"))
        while(get_input(command, tag)) domain->erase_constraint(tag);
    else if(is_equal(object_type, "constraint"))
        while(get_input(command, tag)) domain->erase_constraint(tag);
    else if(is_equal(object_type, "element"))
        while(get_input(command, tag)) domain->erase_element(tag);
    else if(is_equal(object_type, "load"))
        while(get_input(command, tag)) domain->erase_load(tag);
    else if(is_equal(object_type, "material"))
        while(get_input(command, tag)) domain->erase_material(tag);
    else if(is_equal(object_type, "node"))
        while(get_input(command, tag)) domain->erase_node(tag);
    else if(is_equal(object_type, "recorder"))
        while(get_input(command, tag)) domain->erase_recorder(tag);

    return 0;
}

int save_object(const shared_ptr<Bead>& model, istringstream& command) {
    const auto& domain = get_current_domain(model);
    if(domain == nullptr) {
        suanpan_info("erase_object() needs a valid domain.\n");
        return 0;
    }

    string object_id;
    if(!get_input(command, object_id)) {
        suanpan_info("save_object() needs a valid object type.\n");
        return 0;
    }

    if(is_equal(object_id, "Recorder")) {
        unsigned tag;
        while(get_input(command, tag))
            if(domain->find_recorder(tag)) domain->get_recorder(tag)->save();
    } else if(is_equal(object_id, "Stiffness")) {
        string name;
        if(!command.eof() && !get_input(command, name)) name = "K";
        domain->get_factory()->get_stiffness()->save(name.c_str());
    } else if(is_equal(object_id, "Mass")) {
        string name;
        if(!command.eof() && !get_input(command, name)) name = "M";
        domain->get_factory()->get_mass()->save(name.c_str());
    } else if(is_equal(object_id, "Damping")) {
        string name;
        if(!command.eof() && !get_input(command, name)) name = "C";
        domain->get_factory()->get_damping()->save(name.c_str());
    }

    return 0;
}

int create_new_acceleration(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_info("create_new_acceleration() needs a tag.\n");
        return 0;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_info("create_new_acceleration() needs a valid amplitude tag.\n");
        return 0;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_info("create_new_acceleration() needs load magnitude.\n");
        return 0;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_info("create_new_acceleration() needs a valid DoF.\n");
        return 0;
    }

    const auto& step_tag = domain->get_current_step_tag();

    if(!domain->insert(make_shared<Acceleration>(load_id, step_tag, magnitude, dof_id, amplitude_id))) suanpan_error("create_new_acceleration() fails to create new load.\n");

    return 0;
}

int create_new_amplitude(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string amplitude_type;
    if(!get_input(command, amplitude_type)) {
        suanpan_info("create_new_amplitude() needs a valid amplitude type.\n");
        return 0;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_info("create_new_amplitude() needs a valid amplitude type.\n");
        return 0;
    }

    const auto step_tag = domain->get_current_step_tag();

    if(is_equal(amplitude_type, "Tabular")) {
        string file_name;
        if(!get_input(command, file_name)) {
            suanpan_info("create_new_amplitude() needs a valid file.\n");
            return 0;
        }
        domain->insert(make_shared<Tabular>(tag, file_name.c_str(), step_tag));
    } else if(is_equal(amplitude_type, "Decay")) {
        double A;
        if(!get_input(command, A)) {
            suanpan_info("create_new_amplitude() needs a A.\n");
            return 0;
        }
        double TD;
        if(!get_input(command, TD)) {
            suanpan_info("create_new_amplitude() needs a TD.\n");
            return 0;
        }
        domain->insert(make_shared<Decay>(tag, A, TD, step_tag));
    } else if(is_equal(amplitude_type, "Modulated") || is_equal(amplitude_type, "Sine") || is_equal(amplitude_type, "Cosine")) {
        double A;
        if(!get_input(command, A)) {
            suanpan_info("create_new_amplitude() needs a magnitude.\n");
            return 0;
        }

        double omega;
        vector<double> W;
        while(get_input(command, omega)) W.emplace_back(omega);

        if(is_equal(amplitude_type, "Modulated"))
            domain->insert(make_shared<Modulated>(tag, A, W, step_tag));
        else if(is_equal(amplitude_type, "Sine"))
            domain->insert(make_shared<Sine>(tag, A, W, step_tag));
        else if(is_equal(amplitude_type, "Cosine"))
            domain->insert(make_shared<Cosine>(tag, A, W, step_tag));
    }

    return 0;
}

int create_new_bc(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned bc_id;
    if(!get_input(command, bc_id)) {
        suanpan_info("create_new_bc() needs BC type.\n");
        return 0;
    }

    string dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_info("create_new_bc() needs valid DoFs.\n");
        return 0;
    }

    unsigned node;
    vector<uword> node_tag;
    while(get_input(command, node)) node_tag.push_back(node);

    const auto& step_tag = domain->get_current_step_tag();

    const auto bc_type = suanpan::to_lower(dof_id[0]);
    if(is_equal(bc_type, 'p'))
        domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), "PINNED"));
    else if(is_equal(bc_type, 'e'))
        domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), "ENCASTRE"));
    else if(is_equal(bc_type, 'x'))
        domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), "XSYMM"));
    else if(is_equal(bc_type, 'y'))
        domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), "YSYMM"));
    else if(is_equal(bc_type, 'z'))
        domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), "ZSYMM"));
    else if(is_equal(bc_type, '1'))
        domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), 1));
    else if(is_equal(bc_type, '2'))
        domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), 2));
    else if(is_equal(bc_type, '3'))
        domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), 3));
    else if(is_equal(bc_type, '4'))
        domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), 4));
    else if(is_equal(bc_type, '5'))
        domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), 5));
    else if(is_equal(bc_type, '6'))
        domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), 6));

    return 0;
}

int create_new_cload(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_info("create_new_cload() needs a tag.\n");
        return 0;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_info("create_new_cload() needs a valid amplitude tag.\n");
        return 0;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_info("create_new_cload() needs load magnitude.\n");
        return 0;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_info("create_new_cload() needs a valid DoF.\n");
        return 0;
    }

    unsigned node;
    vector<uword> node_tag;
    while(get_input(command, node)) node_tag.push_back(node);

    if(!domain->insert(make_shared<CLoad>(load_id, domain->get_current_step_tag(), magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_cload() fails to create new load.\n");

    return 0;
}

int create_new_converger(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string converger_id;
    if(!get_input(command, converger_id)) {
        suanpan_info("create_new_converger() requires converger type.\n");
        return 0;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_info("create_new_converger() requires a tag.\n");
        return 0;
    }

    auto tolerance = 1E-6;
    if(!command.eof() && !get_input(command, tolerance)) {
        suanpan_info("create_new_converger() reads wrong tolerance.\n");
        return 0;
    }

    auto max_iteration = 10;
    if(!command.eof() && !get_input(command, max_iteration)) {
        suanpan_info("create_new_converger() reads wrong max iteration.\n");
        return 0;
    }

    string print_flag = "false";
    if(!command.eof() && !get_input(command, print_flag)) {
        suanpan_info("create_new_converger() reads wrong print flag.\n");
        return 0;
    }

    auto code = 0;
    if(is_equal(converger_id, "AbsResidual")) {
        if(domain->insert(make_shared<AbsResidual>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
    } else if(is_equal(converger_id, "RelResidual")) {
        if(domain->insert(make_shared<RelResidual>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
    } else if(is_equal(converger_id, "AbsIncreDisp")) {
        if(domain->insert(make_shared<AbsIncreDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
    } else if(is_equal(converger_id, "RelIncreDisp")) {
        if(domain->insert(make_shared<RelIncreDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
    } else if(is_equal(converger_id, "AbsIncreEnergy")) {
        if(domain->insert(make_shared<AbsIncreEnergy>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
    } else if(is_equal(converger_id, "RelIncreEnergy")) {
        if(domain->insert(make_shared<RelIncreEnergy>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
    } else
        suanpan_info("create_new_converger() cannot identify the converger type.\n");

    if(code == 1) {
        if(domain->get_current_step_tag() != 0) domain->get_current_step()->set_converger_tag(tag);
        domain->set_current_converger_tag(tag);
    } else
        suanpan_info("create_new_converger() fails to create the new converger.\n");

    return 0;
}

int create_new_criterion(const shared_ptr<DomainBase>& domain, istringstream& command) {
    const auto& step_tag = domain->get_current_step_tag();
    if(step_tag == 0) {
        suanpan_info("create_new_criterion() needs a valid step.\n");
        return 0;
    }

    string criterion_type;
    if(!get_input(command, criterion_type)) {
        suanpan_info("create_new_criterion() need a criterion type.\n");
        return 0;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_info("create_new_criterion() requires a tag.\n");
        return 0;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_info("create_new_criterion() requires a node.\n");
        return 0;
    }

    unsigned dof;
    if(!get_input(command, dof)) {
        suanpan_info("create_new_criterion() requires a dof.\n");
        return 0;
    }

    double limit;
    if(!get_input(command, limit)) {
        suanpan_info("create_new_criterion() requires a limit.\n");
        return 0;
    }

    if(is_equal(criterion_type, "MaxDisplacement"))
        domain->insert(make_shared<MaxDisplacement>(tag, step_tag, node, dof, limit));
    else if(is_equal(criterion_type, "MinDisplacement"))
        domain->insert(make_shared<MinDisplacement>(tag, step_tag, node, dof, limit));

    return 0;
}

int create_new_displacement(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_info("create_new_displacement() needs a tag.\n");
        return 0;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_info("create_new_displacement() needs a valid amplitude tag.\n");
        return 0;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_info("create_new_displacement() needs load magnitude.\n");
        return 0;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_info("create_new_displacement() needs a valid DoF.\n");
        return 0;
    }

    unsigned node;
    vector<uword> node_tag;
    while(get_input(command, node)) node_tag.push_back(node);

    const auto& step_tag = domain->get_current_step_tag();

    if(!domain->insert(make_shared<Displacement>(load_id, step_tag, magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_displacement() fails to create new load.\n");

    return 0;
}

int create_new_external_module(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string library_name;

    if(!get_input(command, library_name)) {
        suanpan_info("create_new_external_module() needs module name.\n");
        return 0;
    }

    auto code = 0;
    for(const auto& I : domain->get_external_module_pool())
        if(I->library_name == library_name) {
            code = 1;
            break;
        }

    if(code == 0) domain->insert(make_shared<ExternalModule>(library_name));

    return 0;
}

int create_new_integrator(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string integrator_type;
    if(!get_input(command, integrator_type)) {
        suanpan_error("create_new_integrator() needs a valid integrator type.\n");
        return 0;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_integrator() needs a valid tag.\n");
        return 0;
    }

    auto code = 0;
    if(is_equal(integrator_type, "Newmark")) {
        auto alpha = .25, beta = .5;
        if(!command.eof()) {
            if(!get_input(command, alpha)) {
                suanpan_error("create_new_integrator() needs a valid alpha.\n");
                return 0;
            }
            if(!get_input(command, beta)) {
                suanpan_error("create_new_integrator() needs a valid beta.\n");
                return 0;
            }
            if(domain->insert(make_shared<Newmark>(tag, alpha, beta))) code = 1;
        }
    } else if(is_equal(integrator_type, "GeneralizedAlpha")) {
        auto alpha_m = .0, alpha_f = .0;
        if(!command.eof()) {
            if(!get_input(command, alpha_m)) {
                suanpan_error("create_new_integrator() needs a valid alpha_m.\n");
                return 0;
            }
            if(!get_input(command, alpha_f)) {
                suanpan_error("create_new_integrator() needs a valid alpha_f.\n");
                return 0;
            }
            if(domain->insert(make_shared<GeneralizedAlpha>(tag, alpha_m, alpha_f))) code = 1;
        }
    } else if(is_equal(integrator_type, "CentralDifference"))
        if(domain->insert(make_shared<CentralDifference>(tag))) code = 1;

    if(code == 1) {
        if(domain->get_current_step_tag() != 0) domain->get_current_step()->set_integrator_tag(tag);
        domain->set_current_integrator_tag(tag);
    } else
        suanpan_info("create_new_integrator() fails to create the new integrator.\n");

    return 0;
}

int create_new_mass(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_debug("create_new_mass() needs a valid tag.\n");
        return 0;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_debug("create_new_mass() needs one valid node.\n");
        return 0;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_debug("create_new_mass() needs a valid magnitude.\n");
        return 0;
    }

    unsigned dof;
    vector<uword> dof_tag;
    while(get_input(command, dof)) dof_tag.push_back(dof);

    domain->insert(make_shared<Mass>(tag, node, magnitude, uvec(dof_tag)));

    return 0;
}

int create_new_node(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned node_id;
    if(!get_input(command, node_id)) {
        suanpan_info("create_new_node() needs a tag.\n");
        return 0;
    }

    vector<double> coor;
    double X;
    while(get_input(command, X)) coor.push_back(X);

    if(!domain->insert(make_shared<Node>(node_id, vec(coor)))) suanpan_debug("create_new_node() fails to insert Node %u.\n", node_id);

    return 0;
}

int create_new_recorder(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_info("create_new_recorder() needs a valid tag.\n");
        return 0;
    }

    string file_type;
    if(!get_input(command, file_type)) {
        suanpan_info("create_new_recorder() needs a valid object type.\n");
        return 0;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_info("create_new_recorder() needs a valid object type.\n");
        return 0;
    }

    string variable_type;
    if(!get_input(command, variable_type)) {
        suanpan_info("create_new_recorder() needs a valid recorder type.\n");
        return 0;
    }

    unsigned s_object_tag;
    vector<uword> object_tag;
    while(!command.eof() && get_input(command, s_object_tag)) object_tag.emplace_back(s_object_tag);

    if(is_equal(object_type, "Node") && !domain->insert(make_shared<NodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), true, is_equal(file_type[0], 'h'))))
        suanpan_info("create_new_recorder() fails to create a new node recorder.\n");
    else if(is_equal(object_type, "Element") && !domain->insert(make_shared<ElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), true, is_equal(file_type[0], 'h'))))
        suanpan_info("create_new_recorder() fails to create a new element recorder.\n");

    return 0;
}

int create_new_plainrecorder(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_info("create_new_plainrecorder() needs a valid tag.\n");
        return 0;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_info("create_new_plainrecorder() needs a valid object type.\n");
        return 0;
    }

    string variable_type;
    if(!get_input(command, variable_type)) {
        suanpan_info("create_new_plainrecorder() needs a valid recorder type.\n");
        return 0;
    }

    unsigned s_object_tag;
    vector<uword> object_tag;
    while(!command.eof() && get_input(command, s_object_tag)) object_tag.emplace_back(s_object_tag);

    if(is_equal(object_type, "Node") && !domain->insert(make_shared<NodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), true, false)))
        suanpan_info("create_new_plainrecorder() fails to create a new node recorder.\n");
    else if(is_equal(object_type, "Element") && !domain->insert(make_shared<ElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), true, false)))
        suanpan_info("create_new_plainrecorder() fails to create a new element recorder.\n");

    return 0;
}

int create_new_hdf5recorder(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_info("create_new_hdf5recorder() needs a valid tag.\n");
        return 0;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_info("create_new_hdf5recorder() needs a valid object type.\n");
        return 0;
    }

    string variable_type;
    if(!get_input(command, variable_type)) {
        suanpan_info("create_new_hdf5recorder() needs a valid recorder type.\n");
        return 0;
    }

    unsigned s_object_tag;
    vector<uword> object_tag;
    while(!command.eof() && get_input(command, s_object_tag)) object_tag.emplace_back(s_object_tag);

    if(is_equal(object_type, "Node") && !domain->insert(make_shared<NodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), true, true)))
        suanpan_info("create_new_hdf5recorder() fails to create a new node recorder.\n");
    else if(is_equal(object_type, "Element") && !domain->insert(make_shared<ElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), true, true)))
        suanpan_info("create_new_hdf5recorder() fails to create a new element recorder.\n");

    return 0;
}

int create_new_solver(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string solver_type;
    if((command >> solver_type).fail()) {
        suanpan_info("create_new_solver() requires solver type.\n");
        return 0;
    }

    unsigned tag;
    if((command >> tag).fail()) {
        suanpan_info("create_new_solver() requires a tag.\n");
        return 0;
    }

    auto code = 0;
    if(is_equal(solver_type, "Newton")) {
        if(domain->insert(make_shared<Newton>(tag))) code = 1;
    } else if(is_equal(solver_type, "modifiedNewton")) {
        if(domain->insert(make_shared<Newton>(tag, true))) code = 1;
    } else if(is_equal(solver_type, "mNewton")) {
        if(domain->insert(make_shared<Newton>(tag, true))) code = 1;
    } else if(is_equal(solver_type, "BFGS")) {
        if(domain->insert(make_shared<BFGS>(tag))) code = 1;
    } else if(is_equal(solver_type, "Ramm")) {
        if(domain->insert(make_shared<Ramm>(tag))) code = 1;
    } else if(is_equal(solver_type, "DisplacementControl")) {
        if(domain->insert(make_shared<MPDC>(tag))) code = 1;
    } else if(is_equal(solver_type, "MPDC")) {
        if(domain->insert(make_shared<MPDC>(tag))) code = 1;
    } else
        suanpan_error("create_new_solver() cannot identify solver type.\n");

    if(code == 1) {
        if(domain->get_current_step_tag() != 0) domain->get_current_step()->set_solver_tag(tag);
        domain->set_current_solver_tag(tag);
    } else
        suanpan_error("create_new_solver() cannot create the new solver.\n");

    return 0;
}

int create_new_step(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string step_type;
    if(!get_input(command, step_type)) {
        suanpan_info("create_new_step() requires step type.\n");
        return 0;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_info("create_new_step() requires a tag.\n");
        return 0;
    }

    if(is_equal(step_type, "Static")) {
        auto time = 1.;
        if(!command.eof() && !get_input(command, time)) {
            suanpan_info("create_new_step() reads a wrong time period.\n");
            return 0;
        }
        if(domain->insert(make_shared<Static>(tag, time)))
            domain->set_current_step_tag(tag);
        else
            suanpan_error("create_new_step() cannot create the new step.\n");
    } else if(is_equal(step_type, "Dynamic")) {
        auto time = 1.;
        if(!command.eof() && !get_input(command, time)) {
            suanpan_info("create_new_step() reads a wrong time period.\n");
            return 0;
        }
        if(domain->insert(make_shared<Dynamic>(tag, time)))
            domain->set_current_step_tag(tag);
        else
            suanpan_error("create_new_step() cannot create the new step.\n");
    } else if(is_equal(step_type, "ArcLength")) {
        unsigned node;
        if(!get_input(command, node)) {
            suanpan_info("create_new_step() requires a node.\n");
            return 0;
        }

        unsigned dof;
        if(!get_input(command, dof)) {
            suanpan_info("create_new_step() requires a dof.\n");
            return 0;
        }

        double magnitude;
        if(!get_input(command, magnitude)) {
            suanpan_info("create_new_step() requires a magnitude.\n");
            return 0;
        }

        if(domain->insert(make_shared<ArcLength>(tag, node, dof, magnitude)))
            domain->set_current_step_tag(tag);
        else
            suanpan_error("create_new_step() cannot create the new step.\n");
    } else
        suanpan_info("create_new_step() cannot identify step type.\n");

    return 0;
}

int set_property(const shared_ptr<DomainBase>& domain, istringstream& command) {
    if(domain->get_current_step_tag() == 0) return 0;

    const auto& tmp_step = domain->get_current_step();

    string property_id;
    if(!get_input(command, property_id)) {
        suanpan_info("set_property() need a property type.\n");
        return 0;
    }

    if(is_equal(property_id, "fixed_step_size")) {
        string value;
        get_input(command, value) ? tmp_step->set_fixed_step_size(is_true(value)) : suanpan_info("set_property() need a valid value.\n");
    } else if(is_equal(property_id, "symm_mat")) {
        string value;
        get_input(command, value) ? tmp_step->set_symm(is_true(value)) : suanpan_info("set_property() need a valid value.\n");
    } else if(is_equal(property_id, "band_mat")) {
        string value;
        get_input(command, value) ? tmp_step->set_band(is_true(value)) : suanpan_info("set_property() need a valid value.\n");
    } else if(is_equal(property_id, "sparse_mat")) {
        string value;
        get_input(command, value) ? tmp_step->set_sparse(is_true(value)) : suanpan_info("set_property() need a valid value.\n");
    } else if(is_equal(property_id, "ini_step_size")) {
        double step_time;
        get_input(command, step_time) ? tmp_step->set_ini_step_size(step_time) : suanpan_info("set_property() need a valid value.\n");
    } else if(is_equal(property_id, "min_step_size")) {
        double step_time;
        get_input(command, step_time) ? tmp_step->set_min_step_size(step_time) : suanpan_info("set_property() need a valid value.\n");
    } else if(is_equal(property_id, "max_step_size")) {
        double step_time;
        get_input(command, step_time) ? tmp_step->set_max_step_size(step_time) : suanpan_info("set_property() need a valid value.\n");
    } else if(is_equal(property_id, "max_iteration")) {
        unsigned max_number;
        get_input(command, max_number) ? tmp_step->set_max_substep(max_number) : suanpan_info("set_property() need a valid value.\n");
    }

    return 0;
}

int print_info(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string object_type;
    if((command >> object_type).fail()) {
        suanpan_info("print_info() needs object type.\n");
        return 0;
    }

    unsigned tag;
    if(is_equal(object_type, "node"))
        while(get_input(command, tag)) {
            const auto& tmp_node = get_node(domain, tag);
            if(tmp_node != nullptr) tmp_node->print();
            suanpan_info("\n");
        }
    else if(is_equal(object_type, "element"))
        while(get_input(command, tag)) {
            const auto& tmp_element = get_element(domain, tag);
            if(tmp_element != nullptr) tmp_element->print();
            suanpan_info("\n");
        }
    else if(is_equal(object_type, "material"))
        while(get_input(command, tag)) {
            const auto& tmp_material = get_material(domain, tag);
            if(tmp_material != nullptr) tmp_material->print();
            suanpan_info("\n");
        }
    else if(is_equal(object_type, "constraint"))
        while(get_input(command, tag)) {
            const auto& tmp_constraint = get_constraint(domain, tag);
            if(tmp_constraint != nullptr) tmp_constraint->print();
            suanpan_info("\n");
        }
    else if(is_equal(object_type, "recorder"))
        while(get_input(command, tag)) {
            const auto& tmp_recorder = get_recorder(domain, tag);
            if(tmp_recorder != nullptr) tmp_recorder->print();
            suanpan_info("\n");
        }

    return 0;
}

void print_command_usage(istringstream& command) {
    string command_id;
    command >> command_id;

    if(is_equal(command_id, "converger")) {
        suanpan_info("\nconverger $type $tag $tolerance [$max_iteration] [$if_print]\n");
        suanpan_info("\t$type --- converger type\n");
        suanpan_info("\t$tag --- converger tag\n");
        suanpan_info("\t$tolerance --- tolerance -> 1E-8\n");
        suanpan_info("\t$max_iteration --- maximum iteration number -> 7\n");
        suanpan_info("\t$if_print --- print error in each iteration -> false\n\n");
    } else if(is_equal(command_id, "step")) {
        suanpan_info("\nstep $type $tag [$time_period]\n");
        suanpan_info("\t$type --- step type\n");
        suanpan_info("\t$tag --- step tag\n");
        suanpan_info("\t$time_period --- step time period -> 1.0\n\n");
    } else if(is_equal(command_id, "Truss2D")) {
        suanpan_info("\nelement Truss2D $tag {$node_tag...} $material_tag $area [$nonlinear_switch] [$constant_area_switch] [$log_strain_switch]\n");
        suanpan_info("\t$tag --- element tag\n");
        suanpan_info("\t$node_tag --- node tag (2)\n");
        suanpan_info("\t$material_tag --- material tag\n");
        suanpan_info("\t$area --- cross section area\n");
        suanpan_info("\t$nonlinear_switch --- if to use corotational formulation -> false\n");
        suanpan_info("\t$constant_area_switch --- if to update area based on constant volume assumption -> false\n");
        suanpan_info("\t$log_strain_switch --- if to use log strain or engineering strain -> false\n\n");
    } else if(is_equal(command_id, "Elastic1D")) {
        suanpan_info("\nmaterial Elastic1D $tag $elastic_modulus [$density]\n");
        suanpan_info("\t$tag --- material tag\n");
        suanpan_info("\t$elastic_modulus --- elastic modulus\n");
        suanpan_info("\t$density --- density -> 0.0\n\n");
    } else if(is_equal(command_id, "Bilinear1D")) {
        suanpan_info("\nmaterial Bilinear1D $tag $elastic_modulus $yield_stress [$hardening_ratio] [$beta] [$density]\n");
        suanpan_info("\t$tag --- material tag\n");
        suanpan_info("\t$elastic_modulus --- elastic modulus\n");
        suanpan_info("\t$yield_stress --- yield stress\n");
        suanpan_info("\t$hardening_ratio --- hardening ratio -> 0.0\n");
        suanpan_info("\t$beta --- mixed hardening 0.0 for isotropic hardening 1.0 for kinematic hardening -> 0.0\n");
        suanpan_info("\t$density --- density -> 0.0\n\n");
    }
}
