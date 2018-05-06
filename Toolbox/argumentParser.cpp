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

#include "argumentParser.h"
#include "commandParser.h"
#include <Step/Bead.h>
#include <Toolbox/utility.h>
#include <ctime>
#include <suanPan.h>

#ifdef SUANPAN_WIN
#include <Windows.h>
#endif

using std::ofstream;
using std::string;

void argument_parser(const int argc, char** argv) {
    string input_file_name;
    string output_file_name;
    ofstream output_file;
    const auto buffer_backup = cout.rdbuf();

#ifdef SUANPAN_WIN
    CONSOLE_SCREEN_BUFFER_INFO info;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);
    const auto current_attribute = info.wAttributes;
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 10);
#endif

    print_header();

    if(argc > 1) {
        for(auto i = 0; i < argc; ++i) {
            if(is_equal(argv[i], "-v") || is_equal(argv[i], "--version"))
                print_version();
            else if(is_equal(argv[i], "-h") || is_equal(argv[i], "--help"))
                print_helper();
            else if(is_equal(argv[i], "-f") || is_equal(argv[i], "--file"))
                input_file_name = argv[++i];
            else if(is_equal(argv[i], "-o") || is_equal(argv[i], "--output"))
                output_file_name = argv[++i];
        }

        if(!output_file_name.empty()) {
            output_file.open(output_file_name);
            if(output_file.is_open())
                cout.rdbuf(output_file.rdbuf());
            else
                suanpan_error("argumentParser() cannot open the output file.\n");
        }

        if(!input_file_name.empty()) {
            const auto model = make_shared<Bead>();
            if(process_file(model, input_file_name.c_str()) == SUANPAN_EXIT) return;
            cout.rdbuf(buffer_backup);
            cli_mode(model);
        }
    } else {
        const auto model = make_shared<Bead>();
        cli_mode(model);
    }

#ifdef SUANPAN_WIN
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), current_attribute);
#endif
}

void print_header() {
    suanpan_info("+--------------------------------------------------+\n");
    suanpan_info("|   __        __        suanPan is an open source  |\n");
    suanpan_info("|  /  \\      |  \\          FEM framework (%u-bit)  |\n", SUANPAN_ARCH);
    suanpan_info("|  \\__       |__/  __   __          Acrux (0.1.0)  |\n");
    suanpan_info("|     \\ |  | |    |  \\ |  |                        |\n");
    suanpan_info("|  \\__/ |__| |    |__X |  |     maintained by tlc  |\n");
    suanpan_info("|                             all rights reserved  |\n");
    suanpan_info("+--------------------------------------------------+\n\n");
}

void print_version() {
    ostringstream version;
    version << SUANPAN_VERSION;
    auto t_time = time(nullptr);
    suanpan_info("suanPan is an open source FEM framework.\n");
    suanpan_info("\tversion Acrux 0.1.0\n\tcompiled with %s %s\n\tdate %s\n", SUANPAN_COMPILER, version.str().c_str(), asctime(gmtime(&t_time)));
    suanpan_info("[From Wikipedia] Alpha Crucis is a multiple star system located 321 light years from the Sun in the constellation of Crux and part of the asterism known as the Southern Cross.\n\n");
}

void print_helper() {
    suanpan_info("\n");
    suanpan_info("\t-v,  --version\t\tcheck version information\n");
    suanpan_info("\t-h,  --help\t\tprint this helper\n");
    suanpan_info("\t-f,  --file\t\tprocess model file\n");
    suanpan_info("\t-o,  --output\t\tset output file for logging\n\n");
}

void cli_mode(const shared_ptr<Bead>& model) {
    string command_line;
    while(true) {
        suanpan_info("suanPan -> ");
        getline(std::cin, command_line);
        if(!command_line.empty() && command_line[0] != '#') {
            istringstream tmp_str(command_line);
            if(process_command(model, tmp_str) == SUANPAN_EXIT) return;
        }
    }
}
