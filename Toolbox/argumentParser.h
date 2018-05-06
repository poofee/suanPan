/*******************************************************************************
 * Copyright (C) 2017-2018 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
/**
 * @fn argumentParser
 * @brief An argumentParser function.
 *
 * @author tlc
 * @date 21/07/2017
 * @file argumentParser.h
 */

#ifndef ARGUMENTPARSER_H
#define ARGUMENTPARSER_H

#include <Step/Bead.h>

void argument_parser(const int, char**);

void print_header();
void print_version();
void print_helper();

void cli_mode(const shared_ptr<Bead>&);

#endif
