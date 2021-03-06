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

#ifndef SECTIONPARSER_H
#define SECTIONPARSER_H

#include <suanPan.h>

using std::istringstream;

class Section;
class DomainBase;

int create_new_section(const shared_ptr<DomainBase>&, istringstream&);

// 1D
void new_rectangle1d(unique_ptr<Section>&, istringstream&);
void new_circle1d(unique_ptr<Section>&, istringstream&);
void new_trusssection(unique_ptr<Section>&, istringstream&);
void new_fibre1d(unique_ptr<Section>&, istringstream&);

// 2D
void new_rectangle2d(unique_ptr<Section>&, istringstream&);
void new_circle2d(unique_ptr<Section>&, istringstream&);
void new_isection2d(unique_ptr<Section>&, istringstream&);
void new_hsection2d(unique_ptr<Section>&, istringstream&);
void new_fibre2d(unique_ptr<Section>&, istringstream&);

// 3D
void new_fibre3d(unique_ptr<Section>&, istringstream&);

#endif
