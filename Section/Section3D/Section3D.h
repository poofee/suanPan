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
 * @class Section3D
 * @brief A Section3D class.
 * @author tlc
 * @date 27/10/2017
 * @version 0.1.0
 * @file Section3D.h
 * @addtogroup Section-3D
 * @ingroup Section
 * @{
 */

#ifndef SECTION3D_H
#define SECTION3D_H

#include <Section/Section.h>

class Section3D : public Section {
public:
    explicit Section3D(unsigned = 0, unsigned = CT_SECTION, unsigned = 0, double = 0.);
};

#endif

//! @}
