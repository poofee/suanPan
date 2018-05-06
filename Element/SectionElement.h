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
 * @class SectionElement
 * @brief The SectionElement class.
 * @author tlc
 * @date 31/10/2017
 * @version 0.1.0
 * @file SectionElement.h
 * @addtogroup Element
 * @{
 */

#ifndef SECTIONELEMENT_H
#define SECTIONELEMENT_H

#include <Element/Element.h>

class SectionElement : public Element {
public:
    explicit SectionElement(unsigned = 0, // tag
        unsigned = CT_ELEMENT,            // class tag
        unsigned = 0,                     // number of nodes
        unsigned = 0,                     // number of dofs
        const uvec& = {},                 // node encoding
        const uvec& = {},                 // section tags
        bool = false                      // nonlinear geometry switch
    );
    virtual ~SectionElement();
};

#endif

//! @}
