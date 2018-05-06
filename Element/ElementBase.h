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
 * @class ElementBase
 * @brief A ElementBase class.
 * @author tlc
 * @date 21/07/2017
 * @version 0.1.0
 * @file ElementBase.h
 * @addtogroup Element
 * @{
 */

#ifndef ELEMENTBASE_H
#define ELEMENTBASE_H

#include <Domain/Tag.h>

enum class OutputType;

class DomainBase;
class Node;

using std::vector;

class ElementBase : public Tag {
public:
    const bool initialized = false;

    explicit ElementBase(unsigned = 0, // tag
        unsigned = CT_ELEMENTBASE      // class tag
    );
    ElementBase(const ElementBase&) = delete;            // copy forbidden
    ElementBase(ElementBase&&) = delete;                 // move forbidden
    ElementBase& operator=(const ElementBase&) = delete; // assign forbidden
    ElementBase& operator=(ElementBase&&) = delete;      // assign forbidden

    virtual ~ElementBase();

    virtual void initialize(const shared_ptr<DomainBase>&) = 0;

    virtual void update_dof_encoding() = 0;

    virtual const unsigned& get_dof_number() const = 0;
    virtual const unsigned& get_node_number() const = 0;
    virtual const uvec& get_dof_encoding() const = 0;
    virtual const uvec& get_node_encoding() const = 0;

    virtual const vector<weak_ptr<Node>>& get_node_ptr() const = 0;

    virtual const vec& get_resistance() const = 0;

    virtual const mat& get_mass() const = 0;
    virtual const mat& get_damping() const = 0;
    virtual const mat& get_stiffness() const = 0;
    virtual const mat& get_geometry() const = 0;
    virtual const mat& get_secant() const = 0;

    virtual const mat& get_initial_mass() const = 0;
    virtual const mat& get_initial_damping() const = 0;
    virtual const mat& get_initial_stiffness() const = 0;
    virtual const mat& get_initial_geometry() const = 0;
    virtual const mat& get_initial_secant() const = 0;

    virtual Op<vec, op_diagmat> get_diag_mass() const = 0;
    virtual Op<vec, op_diagmat> get_initial_diag_mass() const = 0;

    virtual int update_status() = 0;
    virtual int clear_status() = 0;
    virtual int commit_status() = 0;
    virtual int reset_status() = 0;

    virtual vector<vec> record(const OutputType&) = 0;
};

#endif

//! @}
