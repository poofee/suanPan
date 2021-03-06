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
 * @class Element
 * @brief A Element class.
 * @author tlc
 * @date 21/07/2017
 * @version 0.1.0
 * @file Element.h
 * @addtogroup Element
 * @{
 */

#ifndef ELEMENT_H
#define ELEMENT_H

#include <Domain/Tag.h>

class Node;
class DomainBase;
enum class OutputType;
class Material;
class Section;

using std::array;
using std::vector;

class Element : public Tag {
    const unsigned num_node; /**< number of nodes */
    const unsigned num_dof;  /**< number of DoFs */
protected:
    const uvec node_encoding; /**< node encoding */
    const uvec material_tag;  /**< material tags */
    const uvec section_tag;   /**< sectiopn tags */

    const bool nlgeom = false; /**< nonlinear geometry switch */

    uvec dof_encoding; /**< DoF encoding vector */

    vector<weak_ptr<Node>> node_ptr; /**< node pointers */

    mat initial_mass;      /**< mass matrix */
    mat initial_damping;   /**< damping matrix */
    mat initial_stiffness; /**< stiffness matrix */
    mat initial_geometry;  /**< geometry matrix */

    mat trial_mass;      /**< mass matrix */
    mat trial_damping;   /**< damping matrix */
    mat trial_stiffness; /**< stiffness matrix */
    mat trial_geometry;  /**< geometry matrix */

    mat current_mass;      /**< mass matrix */
    mat current_damping;   /**< damping matrix */
    mat current_stiffness; /**< stiffness matrix */
    mat current_geometry;  /**< geometry matrix */

    vec trial_resistance;   /**< resistance vector */
    vec current_resistance; /**< resistance vector */
public:
    const bool initialized = false;

    explicit Element(unsigned = 0, // tag
        unsigned = CT_ELEMENT,     // class tag
        unsigned = 0,              // number of nodes
        unsigned = 0,              // number of dofs
        uvec = {},                 // node encoding
        uvec = {},                 // material tags
        uvec = {},                 // section tags
        bool = false               // nonlinear geometry switch
    );
    Element(const Element&) = delete;            // copy forbidden
    Element(Element&&) = delete;                 // move forbidden
    Element& operator=(const Element&) = delete; // assign forbidden
    Element& operator=(Element&&) = delete;      // assign forbidden

    virtual ~Element();

    virtual void initialize(const shared_ptr<DomainBase>&) = 0;

    void update_dof_encoding();

    const unsigned& get_dof_number() const;
    const unsigned& get_node_number() const;
    const uvec& get_dof_encoding() const;
    const uvec& get_node_encoding() const;

    const vector<weak_ptr<Node>>& get_node_ptr() const;

    virtual const vec& get_resistance() const;

    virtual const mat& get_mass() const;
    virtual const mat& get_damping() const;
    virtual const mat& get_stiffness() const;
    virtual const mat& get_geometry() const;
    virtual const mat& get_secant() const;

    virtual const mat& get_initial_mass() const;
    virtual const mat& get_initial_damping() const;
    virtual const mat& get_initial_stiffness() const;
    virtual const mat& get_initial_geometry() const;
    virtual const mat& get_initial_secant() const;

    virtual Op<vec, op_diagmat> get_diag_mass() const;
    virtual Op<vec, op_diagmat> get_initial_diag_mass() const;

    virtual int update_status() = 0;
    virtual int clear_status() = 0;
    virtual int commit_status() = 0;
    virtual int reset_status() = 0;

    virtual vector<vec> record(const OutputType&);
};

#endif

//! @}
