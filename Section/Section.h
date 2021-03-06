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
 * @class Section
 * @brief A Section class.
 * @author tlc
 * @date 27/07/2017
 * @version 0.1.0
 * @file Section.h
 * @addtogroup Section
 * @{
 */

#ifndef SECTION_H
#define SECTION_H

#include <Domain/Tag.h>
#include <Section/ParameterType.h>

enum class SectionType : unsigned { D0 = 0, D1 = 1, D2 = 2, D3 = 3 };

class DomainBase;
class Material;

using std::vector;

class Section : public Tag {
protected:
    const unsigned material_tag; /**< material tags */

    const double area;

    vec eccentricity;

    vec trial_deformation;   /**< trial deformation */
    vec current_deformation; /**< current deformation */

    vec trial_deformation_rate;   /**< trial deformation rate */
    vec current_deformation_rate; /**< current deformation rate */

    vec trial_resistance;   /**< trial resistance */
    vec current_resistance; /**< current resistance */

    mat initial_stiffness; /**< initial stiffness matrix */
    mat current_stiffness; /**< stiffness matrix */
    mat trial_stiffness;   /**< stiffness matrix */
public:
    const SectionType section_type;

    explicit Section(unsigned T = 0, unsigned CT = CT_SECTION, SectionType = SectionType::D0, unsigned MT = 0, double = 0.);

    virtual ~Section();

    virtual void initialize(const shared_ptr<DomainBase>&) = 0;

    virtual const vec& get_deformation() const;
    virtual const vec& get_deformation_rate() const;
    virtual const vec& get_resistance() const;
    virtual const mat& get_stiffness() const;
    virtual const mat& get_initial_stiffness() const;

    virtual unique_ptr<Section> get_copy() = 0;

    virtual double get_parameter(const ParameterType& = ParameterType::NONE);

    virtual int update_incre_status(const vec&);
    virtual int update_incre_status(const vec&, const vec&);
    virtual int update_trial_status(const vec&);
    virtual int update_trial_status(const vec&, const vec&);

    virtual int clear_status() = 0;
    virtual int commit_status() = 0;
    virtual int reset_status() = 0;
};

namespace suanpan {
unique_ptr<Section> make_copy(const shared_ptr<Section>&);
unique_ptr<Section> make_copy(const unique_ptr<Section>&);
}

#endif

//! @}
