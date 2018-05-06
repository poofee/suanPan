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
 * @class S4
 * @brief A S4 class.
 *
 * @author tlc
 * @date 08/02/2018
 * @version 0.1.0
 * @file S4.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef S4_H
#define S4_H

#include <Element/MaterialElement.h>

class S4 : public MaterialElement {
    struct IntegrationPoint {
        vec coor;
        const double factor;
        unique_ptr<Material> s_material;
        mat BM, BP;
        IntegrationPoint(vec, double, unique_ptr<Material>&&);
    };

    static const unsigned s_node, s_dof, s_size;

    static const uvec m_dof, p_dof;

    static mat iso_mapping;

    mat trans_mat;

    const double thickness;

    const char int_scheme;

    vector<IntegrationPoint> int_pt;

    mat direction_cosine();
    static void transform_from_local_to_global(vec&, const mat&);
    static void transform_from_local_to_global(mat&, const mat&);
    static vec reshuffle(const vec&, const vec&);
    static mat reshuffle(const mat&, const mat&);

public:
    S4(unsigned, const uvec&, unsigned, double = 1., char = 'I');

    void initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    vector<vec> record(const OutputType&) override;

    void print() override;
};

#endif

//! @}
