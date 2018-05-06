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
 * @class CP8
 * @brief The CP8 class handles CPS8, CPE8, CPS8R and CPE8R elements. It is a four node constant strain membrane element with optional reduced integration for both plane stress and plane strain problems.
 * @author tlc
 * @date 12/09/2017
 * @version 0.1.2
 * @file CP8.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef CP8_H
#define CP8_H

#include <Element/MaterialElement.h>

class CP8 final : public MaterialElement {
    struct IntegrationPoint {
        vec coor;
        double weight, jacob_det;
        unique_ptr<Material> m_material;
        mat pn_pxy, strain_mat;
        IntegrationPoint(vec, double, double, unique_ptr<Material>&&, mat);
    };

    static const unsigned m_node, m_dof, m_size;

    const double thickness;

    const bool reduced_scheme;

    vector<IntegrationPoint> int_pt;

public:
    CP8(unsigned,      // tag
        const uvec&,   // node tags
        unsigned,      // material tag
        double = 1.,   // thickness
        bool = false,  // reduced integration
        bool = false); // nonlinear geometry switch

    void initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
