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
 * @class RebarLayer
 * @brief A RebarLayer class.
 *
 * @author tlc
 * @date 10/01/2018
 * @version 0.1.0
 * @file RebarLayer.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef REBARLAYER_H
#define REBARLAYER_H

#include <Element/MaterialElement.h>

class RebarLayer : public MaterialElement {
    struct IntegrationPoint {
        vec coor;
        double weight, jacob_det;
        unique_ptr<Material> m_material;
        mat pn_pxy;
        mat BN, BG;
        IntegrationPoint(vec, double, double, unique_ptr<Material>&&, mat);
    };

    static const unsigned m_node, m_dof, m_size;

    vector<IntegrationPoint> int_pt;

    const double d_pos_a, d_pos_b, d_rho_a, d_rho_b, d_area_a, d_area_b;

public:
    RebarLayer(unsigned, const uvec&, unsigned, double = 1., bool = true);

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
