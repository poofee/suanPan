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
 * @class QE2
 * @brief A QE2 class.
 *
 * [10.1002/nme.1620381102](http://onlinelibrary.wiley.com/doi/10.1002/nme.1620381102/full)
 *
 * @author tlc
 * @date 13/09/2017
 * @version 0.2.0
 * @file QE2.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef QE2_H
#define QE2_H

#include <Element/MaterialElement.h>

class QE2 : public MaterialElement {
    struct IntegrationPoint {
        vec coor;
        const double factor;
        unique_ptr<Material> m_material;
        mat P, A, B, BI;
        IntegrationPoint(vec, double, unique_ptr<Material>&&);
    };

    static const unsigned m_node, m_dof, m_size;

    static mat mapping;

    const double thickness;

    vector<IntegrationPoint> int_pt;

    mat trans_mat; // temporaty factor matrix used to recover stress and strain

    mat HT, HIL, HILI, L, LI; // constant matrices

    mat initial_qtitt, trial_qtitt, current_qtitt; // eq. 65
    vec trial_qtifi, current_qtifi;                // eq. 65
    vec trial_lambda, current_lambda;              // enhanced strain
    vec trial_alpha, current_alpha;                // strain
    vec trial_beta, current_beta;                  // stress

    vec pre_disp;

public:
    QE2(unsigned, const uvec&, unsigned, double = 1.);

    void initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(const OutputType&) override;

    void print() override;
};

#endif

//! @}
