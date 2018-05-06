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
 * @class GCMT
 * @brief A GCMT class.
 *
 * @author tlc
 * @date 01/03/2018
 * @version 0.4.0
 * @file GCMT.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef GCMT_H
#define GCMT_H

#include <Element/MaterialElement.h>
#include <Toolbox/IntegrationPlan.h>

class GCMT : public MaterialElement {
    class ResultantConverter {
        mat converter_a, converter_b;
        mat direction_cosine;

    public:
        ResultantConverter(unsigned, double, const vector<weak_ptr<Node>>&, const IntegrationPlan&, const mat&);
        double F(const vec&) const;
        double V(const vec&) const;
        double M(const vec&) const;
    };

    struct IntegrationPoint {
        vec coor;
        const double factor;
        unique_ptr<Material> m_material;
        mat poly_stress, poly_strain;
        IntegrationPoint(vec, double, unique_ptr<Material>&&);
    };

    static const unsigned m_node, m_dof, m_size;

    static mat mapping;

    const double thickness;

    const char int_scheme;

    vector<IntegrationPoint> int_pt;
    vector<ResultantConverter> edge;

    mat trans_mat; // temporaty factor matrix used to recover stress and strain

    mat HT, NT, MT, N, M; // constant matrices

    mat initial_viwt, trial_viwt, current_viwt;
    vec trial_vif, current_vif;
    vec trial_alpha, current_alpha; // stress
    vec trial_beta, current_beta;   // strain
    vec trial_zeta, current_zeta;   // enhanced strain

    vec pre_disp;

    static mat form_interpolation_displacement(const mat&, const mat&);
    static mat form_interpolation_enhanced_strain(const mat&);

public:
    GCMT(unsigned, const uvec&, unsigned, double = 1., char = 'I');

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
