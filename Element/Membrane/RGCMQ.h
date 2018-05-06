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
 * @class RGCMQ
 * @brief A RGCMQ class.
 *
 * @author tlc
 * @date 08/02/2018
 * @version 0.1.0
 * @file RGCMQ.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef RGCMQ_H
#define RGCMQ_H

#include <Element/MaterialElement.h>
#include <Toolbox/IntegrationPlan.h>

class RGCMQ : public MaterialElement {
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
        unique_ptr<Material> rebar_x = nullptr, rebar_y = nullptr;
        mat poly_stress, poly_strain;
        IntegrationPoint(vec, double, unique_ptr<Material>&&);
    };

    static const unsigned m_node, m_dof, m_size;

    static mat iso_mapping;

    const double thickness;

    const char int_scheme;

    const double rho_x, rho_y;
    const unsigned mat_x, mat_y;

    const bool reinforced_x = rho_x != 0. && mat_x != 0;
    const bool reinforced_y = rho_y != 0. && mat_y != 0;

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
    RGCMQ(unsigned, const uvec&, unsigned, double = 1., char = 'I', double = 0., unsigned = 0, double = 0., unsigned = 0);

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
