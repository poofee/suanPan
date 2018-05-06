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
 * @class CDPPS
 * @brief The CDPPS class.
 *
 * A 2D plane stress concrete material model that supports stiffness degradation.
 *
 * References:
 * 1. A Plastic-Damage Model for Concrete. https://doi.org/10.1016/0020-7683(89)90050-4
 * 2. Plastic-Damage Model for Cyclic Loading of Concrete Structures. https://doi.org/10.1061/(ASCE)0733-9399(1998)124:8(892)
 * 3. A Plastic-Damage Concrete Model for Earthquake Analysis of Dams. doi:10.1002/(SICI)1096-9845(199809)27:9<937::AID-EQE764>3.0.CO;2-5
 * 4. A Return-Mapping Algorithm for Plastic-Damage Models: 3-D and Plane Stress Formulation. doi:10.1002/1097-0207(20010120)50:2<487::AID-NME44>3.0.CO;2-N
 *
 * @author tlc
 * @date 12/12/2017
 * @version 0.2.0
 * @file CDPPS.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef CDPPS_H
#define CDPPS_H

#include <Material/Material3D/Material3D.h>

class CDPPS : public Material3D {
    static const double sqrt_three_over_two;
    static const mat unit_dev_tensor;
    static const uvec P;

    const double elastic_modulus, poissons_ratio, double_shear;
    const double a_t, cb_t, g_t, f_t, a_c, cb_c, g_c, f_c;
    const double alpha, alpha_p, s0;
    const double factor_a, factor_b, factor_c, tolerance;

    const vec unit_alpha_p;

    podarray<double> t_para, c_para;

    vec trial_full_strain, current_full_strain;
    vec trial_full_stress, current_full_stress;
    vec trial_plastic_strain, current_plastic_strain;

    mat initial_full_stiffness, inv_stiffness;

    bool elastic_flag;
    double r_weight;
    vec qkappah;
    rowvec drweight, dlambdadepsilon;
    mat jacobian, qkappaqsigma;

    static podarray<double> compute_backbone(double, double, double, double);
    static double r(const vec&);
    static rowvec dr(const vec&);
    inline double s(double) const;

    int compute_regular_return(const vec&);
    static mat form_stiffness(const mat&);

public:
    explicit CDPPS(unsigned = 0, // tag
        double = 3E4,            // elastic modulus
        double = .2,             // poissons ratio
        double = 3.,             // crack stress (+)
        double = 30.,            // crush stress (-)
        double = 5E-4,           // normalized crack energy (+)
        double = 1E-4,           // normalized crush energy (+)
        double = .2,             // hardening after crack stress a_t
        double = 5.,             // hardening after crush stress a_c
        double = .5,             // reference damage factor at half crack stress
        double = .5,             // reference damage factor at crush stress
        double = .2,             // dilatancy parameter
        double = 1.16,           // biaxial compression strength ratio
        double = .5,             // stiffness recovery
        double = 2400E-12        // density
    );

    void initialize(const shared_ptr<DomainBase>& = nullptr) override;

    unique_ptr<Material> get_copy() override;

    double get_parameter(const ParameterType&) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(const OutputType&) override;
};

#endif

//! @}
