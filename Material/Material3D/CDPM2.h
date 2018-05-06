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
 * @class CDPM2
 * @brief The CDPM2 class.
 *
 * References:
 * 1. CDPM2: A Damage-Plasticity Approach to Modelling the Failure of Concrete. https://doi.org/10.1016/j.ijsolstr.2013.07.008
 * 2. Damage-Plastic Model for Concrete Failure. https://doi.org/10.1016/j.ijsolstr.2006.06.032
 *
 * @author tlc
 * @date 20/01/2018
 * @version 0.2.0
 * @file CDPM2.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef CDPM2_H
#define CDPM2_H

#include <Material/Material3D/Material3D.h>

class CDPM2 : public Material3D {
    static const double sqrt_six;
    static const mat unit_dev_tensor;
    static const rowvec unit_tensor2;

    const double elastic_modulus, poissons_ratio, double_shear, bulk_modulus;

    const double tolerance = 1E-12;

    const double ft = 3, fc = 30.;

    const double ft_over_fc = ft / fc, fc_fc = fc * fc, sqrt_six_fc = sqrt_six * fc;

    double rbc = 1.16, e = ((fc_fc - ft * ft) * rbc + fc * ft * (rbc * rbc - 1.)) / (2. * rbc * (fc_fc - ft * ft) + fc * ft * (1. - rbc * rbc)), ea = 1. - e * e, eb = 2. * e - 1, ec = eb * eb, ed = 5. * e * e - 4. * e, ah = 8E-2, bh = 3E-3, ch = 2., dh = 1E-6, eh = bh - dh, fh = eh * ch / (ah - bh), gh = (bh - ah) / ch / fc, qh0 = .3, one_minus_qh0 = 1. - qh0, hp = .01, m0 = 3. * (fc / ft - ft / fc) * e / (1. + e), half_m0 = .5 * m0, df = .85, log_df = log((df + 1) / (2. * df - 1.)), sqrt_df = ft / sqrt(1.5 + 3. * df * df), ratio_factor = 3. * (1. - 2. * poissons_ratio) / (1. + poissons_ratio);

    double epsilon0 = ft / elastic_modulus, as = 5., one_minus_as = 1. - as;

    double rft1 = .3, rwf1 = .15, ft1 = rft1 * ft, gft = 5E-2, wf = 2. * (rft1 + rwf1) * gft / ft, wf1 = rwf1 * wf, epsilonft = 5E-5, epsilonfc = 1E-5, h = 100.;

    const rowvec dsigmadepsilon;

    vec trial_plastic_strain, current_plastic_strain;

    mat inv_stiffness;

    podarray<double> compute_plasticity(double, double, double, double) const;
    mat compute_jacobian(const podarray<double>&, double) const;
    vec compute_residual(const podarray<double>&, const vec&, const vec&) const;
    double compute_yield_function(double, double, double, double) const;
    inline double compute_incre_kappap(double, double, double) const;
    inline double compute_kaih(double) const;
    double compute_plastic_potential_ratio(double, double, double) const;
    bool compute_regular_return(vec&, mat&, podarray<double>&, const vec&, double) const;
    bool compute_vertex_return(double&, double, double&) const;
    podarray<double> compute_damage(double, double, double) const;

public:
    explicit CDPM2(unsigned = 0, // tag
        double = 3E4,            // elastic modulus
        double = .2,             // poissons ratio
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
