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
 * @class Bilinear1D
 * @brief A Bilinear1D material class.
 * @author tlc
 * @date 08/08/2017
 * @version 0.1.0
 * @file Bilinear1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef BILINEAR1D_H
#define BILINEAR1D_H

#include <Material/Material1D/Material1D.h>

struct DataBilinear1D {
    double elastic_modulus; /**< elastic modulus */
    double yield_stress;    /**< initial yield stress */
    double hardening_ratio; /**< hardening ratio */
    double beta;            /**< isotropic/kinematic hardening factor */
    double plastic_modulus; /**< plastic modulus */
    double tolerance;       /**< tolerance */
};

class Bilinear1D final : DataBilinear1D, public Material1D {
public:
    explicit Bilinear1D(unsigned = 0, /**< tag */
        double = 2E5,                 /**< elastic modulus */
        double = 400.,                /**< initial yield stress */
        double = .05,                 /**< hardening ratio */
        double = 0.,                  /**< isotropic/kinematic hardening factor */
        double = 0.);                 /**< density */

    void initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
