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
 * @class Bilinear2D
 * @brief A Bilinear2D class.
 * @author tlc
 * @date 04/10/2017
 * @version 0.1.0
 * @file Bilinear2D.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef BILINEAR2D_H
#define BILINEAR2D_H

#include <Material/Material2D/Material2D.h>
#include <Material/Material3D/Bilinear3D.h>

using std::array;

class Bilinear2D : public Material2D {
    static const array<unsigned, 3> F;

    vec trial_full_strain;

    Bilinear3D base;

public:
    explicit Bilinear2D(unsigned = 0, /**< tag */
        double = 2E5,                 /**< elastic modulus */
        double = .25,                 /**< poisson's ratio */
        double = 400.,                /**< initial yield stress */
        double = .05,                 /**< hardening ratio */
        double = 0.,                  /**< isotropic/kinematic hardening factor */
        PlaneType = PlaneType::S,     /**< plane stress or plane strain */
        double = 0.                   /**< density */
    );

    void initialize(const shared_ptr<DomainBase>&) override;

    double get_parameter(const ParameterType& = ParameterType::DENSITY) const override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
