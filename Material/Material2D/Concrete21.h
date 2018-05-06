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
 * @class Concrete21
 * @brief A Concrete21 material class.
 *
 * For plain concrete. Only monotonic load is supported at moment. Cyclic behaviour does not converge in most cases.
 *
 * @author tlc
 * @date 04/04/2018
 * @version 0.1.2
 * @file Concrete21.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef CONCRETE21_H
#define CONCRETE21_H

#include <Material/Material1D/Concrete01.h>
#include <Material/Material2D/Material2D.h>

class Concrete21 : public Material2D {
    Concrete01 concrete_major, concrete_minor;

    const bool degrade = true;

    const double poissons_ratio = .2, peak_strain;

    mat poissons_mat;

    double elastic_modulus = 0.;
    double shear_modulus = 0.;

public:
    Concrete21(unsigned, // tag
        double,          // peak stress in negative
        BackboneType,    // backbone type
        bool = false,    // center oriented or using unloading criterion
        double = 1E-2,   // factrue energy
        double = 0.,     // density
        PlaneType = PlaneType::S);

    void initialize(const shared_ptr<DomainBase>& = nullptr) override;

    unique_ptr<Material> get_copy() override;

    double get_parameter(const ParameterType&) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
