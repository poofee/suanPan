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
 * @class MooneyRivlin
 * @brief The MooneyRivlin class defines a Mooney-Rivlin hyperelastic material for 3-D problems.
 *
 * @author tlc
 * @date 21/01/2018
 * @version 0.1.0
 * @file MooneyRivlin.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef MOONEYRIVLIN_H
#define MOONEYRIVLIN_H

#include <Material/Material3D/Material3D.h>

class MooneyRivlin : public Material3D {
    static const vec I1E;
    static const mat I2EE;

    double K; /**< bulk modulus */
    double A10, A01;

public:
    explicit MooneyRivlin(unsigned = 0, // tag
        double = 2E5,                   // bulk modulus
        double = .2,                    // a10
        double = .2,                    // a01
        double = 0.                     // density
    );

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
