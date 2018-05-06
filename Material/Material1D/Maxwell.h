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
 * @class Maxwell
 * @brief A 1-D Maxwell material class.
 * @author tlc
 * @date 08/03/2018
 * @version 0.2.0
 * @file Maxwell.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef MAXWELL_H
#define MAXWELL_H

#include <Material/Material1D/Material1D.h>
#include <Solver/ODE_Solver/ODE.h>
#include <Solver/ODE_Solver/ODE_Solver.h>

class Maxwell : public Material1D {
    const double elastic_modulus;

    const double* incre_time = nullptr;

    struct Damper : ODE {
        const double elastic_modulus, alpha, damping;
        double c_strain_rate = 0., c_strain_acceleration = 0.;
        Damper(double, double, double);
        unique_ptr<ODE> get_copy() override;
        vec eval(double, const vec&) override;
        void set_current_status(const vec&, const vec&);
    };

    unique_ptr<ODE> viscosity;
    unique_ptr<ODE_Solver> solver;

    Damper* damper_ptr = dynamic_cast<Damper*>(viscosity.get());

public:
    explicit Maxwell(unsigned, // tag
        double,                // elastic modulus
        double,                // alpha
        double                 // damping
    );
    Maxwell(const Maxwell&);

    void initialize(const shared_ptr<DomainBase>& = nullptr) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&, const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
