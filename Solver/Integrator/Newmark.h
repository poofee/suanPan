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
 * @class Newmark
 * @brief A Newmark class defines a solver using Newmark algorithm.
 *
 * `Newmark` algorithm is unconditionally stable if
 * \f{gather}{\alpha\geq\dfrac{1}{4}\left(\dfrac{1}{2}+\beta\right)^2,\qquad\beta\geq\dfrac{1}{2}\f}.
 *
 * There are several choices for solver parameters.
 *
 * Constant Acceleration:
 * \f{gather}{\alpha=\dfrac{1}{4},\qquad\beta=\dfrac{1}{2}\f}.
 *
 * Linear Acceleration:
 * \f{gather}{\alpha=\dfrac{1}{6},\qquad\beta=\dfrac{1}{2}\f}.
 *
 * @author tlc
 * @date 25/08/2017
 * @version 0.1.1
 * @file Newmark.h
 * @addtogroup Integrator
 * @{
 */

#ifndef NEWMARK_H
#define NEWMARK_H

#include "Integrator.h"

class Newmark final : public Integrator {
    const double alpha; /**< parameter = .25 */
    const double beta;  /**< parameter = .5 */

    double DT = 0.; /**< previous incremental time */

    double C0 = 0., C1 = 0., C2 = 0., C3 = 0., C4 = 0., C5 = 0., C6 = 0., C7 = 0.; /**< parameters */

    void update_parameter();

public:
    explicit Newmark(const unsigned& = 0, const double& = .25, const double& = .5);

    void assemble_resistance() override;
    void assemble_matrix() override;

    void commit_status() const override;

    void print() override;
};

#endif

//! @}
