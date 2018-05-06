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
 * @fn MatrixModifier
 * @brief The MatrixModifier class.
 *
 * @author tlc
 * @date 27/10/2017
 * @version 0.1.0
 * @file MatrixModifier.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef MATRIXMODIFIER_H
#define MATRIXMODIFIER_H

#include <Element/Element.h>
#include <suanPan.h>

namespace suanpan {
namespace mass {
    struct lumped_simple {
        template <typename T> static Op<Col<T>, op_diagmat> apply(const Mat<T>&);
    };

    template <typename T> Op<Col<T>, op_diagmat> lumped_simple::apply(const Mat<T>& mass) { return diagmat(sum(mass)); }

    struct lumped_scale {
        template <typename T> static Op<Col<T>, op_diagmat> apply(const Mat<T>&, unsigned);
    };

    template <typename T> Op<Col<T>, op_diagmat> lumped_scale::apply(const Mat<T>& mass, const unsigned dim) {
        Col<T> diag_mass(mass.n_rows);

        for(unsigned I = 0; I < dim; ++I) {
            auto total_mass = 0.;
            auto true_mass = 0.;
            for(auto J = I; J < diag_mass.n_elem; J += dim) {
                true_mass += mass(J, J);
                total_mass += sum(mass.row(J));
            }
            auto factor = total_mass / true_mass;
            for(auto J = I; J < diag_mass.n_elem; J += dim) diag_mass(J) = mass(J, J) * factor;
        }

        return diagmat(diag_mass);
    }
}
namespace damping {
    struct rayleigh {
        template <typename T> static void apply(const shared_ptr<Element>&, double, double);
    };

    template <typename T> void rayleigh::apply(const shared_ptr<Element>& element_obj, const double alpha, const double beta) { access::rw(element_obj->get_damping()) = alpha * element_obj->get_stiffness() + beta * element_obj->get_mass(); }
}
}

#endif

//! @}
