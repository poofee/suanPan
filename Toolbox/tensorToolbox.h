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

#ifndef TENSORTOOLBOX_H
#define TENSORTOOLBOX_H

#include <suanPan.h>

namespace tensor {
vec unit_tensor2();
mat unit_deviatoric_tensor4();
mat unit_deviatoric_tensor4v2();
mat unit_symmetric_tensor4();

// applies to principal tensor
namespace stress {
    double invariant1(const vec&);
    double invariant2(const vec&);
    double invariant3(const vec&);

    double lode(const vec&);
}
namespace strain {
    double invariant1(const vec&);
    double invariant2(const vec&);
    double invariant3(const vec&);

    double lode(const vec&);
}
double trace(const vec&);
double mean(const vec&);
vec dev(const vec&);

mat dev(const mat&);

namespace strain {
    mat to_tensor(const vec&);
    vec to_voigt(const mat&);
    double norm(const vec&);
}
namespace stress {
    mat to_tensor(const vec&);
    vec to_voigt(const mat&);
    double norm(const vec&);
}
}

namespace transform {
double atan2(const vec&);
mat compute_jacobian_nominal_to_principal(const mat&);
mat compute_jacobian_principal_to_nominal(const mat&);
vec haigh_westergaard_to_principal(double, double, double);
namespace strain {
    double angle(const vec&);
    mat trans(double);
    vec principal(const vec&);
    vec rotate(const vec&, double);
}
namespace stress {
    double angle(const vec&);
    mat trans(double);
    vec principal(const vec&);
    vec rotate(const vec&, double);
}
namespace beam {
    mat global_to_local(double, double, double, double);
    mat global_to_local(double, double, double);
    mat global_to_local(const vec&, const double);
}
}

namespace suanpan {
template <typename T> T ramp(const T in) { return in > T(0) ? in : T(0); }
}

#endif
