////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2018 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "MooneyRivlin.h"
#include <Toolbox/tensorToolbox.h>

const vec MooneyRivlin::I1E{ std::initializer_list<double>{ 2., 2., 2., 0., 0., 0. } };
const mat MooneyRivlin::I2EE;

MooneyRivlin::MooneyRivlin(const unsigned T, const double KK, const double AA, const double AB, const double R)
    : Material3D(T, MT_MOONEYRIVLIN, R)
    , K(KK)
    , A10(AA)
    , A01(AB) {}

void MooneyRivlin::initialize(const shared_ptr<DomainBase>&) {
    if(I2EE.is_empty()) {
        mat TI2EE(6, 6, fill::zeros);
        TI2EE(span(0, 2), span(0, 2)) = 4.;
        for(uword I = 0; I < 3; ++I) TI2EE(I) = 0.;
        for(uword I = 3; I < 6; ++I) TI2EE(I) = -2.;
        access::rw(I2EE) = TI2EE;
    }

    const vec zero_strain{ std::initializer_list<double>{ 1., 0., 0., 0., 1., 0., 0., 0., 1. } };
    update_trial_status(zero_strain);
    current_stiffness = initial_stiffness = trial_stiffness;
}

unique_ptr<Material> MooneyRivlin::get_copy() { return make_unique<MooneyRivlin>(*this); }

int MooneyRivlin::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;

    const mat gradient(trial_strain.memptr(), 3, 3, false);

    const mat deformation_tensor = gradient.t() * gradient;

    const auto& C1 = deformation_tensor(0, 0);
    const auto& C2 = deformation_tensor(1, 1);
    const auto& C3 = deformation_tensor(2, 2);
    const auto& C4 = deformation_tensor(0, 1);
    const auto& C5 = deformation_tensor(1, 2);
    const auto& C6 = deformation_tensor(0, 2);

    const auto I1 = C1 + C2 + C3;
    const auto I2 = C1 * C2 + C1 * C3 + C2 * C3 - C4 * C4 - C5 * C5 - C6 * C6;
    const auto I3 = det(deformation_tensor);

    const auto J3 = sqrt(I3);
    const auto J3M1 = J3 - 1.;

    vec I2E{ std::initializer_list<double>{ C2 + C3, C3 + C1, C1 + C2, -C4, -C5, -C6 } };
    I2E *= 2.;
    vec I3E{ std::initializer_list<double>{ C2 * C3 - C5 * C5, C3 * C1 - C6 * C6, C1 * C2 - C4 * C4, C5 * C6 - C3 * C4, C6 * C4 - C1 * C5, C4 * C5 - C2 * C6 } };
    I3E *= 2.;

    auto W1 = pow(I3, -1. / 3.);
    auto W2 = I1 * pow(I3, -4. / 3.) / 3.;
    auto W3 = pow(I3, -2. / 3.);
    auto W4 = 2. / 3. * I2 * pow(I3, -5. / 3.);
    auto W5 = .5 * pow(I3, -.5);

    const auto J1E = W1 * I1E - W2 * I3E;
    const auto J2E = W3 * I2E - W4 * I3E;
    const auto J3E = W5 * I3E;

    trial_stress = A10 * J1E + A01 * J2E + K * J3M1 * J3E;

    mat I3EE(6, 6, fill::zeros);
    I3EE(0, 1) = I3EE(1, 0) = -2. * (I3EE(3, 3) = -2. * C3);
    I3EE(0, 2) = I3EE(2, 0) = -2. * (I3EE(5, 5) = -2. * C2);
    I3EE(1, 2) = I3EE(2, 1) = -2. * (I3EE(4, 4) = -2. * C1);
    I3EE(0, 4) = I3EE(4, 0) = -2. * (I3EE(3, 5) = I3EE(5, 3) = 2. * C5);
    I3EE(1, 5) = I3EE(5, 1) = -2. * (I3EE(3, 4) = I3EE(4, 3) = 2. * C6);
    I3EE(2, 3) = I3EE(3, 2) = -2. * (I3EE(4, 5) = I3EE(5, 4) = 2. * C4);

    const auto W6 = pow(I3, -2. / 3.);
    const auto W8 = pow(I3, -.5);
    const auto W9 = .5 * W8;
    const auto W7 = .75 * (W5 = 8. / 9. * I2 * pow(I3, -5. / 3.));
    W4 = 2. * (W1 = 2. / 3. * W8);
    W3 = .375 * (W2 = 8. / 9. * I1 * pow(I3, -4. / 3.));

    const vec T1 = (A10 * W1 * J1E + A01 * W4 * J2E) * J3E.t();
    trial_stiffness = -T1 - T1.t() + A01 * W6 * I2EE + (A01 * W5 + A10 * W2 + K - K * J3M1 * W8) * J3E * J3E.t() + (K * J3M1 * W9 - A10 * W3 - A01 * W7) * I3EE;

    return 0;
}

int MooneyRivlin::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    trial_strain.zeros();
    trial_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return 0;
}

int MooneyRivlin::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return 0;
}

int MooneyRivlin::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return 0;
}

void MooneyRivlin::print() { suanpan_info("Mooney-Rivlin Material.\n"); }
