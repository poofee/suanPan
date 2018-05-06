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

#include "Concrete02.h"
#include <Toolbox/utility.h>

const double Concrete02::crack_strain = 1E-4;
const double Concrete02::one_over_six = 1. / 6.;

Concrete02::Concrete02(const unsigned T, const double SP, const double R)
    : Material1D(T, MT_CONCRETE02, R)
    , peak_stress(SP > 0. ? -SP : SP)
    , peak_strain(-8.5E-4 * pow(-peak_stress, .25) + datum::eps * 1E10) // to avoid zero stiffness
    , tension_origin(0.) {}

void Concrete02::initialize(const shared_ptr<DomainBase>&) {
    initial_stiffness = current_stiffness = trial_stiffness = compute_compression_backbone(0.)(1);

    crack_stress = .33 * sqrt(-peak_stress);
    NT = 1.1;
    MT = trial_stiffness(0) * crack_strain / crack_stress;
    MC = 1. - 17.9 / peak_stress;
    NC = std::max(1., -peak_stress / 6.68 - 1.85);

    trial_history = current_history.zeros(18);
    current_flag.set_size(2);
    current_flag.fill(true);
    trial_flag = current_flag;
}

unique_ptr<Material> Concrete02::get_copy() { return make_unique<Concrete02>(*this); }

double Concrete02::get_parameter(const ParameterType& P) const {
    switch(P) {
    case ParameterType::DENSITY:
        return density;
    case ParameterType::ELASTICMODULUS:
    case ParameterType::YOUNGSMODULUS:
    case ParameterType::E:
        return initial_stiffness(0);
    default:
        return 0.;
    }
}

int Concrete02::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;
    incre_strain = trial_strain - current_strain;

    if(abs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& max_c_strain = trial_history(0);         // maximum compression strain recorded
    auto& max_t_strain = trial_history(1);         // maximum tension strain recorded
    auto& reverse_c_strain = trial_history(2);     // unloading point strain compression side
    auto& reverse_c_stress = trial_history(3);     // unloading point stress compression side
    auto& reverse_t_strain = trial_history(4);     // unloading point strain tension side
    auto& reverse_t_stress = trial_history(5);     // unloading point stress tension side
    auto& residual_c_strain = trial_history(6);    // residual strain in compression unloading path
    auto& residual_c_stiffness = trial_history(7); // tangent stiffness at residual
    auto& residual_t_strain = trial_history(8);    // residual strain in compression unloading path
    auto& residual_t_stiffness = trial_history(9); // tangent stiffness at residual
    auto& reload_c_strain = trial_history(10);
    auto& reload_c_stress = trial_history(11);
    auto& reload_c_stiffness = trial_history(12);
    auto& reload_t_strain = trial_history(13);
    auto& reload_t_stress = trial_history(14);
    auto& reload_t_stiffness = trial_history(15);
    auto& reverse_below_c_stress = trial_history(16); // unloading point stress compression side
    auto& reverse_below_t_stress = trial_history(17); // unloading point stress compression side

    trial_flag = current_flag;
    auto& on_compression_backbone = trial_flag(0);
    auto& on_tension_backbone = trial_flag(1);

    const auto strain_a = trial_strain(0) - residual_c_strain;
    const auto side = int(suanpan::sign(strain_a));
    const auto load_direction = int(suanpan::sign(incre_strain(0)));

    // step 1: is it tension or compression zone?
    // step 2: is it on backbone?
    // step 3: is it loading or unloading?

    auto update_residual = false;

    if(side < 0) {
        // the trial position is in compression zone
        // if current position is on backbone
        if(on_compression_backbone)
            // yes on backbone
            if(load_direction < 0) {
                // loading
                update_residual = true;
                const auto status = compute_compression_backbone(max_c_strain = trial_strain(0));
                trial_stress = status(0);
                trial_stiffness = status(1);
            } else {
                // unloading
                on_compression_backbone = false;
                // update reversing point
                reverse_c_strain = current_strain(0);
                reverse_c_stress = current_stress(0);
                trial_stiffness = reverse_c_stress / (reverse_c_strain - residual_c_strain);
                trial_stress = trial_stiffness * strain_a;
                // const auto status = compute_transition(trial_strain(0), reverse_c_strain, reverse_c_stress, initial_stiffness(0), residual_c_strain, 0., residual_c_stiffness);
                // trial_stress = status(0);
                // trial_stiffness = status(1);
            }
        else if(trial_strain(0) >= reverse_c_strain) {
            // still inside backbone
            trial_stiffness = reverse_c_stress / (reverse_c_strain - residual_c_strain);
            trial_stress = trial_stiffness * strain_a;
            // const auto status = compute_transition(trial_strain(0), reverse_c_strain, reverse_c_stress, initial_stiffness(0), residual_c_strain, 0., residual_c_stiffness);
            // trial_stress = status(0);
            // trial_stiffness = status(1);
        } else {
            // reload to backbone
            on_compression_backbone = true;
            const auto status = compute_compression_backbone(max_c_strain = trial_strain(0));
            trial_stress = status(0);
            trial_stiffness = status(1);
            update_residual = true;
        }

        if(update_residual) {
            // equation (3-150)
            const auto secant_stiffness = initial_stiffness(0) * (trial_stress(0) / initial_stiffness(0) / peak_strain + .57) / (trial_strain(0) / peak_strain + .57);
            // equation (3-154)
            residual_c_strain = trial_strain(0) - trial_stress(0) / secant_stiffness;
            // equation (3-151)
            residual_c_stiffness = .1 * initial_stiffness(0) * exp(-2. * trial_strain(0) / peak_strain);
            // equation (3-153)
            const auto status = compute_compression_backbone(reload_c_strain = trial_strain(0) + trial_strain(0) / (1.15 + 2.75 * trial_strain(0) / peak_strain));
            reload_c_stress = status(0);
            reload_c_stiffness = status(1);
        }
    } else {
        // if(tension_origin == 0.) tension_origin = residual_c_strain;

        // the trial position is in tension zone
        if(on_tension_backbone)
            // yes on backbone
            if(load_direction > 0) {
                // loading
                update_residual = true;
                const auto status = compute_tension_backbone(max_t_strain = trial_strain(0));
                trial_stress = status(0);
                trial_stiffness = status(1);
            } else {
                // unloading
                on_tension_backbone = false;
                reverse_t_strain = current_strain(0);
                reverse_t_stress = current_stress(0);
                trial_stiffness = reverse_t_stress / (reverse_t_strain - residual_c_strain);
                trial_stress = trial_stiffness * (trial_strain(0) - residual_c_strain);
            }
        else if(trial_strain(0) <= reverse_t_strain) {
            // not on backbone but still inside of backbone
            trial_stiffness = reverse_t_stress / (reverse_t_strain - residual_c_strain);
            trial_stress = trial_stiffness * strain_a;
        } else {
            // reloading from unloading path
            on_tension_backbone = true;
            const auto status = compute_tension_backbone(max_t_strain = trial_strain(0));
            trial_stress = status(0);
            trial_stiffness = status(1);
            update_residual = true;
        }

        if(update_residual) {
            const auto offset = trial_strain(0) - tension_origin;
            // equation (3-160)
            const auto secant_stiffness = initial_stiffness(0) * (trial_stress(0) / initial_stiffness(0) / crack_strain + .67) / (offset / crack_strain + .67);
            // equation (3-164)
            residual_t_strain = trial_strain(0) - trial_stress(0) / secant_stiffness;
            // equation (3-161)
            residual_t_stiffness = initial_stiffness(0) / (pow(abs(offset / crack_strain), 1.1) + 1.);
            // equation (3-163)
            const auto status = compute_tension_backbone(reload_t_strain = 1.22 * trial_strain(0));
            reload_t_stress = status(0);
            reload_t_stiffness = status(1);
        }
    }

    // if the stiffness is too small, set to a larger value to avoid instability
    if(abs(trial_stiffness(0)) < 1E-12) trial_stiffness = 1E-8 * suanpan::sign(trial_stiffness(0));

    return SUANPAN_SUCCESS;
}

int Concrete02::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
    if(update_trial_status(t_strain) == SUANPAN_SUCCESS) {
        const auto factor = 3.5E-2 * peak_stress * peak_stress;
        const auto magnification_factor = (1. + pow(abs(t_strain_rate(0) / factor), one_over_six)) / (1. + pow(abs(1E-5 / factor), one_over_six));
        trial_stress *= magnification_factor;
        trial_stiffness *= magnification_factor;
        return SUANPAN_SUCCESS;
    }
    return SUANPAN_FAIL;
}

int Concrete02::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history.zeros();
    current_flag.fill(true);
    current_stiffness = initial_stiffness;
    return reset_status();
}

int Concrete02::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_flag = trial_flag;
    current_stiffness = trial_stiffness;
    return 0;
}

int Concrete02::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_flag = current_flag;
    trial_stiffness = current_stiffness;
    return 0;
}

podarray<double> Concrete02::compute_compression_backbone(const double t_strain) const {
    podarray<double> status(2);

    const auto normal_strain = t_strain / peak_strain;

    const auto tmp_a = pow(normal_strain, NC);
    const auto tmp_b = NC == 1. ? 1. + (MC - 1. + log(normal_strain)) * normal_strain : 1. + (MC - NC / (NC - 1.)) * normal_strain + tmp_a / (NC - 1.);
    status(0) = peak_stress * MC * normal_strain / tmp_b;
    status(1) = peak_stress / peak_strain * MC * (1. - tmp_a) / tmp_b / tmp_b;

    return status;
}

podarray<double> Concrete02::compute_tension_backbone(const double t_strain) const {
    podarray<double> status(2);

    const auto normal_strain = (t_strain - tension_origin) / crack_strain;

    const auto tmp_a = pow(normal_strain, NT);
    const auto tmp_b = NT == 1. ? 1. + (MT - 1. + log(normal_strain)) * normal_strain : 1. + (MT - NT / (NT - 1.)) * normal_strain + tmp_a / (NT - 1.);
    status(0) = crack_stress * MT * normal_strain / tmp_b;
    status(1) = crack_stress / crack_strain * MT * (1. - tmp_a) / tmp_b / tmp_b;

    return status;
}

/**
 * \brief compute the transition curve between two given points
 * \param TX new_x
 * \param XS x_start
 * \param YS y_start
 * \param ES tagent_start
 * \param XF x_final
 * \param YF y_final
 * \param EF tagent_final
 * \return response[new_y,new_tagent]
 */
podarray<double> Concrete02::compute_transition(const double TX, const double XS, const double YS, const double ES, const double XF, const double YF, const double EF) {
    podarray<double> response(2);

    const auto TA = XF - XS;
    const auto ESEC = (YF - YS) / TA; // eq. 3-119
    const auto TB = ESEC - ES;
    const auto R = (EF - ESEC) / TB; // eq. 3-122
    const auto A = TB / pow(abs(TA), R);
    const auto TC = TX - XS;
    const auto TD = A * pow(abs(TC), R);

    response(0) = YS + TC * (ES + TD); // eq. 3-120
    response(1) = ES + (R + 1.) * TD;  // eq. 3-121

    return response;
}
