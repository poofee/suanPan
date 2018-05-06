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

#include "Concrete01.h"
#include <Toolbox/utility.h>

const double Concrete01::one_over_six = 1. / 6.;

Concrete01::Concrete01(const unsigned T, const double SP, const BackboneType TP, const bool CO, const bool NT, const TensionType TT, const double FE, const double R)
    : Material1D(T, MT_CONCRETE01, R)
    , backbone_type(TP)
    , tension_type(TT)
    , no_residual(CO)
    , no_tension(NT)
    , peak_stress(SP > 0. ? -SP : SP)
    , peak_strain(-8.5E-4 * pow(-peak_stress, .25) + datum::eps * 1E10) // to avoid zero stiffness
    , facture_energy(FE)
    , tension_origin(0.)
    , confinement(1.)
    , crack_stress(.31 * sqrt(-peak_stress)) {}

void Concrete01::initialize(const shared_ptr<DomainBase>&) {
    switch(backbone_type) {
    case BackboneType::POPOVICS: {
        N = 1. - .05802 * peak_stress;
        break;
    }
    case BackboneType::THORENFELDT: {
        M = .67 - peak_stress / 62.;
        N = .85 - peak_stress / 17.;
        break;
    }
    case BackboneType::TSAI: {
        M = 1. - 17.9 / peak_stress;
        N = std::max(1., -peak_stress / 6.68 - 1.85);
        break;
    }
    case BackboneType::KPS:
        break;
    }

    compute_compression_backbone(0.);

    initial_stiffness = current_stiffness = trial_stiffness;

    trial_history = current_history.zeros(7);
    current_flag.set_size(2);
    current_flag.fill(true);
    trial_flag = current_flag;

    crack_strain = crack_stress / initial_stiffness(0);
    ultimate_strain = 2. * facture_energy / crack_stress;
}

unique_ptr<Material> Concrete01::get_copy() { return make_unique<Concrete01>(*this); }

double Concrete01::get_parameter(const ParameterType& P) const {
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

bool Concrete01::is_cracked() {
    const auto& max_t_strain = trial_history(1);
    return max_t_strain > crack_strain;
}

int Concrete01::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(abs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& max_c_strain = trial_history(0);     // maximum compression strain recorded
    auto& max_t_strain = trial_history(1);     // maximum tension strain recorded
    auto& reverse_c_strain = trial_history(2); // unloading point strain compression side
    auto& reverse_c_stress = trial_history(3); // unloading point stress compression side
    auto& reverse_t_strain = trial_history(4); // unloading point strain tension side
    auto& reverse_t_stress = trial_history(5); // unloading point stress tension side
    auto& residual_strain = trial_history(6);  // residual strain in compression unloading path

    trial_flag = current_flag;
    auto& on_compression_backbone = trial_flag(0);
    auto& on_tension_backbone = trial_flag(1);

    const auto direction = int(suanpan::sign(incre_strain(0)));
    const auto offset = trial_strain(0) - residual_strain;
    const auto side = int(suanpan::sign(offset));

    // step 1: is it tension or compression zone?
    // step 2: is it on backbone?
    // step 3: is it loading or unloading?

    if(side < 0) {
        auto update_residual = false;

        // the trial position is in compression zone
        // if current position is on backbone
        if(on_compression_backbone)
            // yes on backbone
            if(direction < 0) {
                // loading
                update_residual = true;
                compute_compression_backbone(max_c_strain = trial_strain(0));
            } else {
                // unloading
                on_compression_backbone = false;
                // update reversing point
                reverse_c_strain = current_strain(0);
                reverse_c_stress = current_stress(0);
                trial_stiffness = reverse_c_stress / (reverse_c_strain - residual_strain);
                trial_stress = trial_stiffness * offset;
            }
        else if(trial_strain(0) >= reverse_c_strain) {
            // still inside backbone
            trial_stiffness = reverse_c_stress / (reverse_c_strain - residual_strain);
            trial_stress = trial_stiffness * offset;
        } else {
            // reload to backbone
            on_compression_backbone = true;
            compute_compression_backbone(max_c_strain = trial_strain(0));
            update_residual = true;
        }

        if(update_residual && !no_residual) {
            // equation (3-150)
            const auto secant_stiffness = initial_stiffness(0) * (trial_stress(0) / initial_stiffness(0) / peak_strain + .57) / (trial_strain(0) / peak_strain + .57);
            // equation (3-154)
            residual_strain = trial_strain(0) - trial_stress(0) / secant_stiffness;
        }
    } else {
        if(no_tension) {
            trial_stress = 0.;
            trial_stiffness = 0.;
            return SUANPAN_SUCCESS;
        }

        // if(tension_origin == 0.) tension_origin = residual_c_strain;

        // the trial position is in tension zone
        if(on_tension_backbone)
            // yes on backbone
            if(direction > 0) {
                // loading
                compute_tension_backbone(max_t_strain = trial_strain(0));
            } else {
                // unloading
                on_tension_backbone = false;
                reverse_t_strain = current_strain(0);
                reverse_t_stress = current_stress(0);
                trial_stiffness = reverse_t_stress / (reverse_t_strain - residual_strain);
                trial_stress = trial_stiffness * (trial_strain(0) - residual_strain);
            }
        else if(trial_strain(0) <= reverse_t_strain) {
            // not on backbone but still inside of backbone
            trial_stiffness = reverse_t_stress / (reverse_t_strain - residual_strain);
            trial_stress = trial_stiffness * offset;
        } else {
            // reloading from unloading path
            on_tension_backbone = true;
            compute_tension_backbone(max_t_strain = trial_strain(0));
        }
    }

    // keep stiffness as a nonzero number
    if(abs(trial_stiffness(0)) < 1E-10) trial_stiffness(0) = (suanpan::sign(trial_stiffness(0)) == 0. ? 1. : suanpan::sign(trial_stiffness(0))) * 1E-8;

    return SUANPAN_SUCCESS;
}

int Concrete01::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
    if(update_trial_status(t_strain) == SUANPAN_SUCCESS) {
        const auto factor = 3.5E-2 * peak_stress * peak_stress;
        const auto magnification_factor = (1. + pow(abs(t_strain_rate(0) / factor), one_over_six)) / (1. + pow(abs(1E-5 / factor), one_over_six));
        trial_stress *= magnification_factor;
        trial_stiffness *= magnification_factor;
        return SUANPAN_SUCCESS;
    }
    return SUANPAN_FAIL;
}

int Concrete01::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history.zeros();
    current_flag.fill(true);
    current_stiffness = initial_stiffness;
    return reset_status();
}

int Concrete01::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_flag = trial_flag;
    current_stiffness = trial_stiffness;
    return 0;
}

int Concrete01::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_flag = current_flag;
    trial_stiffness = current_stiffness;
    return 0;
}

void Concrete01::compute_compression_backbone(const double t_strain) {
    const auto normal_strain = t_strain / peak_strain;

    switch(backbone_type) {
    case BackboneType::POPOVICS: {
        const auto tmp_a = pow(normal_strain, N) - 1.;
        const auto tmp_b = tmp_a + N;
        trial_stress = peak_stress * normal_strain * N / tmp_b;
        trial_stiffness = peak_stress * N * tmp_a * (1. - N) / peak_strain / tmp_b / tmp_b;
        break;
    }
    case BackboneType::THORENFELDT: {
        const auto tmp_m = normal_strain < 1. ? 1. : M;
        const auto tmp_a = pow(normal_strain, N * tmp_m);
        const auto tmp_b = N - 1. + tmp_a;
        trial_stress = peak_stress * N * normal_strain / tmp_b;
        trial_stiffness = peak_stress / peak_strain * (N * tmp_a * (1. - N * tmp_m) + N * N - N) / tmp_b / tmp_b;
        break;
    }
    case BackboneType::TSAI: {
        const auto tmp_a = pow(normal_strain, N);
        const auto tmp_b = N == 1. ? 1. + (M - 1. + log(normal_strain)) * normal_strain : 1. + (M - N / (N - 1.)) * normal_strain + tmp_a / (N - 1.);
        trial_stress = peak_stress * M * normal_strain / tmp_b;
        trial_stiffness = peak_stress / peak_strain * M * (1. - tmp_a) / tmp_b / tmp_b;
        break;
    }
    case BackboneType::KPS: {
        const auto tmp_a = .5 / ((.29 * peak_stress - 3.) / (145. * peak_stress + 1000.) + peak_strain);
        if(t_strain < peak_strain - .8 / tmp_a) {
            trial_stress = 0.;
            trial_stiffness = 0.;
        } else if(t_strain < peak_strain) {
            trial_stress = confinement * (peak_stress + peak_stress * tmp_a * (t_strain - peak_strain));
            trial_stiffness = confinement * (peak_stress * tmp_a);
        } else {
            trial_stress = confinement * (peak_stress * normal_strain * (2. - normal_strain));
            trial_stiffness = confinement * (2. * peak_stress / peak_strain * (1. - normal_strain));
        }
        break;
    }
    }
}

void Concrete01::compute_tension_backbone(const double t_strain) {
    const auto& residual_strain = trial_history(6); // residual strain

    const auto offset = t_strain - tension_origin;

    if(tension_type == TensionType::EXPONENTIAL) {
        if(offset > crack_strain) {
            trial_stress = crack_stress * pow(crack_strain / offset, .4);
            trial_stiffness = -.4 * trial_stress / offset;
        } else {
            trial_stiffness = crack_stress / (crack_strain - residual_strain);
            trial_stress = trial_stiffness * (offset - residual_strain);
        }
    } else if(tension_type == TensionType::LINEAR) {
        if(offset > ultimate_strain + crack_strain) {
            trial_stress = 0.;
            trial_stiffness = 0.;
        } else if(offset > crack_strain) {
            trial_stiffness = -crack_stress / ultimate_strain;
            trial_stress = crack_stress + trial_stiffness * (offset - crack_strain);
        } else {
            trial_stiffness = crack_stress / (crack_strain - residual_strain);
            trial_stress = trial_stiffness * (offset - residual_strain);
        }
    }
}