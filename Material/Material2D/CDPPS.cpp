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

#include "CDPPS.h"
#include <Recorder/OutputType.h>
#include <Toolbox/tensorToolbox.h>

const double CDPPS::sqrt_three_over_two = sqrt(1.5);
const mat CDPPS::unit_dev_tensor = tensor::unit_deviatoric_tensor4v2();
const uvec CDPPS::P = uvec(std::initializer_list<uword>{ 0, 1, 3 });

podarray<double> CDPPS::compute_backbone(const double f, const double a, const double cb, const double kappa) {
    podarray<double> out(6);

    const auto s_phi = sqrt(1. + a * (a + 2.) * kappa);
    const auto t_phi = (1. + .5 * a) / s_phi;
    const auto b_phi = (1. + a - s_phi) / a;
    const auto p_phi = pow(b_phi, cb);

    out(0) = 1. - p_phi;                                // d
    out(1) = f * s_phi * b_phi;                         // f
    out(2) = out(1) / p_phi;                            // \bar{f}
    out(3) = cb * t_phi * p_phi / b_phi;                // \md{d}
    out(4) = f * t_phi * (1. + a - 2. * s_phi);         // \md{f}
    out(5) = (out(4) + f * t_phi * cb * s_phi) / p_phi; // \md{\bar{f}}

    return out;
}

double CDPPS::r(const vec& in) {
    const auto abs_sum = accu(abs(in));
    const auto r = abs_sum > datum::eps ? .5 + .5 * accu(in) / abs_sum : 0.;
    if(r > 1.) return 1.;
    if(r < 0.) return 0.;
    return r;
}

rowvec CDPPS::dr(const vec& in) {
    const auto abs_sum = accu(abs(in));

    vec out(3);

    out.fill(abs_sum);

    if(abs_sum <= datum::eps) return out.t();

    out -= accu(in) * sign(in);
    out /= 2. * abs_sum * abs_sum;

    return out.t();
}

double CDPPS::s(const double r) const { return s0 + (1. - s0) * r; }

int CDPPS::compute_regular_return(const vec& t_strain) {
    trial_full_strain = t_strain;
    trial_history = current_history;
    trial_plastic_strain = current_plastic_strain;

    const auto& kappa_t = trial_history(0);
    const auto& kappa_c = trial_history(1);

    // predictor \sigma^{tr}
    trial_full_stress = initial_full_stiffness * (trial_full_strain - trial_plastic_strain);

    // principal predictor \hat{\sigma}^{tr}
    vec p_predictor;
    mat trans_mat;
    if(!eig_sym(p_predictor, trans_mat, tensor::stress::to_tensor(trial_full_stress), "std")) return -1;

    // deviatoric principal predictor \hat{s}^{tr}
    const auto d_predictor = tensor::dev(p_predictor);
    // unit deviatoric principal predictor \hat{n}
    const vec u_predictor = normalise(d_predictor);

    // some norm related vectors
    const rowvec c_pfpsigma = sqrt_three_over_two * u_predictor.t() + alpha;
    const vec dsigmadlambda = -double_shear * u_predictor - factor_a;
    const auto dgdsigma_t = (u_predictor(2) + alpha_p) / g_t;
    const auto dgdsigma_c = (u_predictor(0) + alpha_p) / g_c;

    // constant part of yield function throughout the whole iteration loop
    const auto const_yield = alpha * tensor::stress::invariant1(p_predictor) + sqrt(3. * tensor::stress::invariant2(d_predictor));

    t_para = compute_backbone(f_t, a_t, cb_t, kappa_t);
    c_para = compute_backbone(f_c, a_c, cb_c, kappa_c);

    auto tension_flag = p_predictor(2) > datum::eps;

    // update beta
    auto beta = tension_flag ? -c_para(2) / t_para(2) * factor_c - 1. - alpha : 0.;

    auto t_yield = const_yield + factor_c * c_para(2);

    r_weight = r(p_predictor);

    // elastic loading/unloading
    if((tension_flag ? t_yield + beta * p_predictor(2) : t_yield) < tolerance) {
        elastic_flag = true;
        trial_stiffness = initial_full_stiffness;
        return 0;
    }

    vec e_stress, h(2), phpkappa(2), incre_history;
    rowvec pfpsigma = c_pfpsigma, pfpkappa(2);
    mat pqpkappa, pqpsigma;

    if(tension_flag) pfpsigma(2) = c_pfpsigma(2) + beta;

    auto i_lambda = 0.;

    auto counter = 0;
    while(++counter < 20) {
        // store previoud lambda for convergence test
        const auto p_lambda = i_lambda;

        // compute new increment of lambda
        i_lambda = beta == 0. ? t_yield / factor_b : (t_yield + beta * p_predictor(2)) / (factor_b - beta * dsigmadlambda(2));
        const auto abs_error = abs(i_lambda - p_lambda);
        suanpan_extra_debug("CDPPS local iterative loop error: %.5E.\n", abs_error);
        if(counter != 1 && abs_error < tolerance && norm(incre_history) < tolerance) break;

        // update principal effective stress
        e_stress = p_predictor + i_lambda * dsigmadlambda;

        tension_flag = e_stress(2) > datum::eps;

        if(tension_flag) {
            const auto t_factor = factor_c * e_stress(2) / t_para(2);
            pfpkappa(0) = c_para(2) / t_para(2) * t_factor * t_para(5);
            pfpkappa(1) = (factor_c - t_factor) * c_para(5);
        } else {
            pfpkappa(0) = 0.;
            pfpkappa(1) = factor_c * c_para(5);
        }

        r_weight = r(e_stress);

        h.zeros(), phpkappa.zeros();
        if(r_weight != 0.) {
            phpkappa(0) = h(0) = r_weight * dgdsigma_t;
            h(0) *= t_para(1), phpkappa(0) *= t_para(4);
        }
        if(r_weight != 1.) {
            phpkappa(1) = h(1) = (1. - r_weight) * dgdsigma_c;
            h(1) *= c_para(1), phpkappa(1) *= c_para(4);
        }

        pqpkappa = diagmat(i_lambda * phpkappa - 1.);

        drweight = dr(e_stress);

        pqpsigma = i_lambda * vec{ t_para(1) * dgdsigma_t, -c_para(1) * dgdsigma_c } * drweight;

        // compute the increment of damage parameters
        incre_history = solve(mat(pqpkappa - (pqpsigma * dsigmadlambda + h) / dot(pfpsigma, dsigmadlambda) * pfpkappa), trial_history - current_history - i_lambda * h);
        trial_history += incre_history;

        // now update damage parameters
        t_para = compute_backbone(f_t, a_t, cb_t, kappa_t);
        c_para = compute_backbone(f_c, a_c, cb_c, kappa_c);

        // update portion of yield function
        t_yield = const_yield + factor_c * c_para(2);

        // update beta if trial status is in tension
        beta = tension_flag ? -c_para(2) / t_para(2) * factor_c - 1. - alpha : 0.;
        pfpsigma(2) = tension_flag ? c_pfpsigma(2) + beta : c_pfpsigma(2);
    }

    suanpan_debug("CDPPS state determination loop counter: %u.\n", counter);

    if(counter == 20) {
        suanpan_error("CDPPS cannot converge within 20 iterations.\n");
        return -1;
    }

    // effective stress
    trial_full_stress = transform::compute_jacobian_principal_to_nominal(trans_mat) * e_stress;

    const auto d_stress = tensor::dev(trial_full_stress);
    const auto n_norm = tensor::stress::norm(d_stress);
    const vec u_stress = d_stress / n_norm;
    const vec dgdsigma = u_stress + unit_alpha_p;

    trial_plastic_strain += i_lambda * dgdsigma;

    const mat hessian = inv(mat(inv_stiffness + i_lambda / n_norm * (unit_dev_tensor - u_stress * u_stress.t())));
    const vec hessian_dgdsigma = hessian * dgdsigma;

    qkappah = solve(pqpkappa, h);
    qkappaqsigma = solve(pqpkappa, pqpsigma);

    jacobian = transform::compute_jacobian_nominal_to_principal(trans_mat);

    const rowvec dfdsigma = (pfpsigma - pfpkappa * qkappaqsigma) * jacobian;

    dlambdadepsilon = dfdsigma * hessian / (dot(dfdsigma, hessian_dgdsigma) + dot(pfpkappa, qkappah));

    trial_stiffness = hessian - hessian_dgdsigma * dlambdadepsilon;

    return 0;
}

mat CDPPS::form_stiffness(const mat& full_stiffness) {
    mat t_stiffness = full_stiffness(P, P);

    if(abs(full_stiffness(2, 2)) < datum::eps) return t_stiffness;

    for(auto I = 0; I < 3; ++I)
        for(auto J = 0; J < 3; ++J) t_stiffness(I, J) -= full_stiffness(P[I], 2) * full_stiffness(2, P[J]) / full_stiffness(2, 2);

    return t_stiffness;
}

CDPPS::CDPPS(const unsigned T, const double E, const double V, const double ST, const double SC, const double GT, const double GC, const double AT, const double AC, const double DT, const double DC, const double AP, const double BC, const double S, const double R)
    : Material3D(T, MT_CDPPS, R)
    , elastic_modulus(E)
    , poissons_ratio(V < .5 ? V : .2)
    , double_shear(elastic_modulus / (1. + poissons_ratio))
    , a_t(AT < 1. ? AT : .5)
    , cb_t(log(DT < 1. ? DT : .8) / log(.5 * (1. + a_t - sqrt(1. + a_t * a_t)) / a_t))
    , g_t(GT > 0. ? GT : -GT)
    , f_t(ST > 0. ? ST : -ST)
    , a_c(AC > 1. ? AC : 4.)
    , cb_c(log(DC < 1. ? DC : .8) / log(.5 + .5 / a_c))
    , g_c(GC > 0. ? GC : -GC)
    , f_c((SC < 0. ? SC : -SC) * 4. * a_c / pow(1. + a_c, 2.))
    , alpha((BC - 1.) / (2. * BC - 1.))
    , alpha_p(AP)
    , s0(S)
    , factor_a(elastic_modulus / (1. - 2. * poissons_ratio) * alpha_p)
    , factor_b(3. * alpha * factor_a + sqrt_three_over_two * double_shear)
    , factor_c(1. - alpha)
    , tolerance(1E-12)
    , unit_alpha_p(alpha_p * tensor::unit_tensor2())
    , elastic_flag(false)
    , r_weight(0.) {
    // tension
    const auto half_stress = .5 * f_t;
    const auto half_strain = log(1. + a_t + sqrt(1. + a_t * a_t)) / f_t / (1. + .5 * a_t) * g_t + half_stress / elastic_modulus;
    const auto ratio_t = half_stress / half_strain / elastic_modulus;
    if(ratio_t >= DT) {
        suanpan_warning("CDPPS model requires a minimum tension degradation of %.2f, now reset it.\n", ratio_t);
        const auto new_ratio = ceil(100. * ratio_t + 5.) * .01;
        access::rw(cb_t) = log(new_ratio > 1. ? .9 : new_ratio) / log(.5 * (1. + a_t - sqrt(1. + a_t * a_t)) / a_t);
    }
    // compression
    const auto peak_stress = .25 * f_c * pow(1. + a_c, 2.) / a_c;
    const auto peak_strain = log(2. * a_c / (1. + a_c)) / f_c / (1. + .5 * a_c) * g_c + peak_stress / elastic_modulus;
    const auto ratio_c = peak_stress / peak_strain / elastic_modulus;
    if(ratio_c >= DC) {
        suanpan_warning("CDPPS model requires a minimum compression degradation of %.2f, now reset it.\n", ratio_c);
        const auto new_ratio = ceil(100. * ratio_c + 5.) * .01;
        access::rw(cb_c) = log(new_ratio > 1. ? .9 : new_ratio) / log(.5 + .5 / a_c);
    }
}

void CDPPS::initialize(const shared_ptr<DomainBase>&) {
    inv_stiffness.zeros(6, 6);

    vec F(5);
    F(0) = 1. / elastic_modulus;
    F(1) = -poissons_ratio / elastic_modulus;
    F(2) = (2. * poissons_ratio + 2.) / elastic_modulus;
    F(3) = poissons_ratio / (.5 - poissons_ratio) / F(2);
    F(4) = elastic_modulus / (1. - poissons_ratio * poissons_ratio);

    for(auto I = 0; I < 3; ++I)
        for(auto J = 0; J < 3; ++J) inv_stiffness(I, J) = I == J ? F(0) : F(1);
    for(auto I = 3; I < 6; ++I) inv_stiffness(I, I) = F(2);

    initial_stiffness.zeros(3, 3);
    initial_stiffness(0, 1) = initial_stiffness(1, 0) = poissons_ratio * (initial_stiffness(0, 0) = initial_stiffness(1, 1) = F(4));
    initial_stiffness(2, 2) = .5 * F(4) * (1. - poissons_ratio);

    trial_stiffness = current_stiffness = initial_stiffness;

    initial_full_stiffness.zeros(6, 6);
    for(auto I = 0; I < 3; ++I)
        for(auto J = 0; J < 3; ++J) initial_full_stiffness(I, J) = F(3);
    for(auto I = 0; I < 3; ++I) initial_full_stiffness(I, I) += double_shear;
    for(auto I = 3; I < 6; ++I) initial_full_stiffness(I, I) = .5 * double_shear;

    trial_history.zeros(2);
    current_history.zeros(2);
    trial_plastic_strain.zeros(6);
    current_plastic_strain.zeros(6);
    trial_full_strain.zeros(6);
    current_full_strain.zeros(6);
    trial_full_stress.zeros(6);
    current_full_stress.zeros(6);
}

unique_ptr<Material> CDPPS::get_copy() { return make_unique<CDPPS>(*this); }

double CDPPS::get_parameter(const ParameterType& PT) const {
    switch(PT) {
    case ParameterType::DENSITY:
        return density;
    case ParameterType::ELASTICMODULUS:
    case ParameterType::YOUNGSMODULUS:
        return elastic_modulus;
    case ParameterType::POISSONSRATIO:
        return poissons_ratio;
    default:
        return 0.;
    }
}

int CDPPS::update_trial_status(const vec& t_strain) {
    elastic_flag = false;
    trial_strain = t_strain;
    trial_full_strain(P) = trial_strain;

    const auto& stress_zz = trial_full_stress.at(2);
    const auto& stress_yz = trial_full_stress.at(4);
    const auto& stress_zx = trial_full_stress.at(5);

    auto counter = 0;
    auto error = 0.;
    while(++counter < 50) {
        if(compute_regular_return(trial_full_strain) != 0) return -1;
        trial_full_strain(2) -= stress_zz / trial_stiffness(2, 2);
        trial_full_strain(4) -= stress_yz / trial_stiffness(4, 4);
        trial_full_strain(5) -= stress_zx / trial_stiffness(5, 5);
        error = fabs(stress_zz) + fabs(stress_yz) + fabs(stress_zx);
        if(error < 1E-10) break;
    }

    if(counter == 20) suanpan_extra_debug("CDPPS cannot converge in 20 iterations, current error is %.3E.\n", error);

    const auto s_weight = s(r_weight);
    const auto d_a = 1. - c_para(0);
    const auto d_b = 1. - s_weight * t_para(0);

    const auto damage_factor = d_a * d_b;

    if(!elastic_flag) {
        const auto pdpkappa = rowvec{ s_weight * d_a * t_para(3), d_b * c_para(3) };

        mat part_stifness = trial_full_stress * (pdpkappa * qkappaqsigma + t_para(0) * (c_para(0) - 1.) * (1. - s0) * drweight) * jacobian;
        part_stifness.diag() += damage_factor;
        trial_stiffness = part_stifness * trial_stiffness;
        trial_stiffness += dot(pdpkappa, qkappah) * trial_full_stress * dlambdadepsilon;
    } else
        trial_stiffness *= damage_factor;

    trial_stress = trial_full_stress(P);
    trial_stiffness = form_stiffness(trial_stiffness);

    trial_stress *= damage_factor;

    return 0;
}

int CDPPS::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_full_strain.zeros();
    current_full_stress.zeros();
    current_history.zeros();
    current_plastic_strain.zeros();
    current_stiffness = initial_stiffness;
    return reset_status();
}

int CDPPS::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_full_strain = trial_full_strain;
    current_full_stress = trial_full_stress;
    current_history = trial_history;
    current_plastic_strain = trial_plastic_strain;
    current_stiffness = trial_stiffness;
    return 0;
}

int CDPPS::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_full_strain = current_full_strain;
    trial_full_stress = current_full_stress;
    trial_history = current_history;
    trial_plastic_strain = current_plastic_strain;
    trial_stiffness = current_stiffness;
    return 0;
}

vector<vec> CDPPS::record(const OutputType& T) {
    vector<vec> data;
    switch(T) {
    case OutputType::DT:
        data.emplace_back(vec{ t_para(0) });
        break;
    case OutputType::DC:
        data.emplace_back(vec{ c_para(0) });
        break;
    case OutputType::KAPPAT:
        data.emplace_back(vec{ trial_history(0) });
        break;
    case OutputType::KAPPAC:
        data.emplace_back(vec{ trial_history(1) });
        break;
    case OutputType::PE:
        data.emplace_back(trial_plastic_strain);
        break;
    case OutputType::PE11:
        data.emplace_back(vec{ trial_plastic_strain(0) });
        break;
    case OutputType::PE22:
        data.emplace_back(vec{ trial_plastic_strain(1) });
        break;
    case OutputType::PE33:
        data.emplace_back(vec{ trial_plastic_strain(2) });
        break;
    case OutputType::E:
        data.emplace_back(trial_strain);
        break;
    case OutputType::E11:
        data.emplace_back(vec{ trial_strain(0) });
        break;
    case OutputType::E22:
        data.emplace_back(vec{ trial_strain(1) });
        break;
    case OutputType::E33:
        data.emplace_back(vec{ trial_strain(2) });
        break;
    case OutputType::S:
        data.emplace_back(trial_stress);
        break;
    case OutputType::S11:
        data.emplace_back(vec{ trial_stress(0) });
        break;
    case OutputType::S22:
        data.emplace_back(vec{ trial_stress(1) });
        break;
    case OutputType::S33:
        data.emplace_back(vec{ trial_stress(2) });
        break;
    default:
        break;
    }
    return data;
}
