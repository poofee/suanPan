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

#include "CDPM2.h"
#include <Recorder/OutputType.h>
#include <Toolbox/tensorToolbox.h>
#include <Toolbox/utility.h>

const double CDPM2::sqrt_six = sqrt(6.);
const mat CDPM2::unit_dev_tensor = tensor::unit_deviatoric_tensor4v2();
const rowvec CDPM2::unit_tensor2 = tensor::unit_tensor2().t();

CDPM2::CDPM2(const unsigned T, const double E, const double V, const double R)
    : Material3D(T, MT_CDPM2, R)
    , elastic_modulus(E)
    , poissons_ratio(V)
    , double_shear(elastic_modulus / (1. + poissons_ratio))
    , bulk_modulus(elastic_modulus / (3. - 6. * poissons_ratio))
    , dsigmadepsilon(bulk_modulus * unit_tensor2) {}

void CDPM2::initialize(const shared_ptr<DomainBase>&) {
    inv_stiffness.zeros(6, 6);

    vec F(4);
    F(0) = 1. / elastic_modulus;
    F(1) = -poissons_ratio / elastic_modulus;
    F(2) = (2. * poissons_ratio + 2.) / elastic_modulus;
    F(3) = poissons_ratio / (.5 - poissons_ratio) / F(2);

    for(auto I = 0; I < 3; ++I)
        for(auto J = 0; J < 3; ++J) inv_stiffness(I, J) = I == J ? F(0) : F(1);
    for(auto I = 3; I < 6; ++I) inv_stiffness(I, I) = F(2);

    initial_stiffness.zeros(6, 6);

    for(auto I = 0; I < 3; ++I)
        for(auto J = 0; J < 3; ++J) initial_stiffness(I, J) = F(3);

    for(auto I = 0; I < 3; ++I) initial_stiffness(I, I) += double_shear;
    for(auto I = 3; I < 6; ++I) initial_stiffness(I, I) = .5 * double_shear;

    trial_stiffness = current_stiffness = initial_stiffness;

    trial_history.zeros(20);
    current_history.zeros(20);
    trial_plastic_strain.zeros(6);
    current_plastic_strain.zeros(6);
}

unique_ptr<Material> CDPM2::get_copy() { return make_unique<CDPM2>(*this); }

double CDPM2::get_parameter(const ParameterType&) const { return 0.; }

int CDPM2::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;
    trial_history = current_history;

    vec e_strain = trial_strain - current_plastic_strain;
    trial_stress = initial_stiffness * e_strain;

    vec p_stress;
    mat p_axis;
    if(!eig_sym(p_stress, p_axis, tensor::stress::to_tensor(trial_stress), "std")) return SUANPAN_FAIL;

    const auto trans = transform::compute_jacobian_principal_to_nominal(p_axis);

    const auto sigma = mean(p_stress);
    const vec s_stress = p_stress - sigma;
    const auto rho = norm(s_stress);
    const auto theta = tensor::stress::lode(s_stress);

    auto& kappap = trial_history(0);

    if(compute_yield_function(sigma, rho, theta, kappap) < tolerance) {
        trial_stiffness = initial_stiffness;
        return SUANPAN_SUCCESS;
    }

    // auto new_sigma = sigma, new_kappap = kappap;

    // auto flag = false;
    // if the trial point is in pure tension region and close to the hydrostatic axis
    // if(sigma > 0. && rho < 1E-8) flag = compute_vertex_return(new_sigma, rho, new_kappap);

    mat jacobian;
    vec solution;

    // if(!flag) {
    auto status = compute_plasticity(sigma, rho, theta, kappap);

    if(!compute_regular_return(solution, jacobian, status, vec(std::initializer_list<double>{ sigma, rho, kappap, 0. }), theta)) return -1;

    p_stress = transform::haigh_westergaard_to_principal(solution(0), solution(1), theta);

    auto damage = compute_damage(solution(0), solution(1), theta);

    vec damage_vec(3);
    for(uword I = 0; I < damage_vec.n_elem; ++I) damage_vec(I) = 1. - p_stress(I) > 0. ? damage(0) : damage(1);

    // update kappap
    kappap = solution(2);

    //} else {
    //    for(uword I = 0; I < 3; ++I) p_stress(I) = new_sigma;
    //    damage = compute_damage(new_sigma, 0., theta);
    //    kappap = new_kappap;
    //}

    // update plastic strain
    trial_plastic_strain = trial_strain - inv_stiffness * trans * p_stress;
    const vec incre_plastic_strain = trial_plastic_strain - current_plastic_strain;

    // compute stiffness
    const auto& trial_sigma = solution(0);
    const auto& trial_rho = solution(1);
    const auto& trial_lambda = solution(3);

    const auto& pgpsigma = status(4);
    const auto& pgprho = status(5);
    const auto& ppgppsigma = status(6);
    const auto& ppgpprho = status(7);
    const auto& ppgpsigmaprho = status(8);
    const auto& ppgpsigmapkappap = status(9);
    const auto& ppgprhopkappap = status(10);
    const auto& pfptheta = status(21);

    /** \dfrac{\mathrm{d}\theta^{tr}}{\mathrm{d}\mathbold{e}} **/

    // deviatoric trial elastic strain
    e_strain(span(0, 2)) -= mean(e_strain(span(0, 2)));
    // invariant
    const auto invariant2 = tensor::strain::invariant2(e_strain);
    // deviatoric trial elastic strain tensor
    const auto de_strain = tensor::strain::to_tensor(e_strain);
    // according to the definition of the lode angle using strain space
    rowvec dthetadepsilon(6, fill::zeros);
    const auto sin_theta = sin(3. * theta);
    if(abs(sin_theta) > datum::eps) dthetadepsilon = sqrt(.75) * pow(invariant2, -2.5) / sin_theta / double_shear * (1.5 * tensor::strain::invariant3(e_strain) * e_strain.t() - invariant2 * tensor::strain::to_voigt(de_strain * de_strain).t()) * unit_dev_tensor;

    mat history(4, 6);
    history.row(0) = dsigmadepsilon;
    history.row(1) = double_shear / rho * tensor::dev(trial_stress).t() * unit_dev_tensor;
    history.row(2).zeros();
    history.row(3) = -pfptheta * dthetadepsilon;

    // unit deviatoric stress
    const vec u_stress = (p_stress - trial_sigma) / trial_rho;

    mat future(3, 4);
    future.col(0) = -trial_lambda * (ppgpsigmaprho * u_stress + ppgppsigma / 3.);
    future.col(1) = -trial_lambda * (ppgpprho * u_stress + ppgpsigmaprho / 3.);
    future.col(2) = -trial_lambda * (ppgprhopkappap * u_stress + ppgpsigmapkappap / 3.);
    future.col(3) = -pgprho * u_stress - pgpsigma / 3.;

    const mat time_machine = solve(jacobian, history);

    mat t_stiffness = trans * future * time_machine;
    t_stiffness.diag() += 1.;
    // trial_stiffness = diagmat(damage_vec) * transform::compute_jacobian_nominal_to_principal(p_axis) * initial_stiffness * t_stiffness;

    // mat dddepsilon(3, 6, fill::zeros);
    // trial_stiffness += diagmat(p_stress) * dddepsilon;

    // trial_stiffness = trans * trial_stiffness;

    trial_stiffness = initial_stiffness * t_stiffness;

    trial_stress = trans * damage_vec % p_stress;

    suanpan_debug([&]() {
        if(!trial_stress.is_finite() || !trial_stiffness.is_finite()) throw;
    });

    return 0;
}

int CDPM2::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history.zeros();
    current_plastic_strain.zeros();
    current_stiffness = initial_stiffness;
    return reset_status();
}

int CDPM2::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_plastic_strain = trial_plastic_strain;
    current_stiffness = trial_stiffness;
    return 0;
}

int CDPM2::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_plastic_strain = current_plastic_strain;
    trial_stiffness = current_stiffness;
    return 0;
}

vector<vec> CDPM2::record(const OutputType& P) {
    vector<vec> data;
    switch(P) {
    case OutputType::KAPPAP:
        data.emplace_back(vec{ trial_history(0) });
        break;
    case OutputType::DT:
        data.emplace_back(vec{ 0. });
        break;
    case OutputType::DC:
        data.emplace_back(vec{ 0. });
        break;
    default:
        break;
    }
    return data;
}

podarray<double> CDPM2::compute_plasticity(const double sigma, const double rho, const double theta, const double kappap) const {
    podarray<double> status;
    status.zeros(30);

    auto& f = status(0);
    auto& pfpsigma = status(1);
    auto& pfprho = status(2);
    auto& pfpkappap = status(3);
    auto& pgpsigma = status(4);
    auto& pgprho = status(5);
    auto& ppgppsigma = status(6);
    auto& ppgpprho = status(7);
    auto& ppgpsigmaprho = status(8);
    auto& ppgpsigmapkappap = status(9);
    auto& ppgprhopkappap = status(10);
    auto& pkappapplambda = status(11);
    auto& ppkappapplambdapsigma = status(12);
    auto& ppkappapplambdaprho = status(13);
    auto& ppkappapplambdapkappap = status(14);
    auto& qh1 = status(15);
    auto& dqh1 = status(16);
    auto& qh2 = status(17);
    auto& dqh2 = status(18);
    auto& r = status(19);
    auto& norm_m = status(20);
    auto& pfptheta = status(21);

    const auto cos_theta = cos(theta);
    const auto cos_theta2 = 4. * cos_theta * cos_theta;
    const auto cos_theta3 = cos_theta2 * ea;
    const auto rc = sqrt(cos_theta3 + ed);
    const auto rb = 2. * ea * cos_theta + eb * rc;

    // deviatoric section
    r = (cos_theta3 + ec) / rb; // eq. 19
    const auto dr = (r + (eb * r / rc - 2.) * 2. * cos_theta) * 2. * ea / rb * sin(theta);

    const auto rh = -sigma / fc - 1. / 3.;                                                        // eq. 34
    const auto kaih = rh >= datum::eps ? ah + (bh - ah) * exp(-rh / ch) : eh * exp(rh / fh) + dh; // eq. 33
    const auto dkaih = exp(rh / (rh >= datum::eps ? -ch : fh)) * gh;

    const auto quad_kaih = cos_theta2 / kaih;

    double pmpsigma, pmprho, pmpkappap;

    if(kappap < 1.) {
        if(kappap >= 0.) {
            const auto kappap2 = kappap * kappap;
            const auto kappap3 = kappap * kappap2;
            const auto kappapa = kappap3 - 3. * kappap2 + 2. * kappap;
            const auto kappapb = 3. * kappap2 - 6. * kappap + 2.;

            qh1 = qh0 + one_minus_qh0 * (kappapa + kappap) - hp * kappapa; // eq. 30
            if(qh1 > 1.) qh1 = 1.;

            dqh1 = one_minus_qh0 * (kappapb + 1.) - hp * kappapb;
        } else
            qh1 = qh0;

        qh2 = 1.;
        // dqh2 = 0.;

        const auto one_minus_qh1 = 1. - qh1;
        const auto qh1_qh1 = qh1 * qh1;

        const auto ag = 3. * ft_over_fc + half_m0; // eq. 28
        const auto dg = log(ag / (3. + half_m0)) + log_df;
        const auto bg = (1. + ft_over_fc) / 3. / dg; // eq. 29
        const auto eg = (sigma - ft / 3.) / bg / fc;
        const auto cg = exp(eg); // eq. 23

        // const auto mg = fc * ag * bg * cg;

        const auto pmgpsigma = ag * cg;
        // const auto pmgpkappap = 0.;

        const auto af = rho / sqrt_six_fc;
        const auto bf = sigma / fc;
        const auto cf = af + bf;
        const auto ef = one_minus_qh1 * cf * cf + 3. * af;
        const auto ef2 = 2. * ef;
        const auto cf2 = 2. * cf;

        const auto pefpsigma = one_minus_qh1 * cf2 / fc;
        const auto pefprho = (one_minus_qh1 * cf2 + 3.) / sqrt_six_fc;
        const auto pefpkappap = -cf * cf * dqh1;
        const auto ppefppsigma = one_minus_qh1 * 2. / fc_fc;
        // const auto ppefpprho = one_minus_qh1 / 3. / fc_fc;
        const auto ppefpprho = ppefppsigma / 6.;
        // const auto ppefpsigmaprho = one_minus_qh1 / sqrt_three_over_two / fc_fc;
        const auto ppefpsigmaprho = ppefppsigma / sqrt_six;
        const auto ppefpsigmapkappap = -cf2 / fc * dqh1;
        // const auto ppefprhopkappap = -cf2 / sqrt_six_fc * dqh1;
        const auto ppefprhopkappap = ppefpsigmapkappap / sqrt_six;

        f = ef * ef + (m0 * (r * af + bf) - 1.) * qh1_qh1; // eq. 18

        pfpsigma = ef2 * pefpsigma + qh1_qh1 * m0 / fc;
        pfprho = ef2 * pefprho + qh1_qh1 * m0 * r / sqrt_six_fc;
        pfpkappap = ef2 * pefpkappap + 2. * qh1 * dqh1 * (m0 * (r * af + bf) - 1.);
        pfptheta = m0 * qh1_qh1 * af * dr;

        // const auto g = ef * ef + qh1_qh1 * (m0 * af + mg / fc);

        pgpsigma = ef2 * pefpsigma + qh1_qh1 * pmgpsigma / fc;
        pgprho = ef2 * pefprho + qh1_qh1 * m0 / sqrt_six_fc;

        ppgppsigma = 2. * (pefpsigma * pefpsigma + ef * ppefppsigma) + qh1_qh1 * pmgpsigma / bg / fc_fc;
        ppgpprho = 2. * (pefprho * pefprho + ef * ppefpprho);
        ppgpsigmaprho = 2. * (pefprho * pefpsigma + ef * ppefpsigmaprho);

        ppgpsigmapkappap = 2. * (pefpkappap * pefpsigma + ef * ppefpsigmapkappap + qh1 * dqh1 * pmgpsigma / fc);
        ppgprhopkappap = 2. * (pefpkappap * pefprho + ef * ppefprhopkappap + qh1 * dqh1 * m0 / sqrt_six_fc);

        norm_m = sqrt(pgpsigma * pgpsigma / 3. + pgprho * pgprho);

        pmpsigma = (pgpsigma * ppgppsigma / 3. + pgprho * ppgpsigmaprho) / norm_m;
        pmprho = (pgpsigma * ppgpsigmaprho / 3. + pgprho * ppgpprho) / norm_m;
        pmpkappap = (pgpsigma * ppgpsigmapkappap / 3. + pgprho * ppgprhopkappap) / norm_m;
    } else {
        qh1 = 1.;
        // dqh1 = 0.;

        qh2 = 1. + hp * (kappap - 1.);

        dqh2 = hp;

        const auto ag = 3. * ft_over_fc * qh2 + half_m0;
        const auto dg = log(ag / (3. * qh2 + half_m0)) + log_df;
        const auto bg = qh2 * (1. + ft_over_fc) / 3. / dg;
        const auto eg = (sigma - qh2 * ft / 3.) / bg / fc;
        const auto cg = exp(eg);
        // const auto mg = fc * ag * bg * cg;

        const auto dagdkappap = 3. * ft_over_fc * dqh2;
        const auto ddgdkappap = dagdkappap / ag - 3. * dqh2 / (3. * qh2 + half_m0);
        const auto dbgdkappap = (1. + ft_over_fc) / 3. * (dqh2 * dg - qh2 * ddgdkappap) / dg / dg;
        const auto pegpkappap = (-dqh2 * ft / 3. - (sigma - qh2 * ft / 3.) * dbgdkappap / bg) / bg / fc;
        const auto pcgpsigma = cg / bg / fc;
        const auto pcgpkappap = cg * pegpkappap;

        const auto pmgpsigma = ag * cg;
        // const auto pmgpkappap = fc * (dagdkappap * bg * cg + ag * dbgdkappap * cg + ag * bg * pcgpkappap);

        const auto af = rho / sqrt_six_fc;
        const auto bf = sigma / fc;
        const auto cf = r * af + bf;

        f = 9. * af * af + m0 * qh2 * cf - qh2 * qh2;

        pfpsigma = m0 * qh2 / fc;
        pfprho = (18. * af + m0 * qh2 * r) / sqrt_six_fc;
        pfpkappap = (m0 * cf - 2. * qh2) * dqh2;
        pfptheta = m0 * qh2 * af * dr;

        // const auto g = 1.5 * rho * rho / fc2 + m0 * rho / sqrt_six_fc + mg / fc;
        pgpsigma = pmgpsigma / fc;
        pgprho = 3. * rho / fc_fc + m0 / sqrt_six_fc;

        ppgppsigma = ag / fc * pcgpsigma;
        ppgpprho = 3. / fc_fc;

        ppgpsigmapkappap = (dagdkappap * cg + ag * pcgpkappap) / fc;

        norm_m = sqrt(pgpsigma * pgpsigma / 3. + pgprho * pgprho);

        pmpsigma = pgpsigma * ppgppsigma / norm_m / 3.;
        pmprho = pgprho * ppgpprho / norm_m;
        pmpkappap = pgpsigma * ppgpsigmapkappap / norm_m / 3.;
    }

    pkappapplambda = quad_kaih * norm_m;

    ppkappapplambdapsigma = quad_kaih * (pmpsigma - norm_m * dkaih / kaih);
    ppkappapplambdaprho = quad_kaih * pmprho;
    ppkappapplambdapkappap = quad_kaih * pmpkappap;

    return status;
}

mat CDPM2::compute_jacobian(const podarray<double>& status, const double lambda) const {
    mat jacobian(4, 4);

    const auto& pfpsigma = status(1);
    const auto& pfprho = status(2);
    const auto& pfpkappap = status(3);
    const auto& pgpsigma = status(4);
    const auto& pgprho = status(5);
    const auto& ppgppsigma = status(6);
    const auto& ppgpprho = status(7);
    const auto& ppgpsigmaprho = status(8);
    const auto& ppgpsigmapkappap = status(9);
    const auto& ppgprhopkappap = status(10);
    const auto& pkappapplambda = status(11);
    const auto& ppkappapplambdapsigma = status(12);
    const auto& ppkappapplambdaprho = status(13);
    const auto& ppkappapplambdapkappap = status(14);

    const auto bulk_lambda = lambda * bulk_modulus;
    const auto shear_lambda = lambda * double_shear;

    jacobian(0, 0) = bulk_lambda * ppgppsigma + 1.;
    jacobian(0, 1) = bulk_lambda * ppgpsigmaprho;
    jacobian(0, 2) = bulk_lambda * ppgpsigmapkappap;
    jacobian(0, 3) = bulk_modulus * pgpsigma;

    jacobian(1, 0) = shear_lambda * ppgpsigmaprho;
    jacobian(1, 1) = shear_lambda * ppgpprho + 1.;
    jacobian(1, 2) = shear_lambda * ppgprhopkappap;
    jacobian(1, 3) = double_shear * pgprho;

    jacobian(2, 0) = -lambda * ppkappapplambdapsigma;
    jacobian(2, 1) = -lambda * ppkappapplambdaprho;
    jacobian(2, 2) = -lambda * ppkappapplambdapkappap + 1.;
    jacobian(2, 3) = -pkappapplambda;

    jacobian(3, 0) = pfpsigma;
    jacobian(3, 1) = pfprho;
    jacobian(3, 2) = pfpkappap;
    jacobian(3, 3) = 0.;

    return jacobian;
}

vec CDPM2::compute_residual(const podarray<double>& status, const vec& solution, const vec& target) const {
    const auto& f = status(0);
    const auto& pgpsigma = status(4);
    const auto& pgprho = status(5);
    const auto& pkappapplambda = status(11);

    const auto& sigma = solution(0);
    const auto& rho = solution(1);
    const auto& kappap = solution(2);
    const auto& lambda = solution(3);

    vec estimation(4);

    estimation(0) = sigma + lambda * bulk_modulus * pgpsigma;
    estimation(1) = rho + lambda * double_shear * pgprho;
    estimation(2) = kappap - lambda * pkappapplambda;
    estimation(3) = f;

    return estimation - target;
}

/**
 * \brief Simple version of computing the yield function.
 * \param sigma
 * \param rho
 * \param theta
 * \param kappap
 * \return yield_function
 */
double CDPM2::compute_yield_function(const double sigma, const double rho, const double theta, const double kappap) const {
    double f;

    const auto cos_theta = cos(theta);
    const auto cos_theta2 = 4. * cos_theta * cos_theta;
    const auto cos_theta3 = cos_theta2 * ea;

    const auto r = (cos_theta3 + ec) / (2. * ea * cos_theta + eb * sqrt(cos_theta3 + ed)); // eq. 19

    if(kappap < 1.) {
        double qh1;
        if(kappap >= 0.) {
            const auto kappap2 = kappap * kappap;
            const auto kappapa = kappap * kappap2 - 3. * kappap2 + 2. * kappap;
            qh1 = qh0 + one_minus_qh0 * (kappapa + kappap) - hp * kappapa; // eq. 30
            if(qh1 > 1.) qh1 = 1.;
        } else
            qh1 = qh0;

        const auto af = rho / sqrt_six_fc;
        const auto bf = sigma / fc;
        const auto cf = af + bf;
        const auto ef = (1. - qh1) * cf * cf + 3. * af;

        f = ef * ef + (m0 * (r * af + bf) - 1.) * qh1 * qh1; // eq. 18
    } else {
        const auto qh2 = 1. + hp * (kappap - 1.);

        const auto af = rho / sqrt_six_fc;

        f = 9. * af * af + m0 * qh2 * (r * af + sigma / fc) - qh2 * qh2;
    }

    return f;
}

/**
 * \brief Compute the increment of kappap when returning from the given point to another point on the hydrostatic axis.
 * \param sigma
 * \param rho
 * \param new_sigma
 * \return incre_kappap
 */
double CDPM2::compute_incre_kappap(const double sigma, const double rho, const double new_sigma) const { return sqrt(pow((new_sigma - sigma) / 3. / bulk_modulus, 2.) + pow(rho / double_shear, 2.)) / compute_kaih(new_sigma); }

double CDPM2::compute_kaih(const double sigma) const {
    const auto rh = -sigma / fc - 1. / 3.;                                             // eq. 34
    return rh >= datum::eps ? ah + (bh - ah) * exp(-rh / ch) : eh * exp(rh / fh) + dh; // eq. 33
}

double CDPM2::compute_plastic_potential_ratio(const double sigma, const double rho, const double kappap) const {
    double pgpsigma, pgprho;

    if(kappap < 1.) {
        double qh1;
        if(kappap >= 0.) {
            const auto kappap2 = kappap * kappap;
            const auto kappapa = kappap * kappap2 - 3. * kappap2 + 2. * kappap;

            qh1 = qh0 + one_minus_qh0 * (kappapa + kappap) - hp * kappapa; // eq. 30
            if(qh1 > 1.) qh1 = 1.;
        } else
            qh1 = qh0;

        const auto one_minus_qh1 = 1. - qh1;
        const auto qh1_qh1 = qh1 * qh1;

        const auto ag = 3. * ft_over_fc + half_m0; // eq. 28
        const auto dg = log(ag / (3. + half_m0)) + log_df;
        const auto bg = (1. + ft_over_fc) / 3. / dg; // eq. 29

        const auto af = rho / sqrt_six_fc;
        const auto cf = af + sigma / fc;
        const auto ef2 = 2. * (one_minus_qh1 * cf * cf + 3. * af);

        const auto pefpsigma = one_minus_qh1 * 2. * cf / fc;
        const auto pefprho = (one_minus_qh1 * 2. * cf + 3.) / sqrt_six_fc;

        pgpsigma = ef2 * pefpsigma + qh1_qh1 * ag * exp((sigma - ft / 3.) / bg / fc) / fc;
        pgprho = ef2 * pefprho + qh1_qh1 * m0 / sqrt_six_fc;
    } else {
        const auto qh2 = 1. + hp * (kappap - 1.);

        const auto ag = 3. * ft_over_fc * qh2 + half_m0;
        const auto dg = log(ag / (3. * qh2 + half_m0)) + log_df;
        const auto bg = qh2 * (1. + ft_over_fc) / 3. / dg;

        pgpsigma = ag * exp((sigma - qh2 * ft / 3.) / bg / fc) / fc;
        pgprho = 3. * rho / fc_fc + m0 / sqrt_six_fc;
    }

    return pgprho / pgpsigma * ratio_factor;
}

bool CDPM2::compute_regular_return(vec& solution, mat& jacobian, podarray<double>& status, const vec& target, const double theta) const {
    solution = target;

    const auto& new_sigma = solution(0);
    const auto& new_rho = solution(1);
    const auto& new_kappap = solution(2);
    const auto& new_lambda = solution(3);

    auto counter = 0;
    auto error = 1.;
    while(++counter < 20) {
        jacobian = compute_jacobian(status, new_lambda);

        const vec increment = solve(jacobian, compute_residual(status, solution, target));

        if((error = norm(increment)) < tolerance) break;
        suanpan_extra_debug("CDPM2 local plasticity iteration error: %.5E.\n", error);

        solution -= increment;

        status = compute_plasticity(new_sigma, new_rho, theta, new_kappap);
    }

    suanpan_debug("CDPM2 state determination loop counter: %u.\n", counter);

    // if the error is not far from the tolerance just continue
    if(counter == 20 && error > 1000. * tolerance) {
        suanpan_error("CDPM2 cannot converge within 20 iterations.\n");
        return false;
    }

    return true;
}

bool CDPM2::compute_vertex_return(double& sigma, const double rho, double& kappap) const {
    auto new_sigma = sigma;
    auto new_kappap = kappap + compute_incre_kappap(sigma, rho, new_sigma);

    const auto f_a = compute_yield_function(new_sigma, 0., 0., new_kappap);
    const auto f_b = compute_yield_function(0., 0., 0., kappap + compute_incre_kappap(sigma, rho, 0.));
    if(f_a < 0. || f_b > 0.) return false;

    // bisection iteration
    auto length = .5 * new_sigma;
    auto counter = 0;
    while(++counter < 100) {
        const auto sigma_middle = new_sigma - length;
        new_kappap = kappap + compute_incre_kappap(sigma, rho, sigma_middle);
        const auto f_middle = compute_yield_function(sigma_middle, 0., 0., new_kappap);
        if(f_middle > -tolerance) {
            new_sigma = sigma_middle;
            if(f_middle < tolerance) break;
        }
        length *= .5;
    }

    suanpan_extra_debug("CDPM2 bisection iteration counter: %u.\n", counter);

    if(counter == 100) return false;

    if(abs(rho / (sigma - new_sigma)) > compute_plastic_potential_ratio(new_sigma, rho, new_kappap)) return false;

    sigma = new_sigma;
    kappap = new_kappap;

    return true;
}

podarray<double> CDPM2::compute_damage(const double sigma, const double rho, const double theta) const {
    podarray<double> damage(2);

    const auto cos_theta = cos(theta);
    const auto cos_theta2 = 4. * cos_theta * cos_theta;
    const auto cos_theta3 = cos_theta2 * ea;

    auto eqv_strain = .5 * m0 * (rho / sqrt_six * (cos_theta3 + ec) / (2. * ea * cos_theta + eb * sqrt(cos_theta3 + ed)) + sigma);
    eqv_strain += sqrt(eqv_strain * eqv_strain + 1.5 * rho * rho);
    eqv_strain *= epsilon0 / fc;

    damage(0) = 0.;
    damage(1) = 0.;

    if(damage(0) < 0.) damage(0) = 0.;
    if(damage(1) < 0.) damage(1) = 0.;

    return damage;
}