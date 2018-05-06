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

#include "Mindlin.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/shapeFunction.hpp>
#include <Toolbox/tensorToolbox.h>
#include <Toolbox/utility.h>
#include <utility>

const unsigned Mindlin::m_node = 4;
const unsigned Mindlin::m_dof = 3;
const unsigned Mindlin::m_size = m_dof * m_node;
mat Mindlin::iso_mapping;

/**
 * \brief create converter for resultant forces
 * \param E edge label
 * \param T element thickness
 * \param NP vector of node pointers
 * \param IP integration plan
 * \param TRANS transformation matrix from parent to global
 */
Mindlin::ResultantConverter::ResultantConverter(const unsigned E, const double T, const vector<weak_ptr<Node>>& NP, const IntegrationPlan& IP, const mat& TRANS)
    : direction_cosine(2, 3) {
    vec node_i, node_j;

    if(E == 1) {
        node_i = NP[0].lock()->get_coordinate();
        node_j = NP[1].lock()->get_coordinate();
    } else if(E == 2) {
        node_i = NP[1].lock()->get_coordinate();
        node_j = NP[2].lock()->get_coordinate();
    } else if(E == 3) {
        node_i = NP[2].lock()->get_coordinate();
        node_j = NP[3].lock()->get_coordinate();
    } else {
        node_i = NP[3].lock()->get_coordinate();
        node_j = NP[0].lock()->get_coordinate();
    }

    auto cos_edge = node_j(1) - node_i(1);
    auto sin_edge = node_i(0) - node_j(0);

    const auto edge_length = sqrt(cos_edge * cos_edge + sin_edge * sin_edge);

    cos_edge /= edge_length;
    sin_edge /= edge_length;

    if(cos_edge > 1.) cos_edge = 1.;
    if(cos_edge < -1.) cos_edge = -1.;
    if(sin_edge > 1.) sin_edge = 1.;
    if(sin_edge < -1.) sin_edge = -1.;

    direction_cosine(1, 2) = (direction_cosine(0, 0) = cos_edge * cos_edge) - (direction_cosine(0, 1) = sin_edge * sin_edge);
    direction_cosine(0, 2) = -2. * (direction_cosine(1, 0) = -(direction_cosine(1, 1) = cos_edge * sin_edge));

    const auto& X1 = IP(0, 0);
    const auto& X2 = IP(1, 0);

    vec int_pt_a, int_pt_b;

    if(E == 1) {
        int_pt_a = TRANS * vec(std::initializer_list<double>{ 1., X1, -1., -X1 });
        int_pt_b = TRANS * vec(std::initializer_list<double>{ 1., X2, -1., -X2 });
    } else if(E == 2) {
        int_pt_a = TRANS * vec(std::initializer_list<double>{ 1., 1., X1, X1 });
        int_pt_b = TRANS * vec(std::initializer_list<double>{ 1., 1., X2, X2 });
    } else if(E == 3) {
        int_pt_a = TRANS * vec(std::initializer_list<double>{ 1., X1, 1., X1 });
        int_pt_b = TRANS * vec(std::initializer_list<double>{ 1., X2, 1., X2 });
    } else {
        int_pt_a = TRANS * vec(std::initializer_list<double>{ 1., -1., X1, -X1 });
        int_pt_b = TRANS * vec(std::initializer_list<double>{ 1., -1., X2, -X2 });
    }

    auto weight = .5 * edge_length * T;

    if(E == 3 || E == 4) weight = -weight;

    const mat part_a = shape::stress11(int_pt_a) * IP(0, 1) * weight;
    const mat part_b = shape::stress11(int_pt_b) * IP(1, 1) * weight;
    converter_a = part_a + part_b;
    converter_b = part_a * .5 * edge_length * X1 + part_b * .5 * edge_length * X2;
}

double Mindlin::ResultantConverter::F(const vec& alpha) const { return dot(direction_cosine.row(0), converter_a * alpha); }

double Mindlin::ResultantConverter::V(const vec& alpha) const { return dot(direction_cosine.row(1), converter_a * alpha); }

double Mindlin::ResultantConverter::M(const vec& alpha) const { return dot(direction_cosine.row(0), converter_b * alpha); }

Mindlin::IntegrationPoint::IntegrationPoint(vec C, const double F, unique_ptr<Material>&& M)
    : coor(std::move(C))
    , factor(F)
    , m_material(std::move(M)) {}

mat Mindlin::form_interpolation_displacement(const mat& pn_pxy, const mat& pnt_pxy) {
    mat poly_disp(3, m_size, fill::zeros);

    for(unsigned J = 0; J < m_node; ++J) {
        poly_disp(0, m_dof * J) = poly_disp(2, m_dof * J + 1) = pn_pxy(0, J);
        poly_disp(2, m_dof * J) = poly_disp(1, m_dof * J + 1) = pn_pxy(1, J);
        poly_disp(0, m_dof * J + 2) = pnt_pxy(0, J);
        poly_disp(1, m_dof * J + 2) = pnt_pxy(1, J + 4);
        poly_disp(2, m_dof * J + 2) = pnt_pxy(0, J + 4) + pnt_pxy(1, J);
    }

    return poly_disp;
}

mat Mindlin::form_interpolation_enhanced_strain(const mat& pn_pxy) {
    mat poly_enhanced_strain(3, 3);

    poly_enhanced_strain(0, 1) = poly_enhanced_strain(1, 0) = 0.;
    poly_enhanced_strain(2, 1) = poly_enhanced_strain(0, 0) = pn_pxy(0, 0);
    poly_enhanced_strain(2, 0) = poly_enhanced_strain(1, 1) = pn_pxy(1, 0);
    poly_enhanced_strain(0, 2) = pn_pxy(0, 1);
    poly_enhanced_strain(1, 2) = pn_pxy(1, 1);
    poly_enhanced_strain(2, 2) = pn_pxy(0, 1) + pn_pxy(1, 1);

    return poly_enhanced_strain;
}

Mindlin::Mindlin(const unsigned T, const uvec& N, const unsigned M, const double TH, const char IP, const double RHOX, const unsigned MATX, const double RHOY, const unsigned MATY)
    : MaterialElement(T, ET_MINDLIN, m_node, m_dof, N, uvec{ M })
    , thickness(TH)
    , int_scheme(IP)
    , rho_x(RHOX)
    , rho_y(RHOY)
    , mat_x(MATX)
    , mat_y(MATY) {
    if(iso_mapping.is_empty()) {
        mat t_mapping(4, 4);
        t_mapping.fill(.25);
        t_mapping(1, 0) = t_mapping(1, 3) = t_mapping(2, 0) = t_mapping(2, 1) = t_mapping(3, 1) = t_mapping(3, 3) = -.25;
        iso_mapping = t_mapping;
    }
}

void Mindlin::initialize(const shared_ptr<DomainBase>& D) {
    mat ele_coor(m_node, 2);
    for(unsigned I = 0; I < m_node; ++I) {
        auto& t_coor = node_ptr[I].lock()->get_coordinate();
        for(auto J = 0; J < 2; ++J) ele_coor(I, J) = t_coor(J);
    }

    trans_mat = trans(iso_mapping * ele_coor);

    const IntegrationPlan edge_plan(1, 2, IntegrationType::GAUSS);
    edge.clear(), edge.reserve(4);
    for(unsigned I = 0; I < 4; ++I) edge.emplace_back(I + 1, thickness, node_ptr, edge_plan, trans_mat);

    const auto LX1 = ele_coor(1, 1) - ele_coor(0, 1);
    const auto LX2 = ele_coor(2, 1) - ele_coor(1, 1);
    const auto LX3 = ele_coor(3, 1) - ele_coor(2, 1);
    const auto LX4 = ele_coor(0, 1) - ele_coor(3, 1);
    const auto LY1 = ele_coor(0, 0) - ele_coor(1, 0);
    const auto LY2 = ele_coor(1, 0) - ele_coor(2, 0);
    const auto LY3 = ele_coor(2, 0) - ele_coor(3, 0);
    const auto LY4 = ele_coor(3, 0) - ele_coor(0, 0);

    auto& material_proto = D->get_material(unsigned(material_tag(0)));

    if(material_proto->material_type == MaterialType::D2 && std::dynamic_pointer_cast<Material2D>(material_proto)->plane_type == PlaneType::E) suanpan::hacker(thickness) = 1.;

    auto ini_stiffness = material_proto->get_initial_stiffness();

    if(reinforced_x && D->find_material(mat_x))
        ini_stiffness(0, 0) += D->get_material(mat_x)->get_initial_stiffness().at(0) * rho_x;
    else {
        suanpan_error("cannot find the reinforcement material.\n");
        D->disable_element(get_tag());
        return;
    }
    if(reinforced_y && D->find_material(mat_y))
        ini_stiffness(1, 1) += D->get_material(mat_y)->get_initial_stiffness().at(0) * rho_y;
    else {
        suanpan_error("cannot find the reinforcement material.\n");
        D->disable_element(get_tag());
        return;
    }

    const auto poissons_ratio = material_proto->get_parameter(ParameterType::POISSONSRATIO);

    const IntegrationPlan plan(2, int_scheme == 'I' ? 2 : 3, int_scheme == 'I' ? IntegrationType::IRONS : int_scheme == 'L' ? IntegrationType::LOBATTO : IntegrationType::GAUSS);

    mat pnt(2, 8), pne(2, 2), H(11, 11, fill::zeros), HTT(11, 11, fill::zeros);

    vec disp_mode(4);
    disp_mode(0) = 1.;

    N.zeros(11, 12), M.zeros(11, 3);

    int_pt.clear(), int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        vec t_vec{ plan(I, 0), plan(I, 1) };
        const auto pn = shape::quad(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(t_vec, det(jacob) * plan(I, 2) * thickness, material_proto->get_copy());

        if(reinforced_x) int_pt[I].rebar_x = D->get_material(mat_x)->get_copy();
        if(reinforced_y) int_pt[I].rebar_y = D->get_material(mat_y)->get_copy();

        const auto& X = int_pt[I].coor(0);
        const auto& Y = int_pt[I].coor(1);

        disp_mode(1) = X, disp_mode(2) = Y, disp_mode(3) = X * Y;

        const vec coord = trans_mat * disp_mode;

        int_pt[I].poly_stress = shape::stress11(coord);

        int_pt[I].poly_strain = poissons_ratio == 0. ? solve(ini_stiffness, int_pt[I].poly_stress) : shape::strain11(coord, poissons_ratio);

        const auto X2 = 2. * X, Y2 = 2. * Y, XP = X + 1., XM = X - 1., YP = Y + 1., YM = Y - 1.;

        pnt(0, 0) = YM * (LX4 * YP - LX1 * X2);
        pnt(0, 1) = YM * (LX2 * YP + LX1 * X2);
        pnt(0, 2) = YP * (LX3 * X2 - LX2 * YM);
        pnt(0, 3) = -YP * (LX3 * X2 + LX4 * YM);
        pnt(0, 4) = YM * (LY4 * YP - LY1 * X2);
        pnt(0, 5) = YM * (LY2 * YP + LY1 * X2);
        pnt(0, 6) = YP * (LY3 * X2 - LY2 * YM);
        pnt(0, 7) = -YP * (LY3 * X2 + LY4 * YM);
        pnt(1, 0) = XM * (LX4 * Y2 - LX1 * XP);
        pnt(1, 1) = XP * (LX1 * XM + LX2 * Y2);
        pnt(1, 2) = XP * (LX3 * XM - LX2 * Y2);
        pnt(1, 3) = -XM * (LX3 * XP + LX4 * Y2);
        pnt(1, 4) = XM * (LY4 * Y2 - LY1 * XP);
        pnt(1, 5) = XP * (LY1 * XM + LY2 * Y2);
        pnt(1, 6) = XP * (LY3 * XM - LY2 * Y2);
        pnt(1, 7) = -XM * (LY3 * XP + LY4 * Y2);

        pne(0, 0) = X * X - Y * Y + X;
        pne(1, 0) = Y * Y - X * X + Y;
        pne(0, 1) = 6. * X * Y + 3. * Y * Y - 1. + X;
        pne(1, 1) = 6. * X * Y + 3. * X * X - 1. + Y;

        const mat t_mat = int_pt[I].poly_stress.t() * int_pt[I].factor;
        M += t_mat * form_interpolation_enhanced_strain(solve(jacob, pne));
        N += t_mat * form_interpolation_displacement(solve(jacob, pn), solve(jacob, pnt / 16.));
        H += t_mat * int_pt[I].poly_strain;
        HTT += int_pt[I].poly_strain.t() * ini_stiffness * int_pt[I].poly_strain * int_pt[I].factor;
    }

    initial_mass.zeros(m_size, m_size);
    const auto t_density = material_proto->get_parameter();
    if(t_density != 0.) {
        for(const auto& I : int_pt) {
            const auto n_int = shape::quad(I.coor, 0);
            const auto tmp_a = t_density * I.factor;
            for(unsigned J = 0; J < m_node; ++J)
                for(auto K = J; K < m_node; ++K) initial_mass(m_dof * J, m_dof * K) += tmp_a * n_int(J) * n_int(K);
        }
        for(unsigned I = 0; I < m_size; I += m_dof) {
            initial_mass(I + 1, I + 1) = initial_mass(I, I);
            for(auto J = I + m_dof; J < m_size; J += m_dof) initial_mass(J, I) = initial_mass(I + 1, J + 1) = initial_mass(J + 1, I + 1) = initial_mass(I, J);
        }
    }
    trial_mass = current_mass = initial_mass;

    HT = trans(H);

    solve(NT, H, N), solve(MT, H, M);

    const mat T = HTT * MT, W = NT.t() * T;

    trial_stiffness = current_stiffness = initial_stiffness = NT.t() * HTT * NT - W * (trial_viwt = current_viwt = initial_viwt = solve(mat(MT.t() * T), W.t()));

    pre_disp.zeros(m_size);

    current_vif.zeros(3);
    current_zeta.zeros(3);
    current_beta.zeros(11);
    current_alpha.zeros(11);

    trial_vif.zeros(3);
    trial_zeta.zeros(3);
    trial_beta.zeros(11);
    trial_alpha.zeros(11);
}

int Mindlin::update_status() {
    auto idx = 0;
    vec new_disp(m_size);
    for(const auto& t_ptr : node_ptr) {
        auto& t_disp = t_ptr.lock()->get_incre_displacement();
        for(unsigned pos = 0; pos < m_dof; ++pos) new_disp(idx++) = t_disp(pos);
    }

    const vec incre_disp = new_disp - pre_disp;
    const vec incre_zeta = -trial_viwt * incre_disp - trial_vif;

    trial_zeta += incre_zeta;
    trial_beta += NT * incre_disp + MT * incre_zeta;

    vec local_stress(11, fill::zeros);
    mat local_stiffness(11, 11, fill::zeros);
    for(const auto& t_pt : int_pt) {
        const vec t_strain = t_pt.poly_strain * trial_beta;
        if(t_pt.m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        if(!reinforced_x && !reinforced_y) {
            local_stress += t_pt.poly_strain.t() * t_pt.factor * t_pt.m_material->get_trial_stress();
            local_stiffness += t_pt.poly_strain.t() * t_pt.factor * t_pt.m_material->get_trial_stiffness() * t_pt.poly_strain;
        } else {
            auto t_stress = t_pt.m_material->get_trial_stress();
            auto t_stiffness = t_pt.m_material->get_trial_stiffness();
            if(reinforced_x) {
                t_pt.rebar_x->update_trial_status(t_strain(0));
                t_stress(0) += t_pt.rebar_x->get_trial_stress().at(0) * rho_x;
                t_stiffness(0, 0) += t_pt.rebar_x->get_trial_stiffness().at(0) * rho_x;
            }
            if(reinforced_y) {
                t_pt.rebar_y->update_trial_status(t_strain(1));
                t_stress(1) += t_pt.rebar_y->get_trial_stress().at(0) * rho_y;
                t_stiffness(1, 1) += t_pt.rebar_y->get_trial_stiffness().at(0) * rho_y;
            }
            local_stress += t_pt.poly_strain.t() * t_pt.factor * t_stress;
            local_stiffness += t_pt.poly_strain.t() * t_pt.factor * t_stiffness * t_pt.poly_strain;
        }
    }

    const mat T = NT.t() * local_stiffness, V = MT.t() * local_stiffness * MT, W = T * MT;

    if(!solve(trial_alpha, HT, local_stress) || !solve(trial_viwt, V, W.t()) || !solve(trial_vif, V, M.t() * trial_alpha)) return SUANPAN_FAIL;

    trial_resistance = N.t() * trial_alpha - W * trial_vif;
    trial_stiffness = T * (NT - MT * trial_viwt);

    pre_disp = new_disp;

    return SUANPAN_SUCCESS;
}

int Mindlin::commit_status() {
    current_zeta = trial_zeta;
    current_beta = trial_beta;
    current_alpha = trial_alpha;
    current_vif = trial_vif;
    current_viwt = trial_viwt;

    pre_disp.zeros();

    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    if(reinforced_x)
        for(const auto& I : int_pt) code += I.rebar_x->commit_status();
    if(reinforced_y)
        for(const auto& I : int_pt) code += I.rebar_y->commit_status();
    return code;
}

int Mindlin::clear_status() {
    current_zeta.zeros();
    trial_zeta.zeros();
    current_beta.zeros();
    trial_beta.zeros();
    current_alpha.zeros();
    trial_alpha.zeros();
    current_vif.zeros();
    trial_vif.zeros();

    pre_disp.zeros();

    current_viwt = trial_viwt = initial_viwt;

    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    if(reinforced_x)
        for(const auto& I : int_pt) code += I.rebar_x->clear_status();
    if(reinforced_y)
        for(const auto& I : int_pt) code += I.rebar_y->clear_status();
    return code;
}

int Mindlin::reset_status() {
    trial_zeta = current_zeta;
    trial_beta = current_beta;
    trial_alpha = current_alpha;
    trial_vif = current_vif;
    trial_viwt = current_viwt;

    pre_disp.zeros();

    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    if(reinforced_x)
        for(const auto& I : int_pt) code += I.rebar_x->reset_status();
    if(reinforced_y)
        for(const auto& I : int_pt) code += I.rebar_y->reset_status();
    return code;
}

vector<vec> Mindlin::record(const OutputType& T) {
    vector<vec> data;

    switch(T) {
    case OutputType::S:
        for(const auto& I : int_pt) data.emplace_back(I.poly_stress * current_alpha);
        break;
    case OutputType::S11:
        for(const auto& I : int_pt) data.emplace_back(I.poly_stress.row(0) * current_alpha);
        break;
    case OutputType::S22:
        for(const auto& I : int_pt) data.emplace_back(I.poly_stress.row(1) * current_alpha);
        break;
    case OutputType::S12:
        for(const auto& I : int_pt) data.emplace_back(I.poly_stress.row(2) * current_alpha);
        break;
    case OutputType::SP:
        for(const auto& I : int_pt) data.emplace_back(transform::stress::principal(I.poly_stress * current_alpha));
        break;
    case OutputType::SP1:
        for(const auto& I : int_pt) data.emplace_back(vec{ transform::stress::principal(I.poly_stress * current_alpha).at(0) });
        break;
    case OutputType::SP2:
        for(const auto& I : int_pt) data.emplace_back(vec{ transform::stress::principal(I.poly_stress * current_alpha).at(1) });
        break;
    case OutputType::SINT:
        data.emplace_back(current_alpha);
        break;
    case OutputType::E:
        for(const auto& I : int_pt) data.emplace_back(I.poly_strain * current_beta);
        break;
    case OutputType::E11:
        for(const auto& I : int_pt) data.emplace_back(I.poly_strain.row(0) * current_beta);
        break;
    case OutputType::E22:
        for(const auto& I : int_pt) data.emplace_back(I.poly_strain.row(1) * current_beta);
        break;
    case OutputType::E12:
        for(const auto& I : int_pt) data.emplace_back(I.poly_strain.row(2) * current_beta);
        break;
    case OutputType::EP:
        for(const auto& I : int_pt) data.emplace_back(transform::strain::principal(I.poly_strain * current_beta));
        break;
    case OutputType::EP1:
        for(const auto& I : int_pt) data.emplace_back(vec{ transform::strain::principal(I.poly_strain * current_beta).at(0) });
        break;
    case OutputType::EP2:
        for(const auto& I : int_pt) data.emplace_back(vec{ transform::strain::principal(I.poly_strain * current_beta).at(1) });
        break;
    case OutputType::EINT:
        data.emplace_back(current_beta);
        break;
    case OutputType::RESULTANT:
        for(const auto& I : edge) data.emplace_back(vec(std::initializer_list<double>{ I.F(current_alpha), I.V(current_alpha), I.M(current_alpha) }));
        break;
    case OutputType::AXIAL:
        data.emplace_back(vec(std::initializer_list<double>{ edge[0].F(current_alpha), edge[1].F(current_alpha), edge[2].F(current_alpha), edge[3].F(current_alpha) }));
        break;
    case OutputType::SHEAR:
        data.emplace_back(vec(std::initializer_list<double>{ edge[0].V(current_alpha), edge[1].V(current_alpha), edge[2].V(current_alpha), edge[3].V(current_alpha) }));
        break;
    case OutputType::MOMENT:
        data.emplace_back(vec(std::initializer_list<double>{ edge[0].M(current_alpha), edge[1].M(current_alpha), edge[2].M(current_alpha), edge[3].M(current_alpha) }));
        break;
    default:
        for(const auto& I : int_pt)
            for(auto J : I.m_material->record(T)) data.emplace_back(J);
        break;
    }
    return data;
}

void Mindlin::print() {
    suanpan_info("Mindlin mixed quad element %u connects nodes:\n", get_tag());
    node_encoding.t().print();
    suanpan_info("\nMaterial model response:");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("\nIntegration Point %u:\t", I + 1);
        int_pt[I].coor.t().print();
        int_pt[I].m_material->print();
    }
    suanpan_info("\nElement model response:");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("\nIntegration Point %u:\t", I + 1);
        int_pt[I].coor.t().print();
        suanpan_info("Strain:\n");
        (int_pt[I].poly_strain * current_beta).t().print();
        suanpan_info("Stress:\n");
        (int_pt[I].poly_stress * current_alpha).t().print();
    }
}
