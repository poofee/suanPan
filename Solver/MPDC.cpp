#include "MPDC.h"
#include <Converger/Converger.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>

MPDC::MPDC(const unsigned T)
    : Solver(T, CT_MPDC) {}

int MPDC::analyze() {
    auto& C = get_converger();
    auto& G = get_integrator();
    auto& W = G->get_domain().lock()->get_factory();

    if(W->get_storage_scheme() != StorageScheme::BAND) suanpan_warning("band matrix is sugguested for displacement controlled algorithm.\n");

    suanpan_info("current analysis time: %.5f.\n", W->get_trial_time());

    auto& max_iteration = C->get_max_iteration();

    // ninja anchor
    auto& t_ninja = get_ninja(W);

    auto& load_ref = W->get_reference_load();

    // get column index for each nonzero dof
    // uvec load_ref_idx = find(load_ref);
    // for(auto I = 0; I < load_ref_idx.n_elem; ++I) load_ref_idx(I) -= I * load_ref.n_rows;

    auto& load_ref_idx = W->get_reference_dof();

    mat disp_a;

    // iteration counter
    unsigned counter = 0;

    while(true) {
        // assemble resistance
        G->assemble_resistance();
        // assemble stiffness
        G->assemble_matrix();
        // process loads
        G->process_load();
        // process constraints
        G->process_constraint();

        // solve ninja
        auto flag = W->get_stiffness()->solve(t_ninja, load_ref * W->get_trial_load_factor() + W->get_trial_load() - W->get_sushi());
        // make sure lapack solver succeeds
        if(flag != 0) return flag;
        // solve reference displacement
        flag = W->get_stiffness()->solve_trs(disp_a, load_ref);
        // make sure lapack solver succeeds
        if(flag != 0) return flag;

        vec incre_lambda = -solve(mat(disp_a.rows(load_ref_idx)), t_ninja.rows(load_ref_idx));

        if(counter == 0) incre_lambda += solve(mat(disp_a.rows(load_ref_idx)), W->get_incre_settlement().rows(load_ref_idx));

        t_ninja += disp_a * incre_lambda;

        // avoid machine error accumulation
        G->erase_machine_error();
        // update trial load factor
        W->update_trial_load_factor(W->get_trial_load_factor() + incre_lambda);
        // update trial displacement
        W->update_trial_displacement(W->get_trial_displacement() + t_ninja);
        // update for nodes and elements
        if(G->update_trial_status() != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        // exit if converged
        if(C->is_converged()) return SUANPAN_SUCCESS;
        // exit if maximum iteration is hit
        if(++counter > max_iteration) return SUANPAN_FAIL;
    }
}
