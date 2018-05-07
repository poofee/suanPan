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
 * @class Step
 * @brief A Step class.
 * @author tlc
 * @date 27/08/2017
 * @version 0.2.1
 * @file Step.h
 * @addtogroup Step
 * @{
 */

#ifndef STEP_H
#define STEP_H

#include <Domain/Tag.h>

template <typename T> class Factory;
class DomainBase;
class Solver;
class Converger;
class Integrator;

class Step : public Tag {
    bool updated = false;

    bool symm_mat = true;
    bool band_mat = true;
    bool sparse_mat = false;

    double time_period = 1.0; /**< time period */

    double max_step_size = time_period; /**< maximum step size */
    double min_step_size = 1E-8;        /**< minimum step size */

    double ini_step_size = time_period; /**< initial step size */

    unsigned max_substep = 1000; /**< maximum increment number */

    bool fixed_step_size = false; /**< auto-stepping */

    unsigned solver_tag = 0;
    unsigned converger_tag = 0;
    unsigned integrator_tag = 0;

    weak_ptr<DomainBase> database;
    shared_ptr<Factory<double>> factory;
    shared_ptr<Solver> solver;
    shared_ptr<Converger> tester;
    shared_ptr<Integrator> modifier;

public:
    explicit Step(unsigned = 0, unsigned = CT_STEP, double = 1.);
    virtual ~Step();

    bool is_updated() const;

    virtual int initialize();

    virtual int analyze() = 0;

    void set_domain(const weak_ptr<DomainBase>&);
    const weak_ptr<DomainBase>& get_domain() const;

    void set_factory(const shared_ptr<Factory<double>>&);
    const shared_ptr<Factory<double>>& get_factory() const;

    void set_solver_tag(unsigned);
    void set_solver(const shared_ptr<Solver>&);
    const shared_ptr<Solver>& get_solver() const;

    void set_converger_tag(unsigned);
    void set_converger(const shared_ptr<Converger>&);
    const shared_ptr<Converger>& get_converger() const;

    void set_integrator_tag(unsigned);
    void set_integrator(const shared_ptr<Integrator>&);
    const shared_ptr<Integrator>& get_integrator() const;

    void set_time_perid(double);
    double get_time_period() const;

    void set_ini_step_size(double);
    void set_min_step_size(double);
    void set_max_step_size(double);
    void set_max_substep(unsigned);

    double get_ini_step_size() const;
    double get_min_step_size() const;
    double get_max_step_size() const;
    unsigned get_max_substep() const;

    bool is_fixed_step_size() const;
    void set_fixed_step_size(bool);

    bool is_symm() const;
    bool is_band() const;
    bool is_sparse() const;
    void set_symm(bool);
    void set_band(bool);
    void set_sparse(bool);
};

#endif

//! @}
