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
 * @class DomainBase
 * @brief The DomainBase class is a template.
 *
 * The DomainBase is simply an abstraction of the Domain class. It provides all methods signature that are used in Domain class. The purpose is to split the declaration and implementation apart. As the Domain class is widely used in many other classes. The dependency hierarchy is simplified if replaced by the DomainBase.
 *
 * @author tlc
 * @date 01/10/2017
 * @version 0.2.0
 * @file DomainBase.h
 * @addtogroup Domain
 * @{
 */

#ifndef DOMAINBASE_H
#define DOMAINBASE_H

#include <Domain/Tag.h>
#include <unordered_set>

using std::map;
using std::unordered_set;
using std::vector;

template <typename T> class Factory;
class Amplitude;
class Constraint;
class Converger;
class Criterion;
class Element;
class ExternalModule;
class Integrator;
class Load;
class Material;
class Node;
class Recorder;
class Section;
class Solver;
class Step;

using AmplitudeQueue = vector<shared_ptr<Amplitude>>;
using ConstraintQueue = vector<shared_ptr<Constraint>>;
using ConvergerQueue = vector<shared_ptr<Converger>>;
using CriterionQueue = vector<shared_ptr<Criterion>>;
using ElementQueue = vector<shared_ptr<Element>>;
using IntegratorQueue = vector<shared_ptr<Integrator>>;
using LoadQueue = vector<shared_ptr<Load>>;
using MaterialQueue = vector<shared_ptr<Material>>;
using NodeQueue = vector<shared_ptr<Node>>;
using RecorderQueue = vector<shared_ptr<Recorder>>;
using SectionQueue = vector<shared_ptr<Section>>;
using SolverQueue = vector<shared_ptr<Solver>>;
using StepQueue = map<unsigned, shared_ptr<Step>>;

using LongFactory = Factory<double>;

class DomainBase : public Tag {
public:
    explicit DomainBase(unsigned);
    virtual ~DomainBase();

    virtual void set_factory(const shared_ptr<LongFactory>&) = 0;
    virtual const shared_ptr<LongFactory>& get_factory() const = 0;

    virtual bool insert(const shared_ptr<ExternalModule>&) = 0;
    virtual const vector<shared_ptr<ExternalModule>>& get_external_module_pool() const = 0;

    virtual bool insert(const shared_ptr<Amplitude>&) = 0;
    virtual bool insert(const shared_ptr<Constraint>&) = 0;
    virtual bool insert(const shared_ptr<Converger>&) = 0;
    virtual bool insert(const shared_ptr<Criterion>&) = 0;
    virtual bool insert(const shared_ptr<Element>&) = 0;
    virtual bool insert(const shared_ptr<Integrator>&) = 0;
    virtual bool insert(const shared_ptr<Load>&) = 0;
    virtual bool insert(const shared_ptr<Material>&) = 0;
    virtual bool insert(const shared_ptr<Node>&) = 0;
    virtual bool insert(const shared_ptr<Recorder>&) = 0;
    virtual bool insert(const shared_ptr<Section>&) = 0;
    virtual bool insert(const shared_ptr<Solver>&) = 0;
    virtual bool insert(const shared_ptr<Step>&) = 0;

    virtual bool erase_amplitude(unsigned) = 0;
    virtual bool erase_constraint(unsigned) = 0;
    virtual bool erase_converger(unsigned) = 0;
    virtual bool erase_criterion(unsigned) = 0;
    virtual bool erase_element(unsigned) = 0;
    virtual bool erase_integrator(unsigned) = 0;
    virtual bool erase_load(unsigned) = 0;
    virtual bool erase_material(unsigned) = 0;
    virtual bool erase_node(unsigned) = 0;
    virtual bool erase_recorder(unsigned) = 0;
    virtual bool erase_section(unsigned) = 0;
    virtual bool erase_solver(unsigned) = 0;
    virtual bool erase_step(unsigned) = 0;

    virtual void disable_amplitude(unsigned) = 0;
    virtual void disable_constraint(unsigned) = 0;
    virtual void disable_converger(unsigned) = 0;
    virtual void disable_criterion(unsigned) = 0;
    virtual void disable_element(unsigned) = 0;
    virtual void disable_integrator(unsigned) = 0;
    virtual void disable_load(unsigned) = 0;
    virtual void disable_material(unsigned) = 0;
    virtual void disable_node(unsigned) = 0;
    virtual void disable_recorder(unsigned) = 0;
    virtual void disable_section(unsigned) = 0;
    virtual void disable_solver(unsigned) = 0;
    virtual void disable_step(unsigned) = 0;

    virtual void enable_amplitude(unsigned) = 0;
    virtual void enable_constraint(unsigned) = 0;
    virtual void enable_converger(unsigned) = 0;
    virtual void enable_criterion(unsigned) = 0;
    virtual void enable_element(unsigned) = 0;
    virtual void enable_integrator(unsigned) = 0;
    virtual void enable_load(unsigned) = 0;
    virtual void enable_material(unsigned) = 0;
    virtual void enable_node(unsigned) = 0;
    virtual void enable_recorder(unsigned) = 0;
    virtual void enable_section(unsigned) = 0;
    virtual void enable_solver(unsigned) = 0;
    virtual void enable_step(unsigned) = 0;

    virtual const shared_ptr<Amplitude>& get_amplitude(unsigned) const = 0;
    virtual const shared_ptr<Constraint>& get_constraint(unsigned) const = 0;
    virtual const shared_ptr<Converger>& get_converger(unsigned) const = 0;
    virtual const shared_ptr<Criterion>& get_criterion(unsigned) const = 0;
    virtual const shared_ptr<Element>& get_element(unsigned) const = 0;
    virtual const shared_ptr<Integrator>& get_integrator(unsigned) const = 0;
    virtual const shared_ptr<Load>& get_load(unsigned) const = 0;
    virtual const shared_ptr<Material>& get_material(unsigned) const = 0;
    virtual const shared_ptr<Node>& get_node(unsigned) const = 0;
    virtual const shared_ptr<Recorder>& get_recorder(unsigned) const = 0;
    virtual const shared_ptr<Section>& get_section(unsigned) const = 0;
    virtual const shared_ptr<Solver>& get_solver(unsigned) const = 0;
    virtual const shared_ptr<Step>& get_step(unsigned) const = 0;

    virtual const AmplitudeQueue& get_amplitude_pool() const = 0;
    virtual const ConstraintQueue& get_constraint_pool() const = 0;
    virtual const ConvergerQueue& get_converger_pool() const = 0;
    virtual const CriterionQueue& get_criterion_pool() const = 0;
    virtual const ElementQueue& get_element_pool() const = 0;
    virtual const IntegratorQueue& get_integrator_pool() const = 0;
    virtual const LoadQueue& get_load_pool() const = 0;
    virtual const MaterialQueue& get_material_pool() const = 0;
    virtual const NodeQueue& get_node_pool() const = 0;
    virtual const RecorderQueue& get_recorder_pool() const = 0;
    virtual const SectionQueue& get_section_pool() const = 0;
    virtual const SolverQueue& get_solver_pool() const = 0;
    virtual const StepQueue& get_step_pool() const = 0;

    friend shared_ptr<Amplitude>& get_amplitude(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Constraint>& get_constraint(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Converger>& get_converger(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Criterion>& get_criterion(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Element>& get_element(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Integrator>& get_integrator(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Load>& get_load(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Material>& get_material(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Node>& get_node(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Recorder>& get_recorder(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Section>& get_section(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Solver>& get_solver(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Step>& get_step(const shared_ptr<DomainBase>&, unsigned);

    virtual size_t get_amplitude() const = 0;
    virtual size_t get_constraint() const = 0;
    virtual size_t get_converger() const = 0;
    virtual size_t get_criterion() const = 0;
    virtual size_t get_element() const = 0;
    virtual size_t get_integrator() const = 0;
    virtual size_t get_load() const = 0;
    virtual size_t get_material() const = 0;
    virtual size_t get_node() const = 0;
    virtual size_t get_recorder() const = 0;
    virtual size_t get_section() const = 0;
    virtual size_t get_solver() const = 0;
    virtual size_t get_step() const = 0;

    virtual bool find_amplitude(unsigned) const = 0;
    virtual bool find_constraint(unsigned) const = 0;
    virtual bool find_converger(unsigned) const = 0;
    virtual bool find_criterion(unsigned) const = 0;
    virtual bool find_element(unsigned) const = 0;
    virtual bool find_integrator(unsigned) const = 0;
    virtual bool find_load(unsigned) const = 0;
    virtual bool find_material(unsigned) const = 0;
    virtual bool find_node(unsigned) const = 0;
    virtual bool find_recorder(unsigned) const = 0;
    virtual bool find_section(unsigned) const = 0;
    virtual bool find_solver(unsigned) const = 0;
    virtual bool find_step(unsigned) const = 0;

    virtual void set_current_step_tag(unsigned) = 0;
    virtual void set_current_converger_tag(unsigned) = 0;
    virtual void set_current_integrator_tag(unsigned) = 0;
    virtual void set_current_solver_tag(unsigned) = 0;

    virtual unsigned get_current_step_tag() = 0;
    virtual unsigned get_current_converger_tag() = 0;
    virtual unsigned get_current_integrator_tag() = 0;
    virtual unsigned get_current_solver_tag() = 0;

    virtual const shared_ptr<Step>& get_current_step() const = 0;
    virtual const shared_ptr<Converger>& get_current_converger() const = 0;
    virtual const shared_ptr<Integrator>& get_current_integrator() const = 0;
    virtual const shared_ptr<Solver>& get_current_solver() const = 0;

    virtual bool insert_loaded_dof(unsigned) = 0;
    virtual bool insert_restrained_dof(unsigned) = 0;
    virtual bool insert_constrained_dof(unsigned) = 0;

    virtual const unordered_set<unsigned>& get_loaded_dof() const = 0;
    virtual const unordered_set<unsigned>& get_restrained_dof() const = 0;
    virtual const unordered_set<unsigned>& get_constrained_dof() const = 0;

    virtual const bool& is_updated() const = 0;

    virtual int initialize() = 0;
    virtual int initialize_load() = 0;

    virtual int process_load() = 0;
    virtual int process_constraint() = 0;
    virtual int process_criterion() = 0;
    virtual void record() = 0;
    virtual void enable_all() = 0;
    virtual void summary() const = 0;

    virtual void assemble_resistance() const = 0;
    virtual void assemble_damping() const = 0;
    virtual void assemble_initial_stiffness() const = 0;
    virtual void assemble_mass() const = 0;
    virtual void assemble_stiffness() const = 0;

    virtual void erase_machine_error() const = 0;

    virtual int update_current_status() const = 0;
    virtual int update_incre_status() const = 0;
    virtual int update_trial_status() const = 0;

    virtual void commit_status() const = 0;
    virtual void clear_status() const = 0;
    virtual void reset_status() const = 0;
};

#endif

//! @}
