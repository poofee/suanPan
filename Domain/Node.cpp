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

#include "Node.h"
#include <Domain/DomainBase.h>
#include <Recorder/OutputType.h>
#include <utility>

/**
 * \brief Empty constructor.
 * \param T `unique_tag`
 */
Node::Node(const unsigned& T)
    : Tag(T, CT_NODE) {
    suanpan_debug("Node %u ctor() called.\n", T);
}

/**
 * \brief Initialize `coordinate` and set `num_dof` to the size of `coordinate`.
 * \param T `unique_tag`
 * \param C `coordinate`
 */
Node::Node(const unsigned& T, const vec& C)
    : Tag(T, CT_NODE)
    , num_dof(static_cast<unsigned>(C.n_elem))
    , coordinate(C) {
    suanpan_debug("Node %u ctor() called.\n", T);
}

/**
 * \brief Initialize `num_dof` and set the size of `coordinate` to `num_dof`.
 * \param T `unique_tag`
 * \param D `num_dof`
 */
Node::Node(const unsigned& T, const unsigned& D)
    : Tag(T, CT_NODE)
    , num_dof(D)
    , coordinate(num_dof) {
    suanpan_debug("Node %u ctor() called.\n", T);
}

/**
 * \brief Initialize `num_dof` and `coordinate`.
 * \param T `unique_tag`
 * \param D `num_dof`
 * \param C `coordinate`
 */
Node::Node(const unsigned& T, const unsigned& D, vec C)
    : Tag(T, CT_NODE)
    , num_dof(D)
    , coordinate(std::move(C)) {
    suanpan_debug("Node %u ctor() called.\n", T);
}

/**
 * \brief Default destructor.
 */
Node::~Node() { suanpan_debug("Node %u dtor() called.\n", get_tag()); }

/**
 * \brief This method should be called after Element objects are set. Element objects will set the minimum number of DoFs for all related Node objects. This method initialize all member variables with the size of `num_dof` and fill `original_dof` with `-1` to indicated it should be omitted from the system. Finally check if the size of `coordinate` is the same of `num_dof`, if not, resize it to `num_dof`. This will be necessary for beam/plate/shell problems which have more DoFs than coordinates.
 */
void Node::initialize(const shared_ptr<DomainBase>& D) {
    if(initialized || !is_active()) return;

    if(num_dof != 0) {
        original_dof.zeros(num_dof);
        original_dof.fill(static_cast<uword>(-1));

        reordered_dof.reset();

        current_resistance.resize(num_dof);
        current_displacement.resize(num_dof);
        current_velocity.resize(num_dof);
        current_acceleration.resize(num_dof);

        incre_resistance.resize(num_dof);
        incre_displacement.resize(num_dof);
        incre_velocity.resize(num_dof);
        incre_acceleration.resize(num_dof);

        trial_resistance.resize(num_dof);
        trial_displacement.resize(num_dof);
        trial_velocity.resize(num_dof);
        trial_acceleration.resize(num_dof);

        // if(num_dof > coordinate.n_elem) coordinate.resize(num_dof);
    } else {
        suanpan_debug("Node %u is not used in the problem, now disable it.\n", get_tag());
        D->disable_node(get_tag());
    }

    initialized = true;
}

/**
 * \brief Method to set `num_dof`.
 * \param D `num_dof`
 */
void Node::set_dof_number(const unsigned& D) {
    if(num_dof == D) return;

    num_dof = D;
    initialized = false;
}

/**
 * \brief Method to return `num_dof`.
 * \return `num_dof`
 */
const unsigned& Node::get_dof_number() const { return num_dof; }

/**
 * \brief Method to set `original_dof`.
 * \param F Current Index Counter
 */
void Node::set_original_dof(unsigned& F) {
    if(!is_active()) return;

    for(unsigned I = 0; I < num_dof; ++I, ++F)
        if(original_dof(I) != F) {
            original_dof(I) = F;
            initialized = false;
        }
}

/**
 * \brief Method to set `original_dof`.
 * \param D `original_dof`
 */
void Node::set_original_dof(const uvec& D) {
    if(original_dof.size() != D.size()) {
        original_dof = D;
        initialized = false;
    } else {
        auto code = 0;
        for(uword I = 0; I < D.size(); ++I)
            if(original_dof(I) != D(I)) {
                ++code;
                break;
            }
        if(code != 0) {
            original_dof = D;
            initialized = false;
        }
    }
}

/**
 * \brief Method to return `original_dof`.
 * \return `original_dof`
 */
const uvec& Node::get_original_dof() const { return original_dof; }

/**
 * \brief Method to set `reordered_dof`.
 * \param R `reordered_dof`
 */
void Node::set_reordered_dof(const uvec& R) { reordered_dof = R; }

/**
 * \brief Method to return `reordered_dof`.
 * \return `reordered_dof`
 */
const uvec& Node::get_reordered_dof() const {
    if(reordered_dof.is_empty()) return original_dof;
    return reordered_dof;
}

/**
 * \brief Method to set `coordinate`.
 * \param C `coordinate`
 */
void Node::set_coordinate(const vec& C) { coordinate = C; }

/**
 * \brief Method to return `coordinate`.
 * \return `coordinate`
 */
const vec& Node::get_coordinate() const { return coordinate; }

void Node::set_current_resistance(const vec& R) { current_resistance = R; }

/**
 * \brief Method to set variable independently.
 * \param D `current_displacement`
 */
void Node::set_current_displacement(const vec& D) { current_displacement = D; }

/**
 * \brief Method to set variable independently.
 * \param V `current_velocity`
 */
void Node::set_current_velocity(const vec& V) { current_velocity = V; }

/**
 * \brief Method to set variable independently.
 * \param A `current_acceleration`
 */
void Node::set_current_acceleration(const vec& A) { current_acceleration = A; }

void Node::set_incre_resistance(const vec& R) { incre_resistance = R; }

/**
 * \brief Method to set variable independently.
 * \param D `incre_displacement`
 */
void Node::set_incre_displacement(const vec& D) { incre_displacement = D; }

/**
 * \brief Method to set variable independently.
 * \param V `incre_velocity`
 */
void Node::set_incre_velocity(const vec& V) { incre_velocity = V; }

/**
 * \brief Method to set variable independently.
 * \param A `incre_acceleration`
 */
void Node::set_incre_acceleration(const vec& A) { incre_acceleration = A; }

void Node::set_trial_resistance(const vec& R) { trial_resistance = R; }

/**
 * \brief Method to set variable independently.
 * \param D `trial_displacement`
 */
void Node::set_trial_displacement(const vec& D) { trial_displacement = D; }

/**
 * \brief Method to set variable independently.
 * \param V `trial_velocity`
 */
void Node::set_trial_velocity(const vec& V) { trial_velocity = V; }

/**
 * \brief Method to set variable independently.
 * \param A `trial_acceleration`
 */
void Node::set_trial_acceleration(const vec& A) { trial_acceleration = A; }

const vec& Node::get_current_resistance() const { return current_resistance; }

/**
 * \brief Method to return `current_displacement`.
 * \return `current_displacement`
 */
const vec& Node::get_current_displacement() const { return current_displacement; }

/**
 * \brief Method to return `current_velocity`.
 * \return `current_velocity`
 */
const vec& Node::get_current_velocity() const { return current_velocity; }

/**
 * \brief Method to return `current_acceleration`.
 * \return `current_acceleration`
 */
const vec& Node::get_current_acceleration() const { return current_acceleration; }

const vec& Node::get_incre_resistance() const { return incre_resistance; }

/**
 * \brief Method to return `incre_displacement`.
 * \return `incre_displacement`
 */
const vec& Node::get_incre_displacement() const { return incre_displacement; }

/**
 * \brief Method to return `incre_velocity`.
 * \return `incre_velocity`
 */
const vec& Node::get_incre_velocity() const { return incre_velocity; }

/**
 * \brief Method to return `incre_acceleration`.
 * \return `incre_acceleration`
 */
const vec& Node::get_incre_acceleration() const { return incre_acceleration; }

const vec& Node::get_trial_resistance() const { return trial_resistance; }

/**
 * \brief Method to return `trial_displacement`.
 * \return `trial_displacement`
 */
const vec& Node::get_trial_displacement() const { return trial_displacement; }

/**
 * \brief Method to return `trial_velocity`.
 * \return `trial_velocity`
 */
const vec& Node::get_trial_velocity() const { return trial_velocity; }

/**
 * \brief Method to return `trial_acceleration`.
 * \return `trial_acceleration`
 */
const vec& Node::get_trial_acceleration() const { return trial_acceleration; }

void Node::update_current_resistance(const vec& R) {
    trial_resistance = current_resistance = R(reordered_dof);
    incre_resistance.zeros();
}

void Node::update_incre_resistance(const vec& R) {
    incre_resistance = R(reordered_dof);
    trial_resistance = current_resistance + incre_resistance;
}

void Node::update_trial_resistance(const vec& R) {
    trial_resistance = R(reordered_dof);
    incre_resistance = trial_resistance - current_resistance;
}

void Node::update_current_status(const vec& D) {
    trial_displacement = current_displacement = D(reordered_dof);
    incre_displacement.zeros();
}

void Node::update_current_status(const vec& D, const vec& V) {
    trial_velocity = current_velocity = V(reordered_dof);
    incre_velocity.zeros();
    update_current_status(D);
}

void Node::update_current_status(const vec& D, const vec& V, const vec& A) {
    trial_acceleration = current_acceleration = A(reordered_dof);
    incre_acceleration.zeros();
    update_current_status(D, V);
}

/**
 * \brief Method to update displacement.
 * \param D `incre_displacement`
 */
void Node::update_incre_status(const vec& D) {
    incre_displacement = D(reordered_dof);
    trial_displacement = current_displacement + incre_displacement;
}

/**
 * \brief Method to update velocity and call previous method to further update displacement. This is used in Dynamic analysis.
 * \param D `incre_displacement`
 * \param V `incre_velocity`
 */
void Node::update_incre_status(const vec& D, const vec& V) {
    incre_velocity = V(reordered_dof);
    trial_velocity = current_velocity + incre_velocity;
    update_incre_status(D);
}

/**
 * \brief Method to update acceleration and call previous method to further update velocity and displacement. This is used in Dynamic analysis.
 * \param D `incre_displacement`
 * \param V `incre_velocity`
 * \param A `incre_acceleration`
 */
void Node::update_incre_status(const vec& D, const vec& V, const vec& A) {
    incre_acceleration = A(reordered_dof);
    trial_acceleration = current_acceleration + incre_acceleration;
    update_incre_status(D, V);
}

/**
 * \brief Method to update displacement.
 * \param D `trial_displacement`
 */
void Node::update_trial_status(const vec& D) {
    trial_displacement = D(reordered_dof);
    incre_displacement = trial_displacement - current_displacement;
}

/**
 * \brief Method to update velocity and call previous method to further update displacement. This is used in Dynamic analysis.
 * \param D `trial_displacement`
 * \param V `trial_velocity`
 */
void Node::update_trial_status(const vec& D, const vec& V) {
    trial_velocity = V(reordered_dof);
    incre_velocity = trial_velocity - current_velocity;
    update_trial_status(D);
}

/**
 * \brief Method to update acceleration and call previous method to further update velocity and displacement. This is used in Dynamic analysis.
 * \param D `trial_displacement`
 * \param V `trial_velocity`
 * \param A `trial_acceleration`
 */
void Node::update_trial_status(const vec& D, const vec& V, const vec& A) {
    trial_acceleration = A(reordered_dof);
    incre_acceleration = trial_acceleration - current_acceleration;
    update_trial_status(D, V);
}

/**
 * \brief Method to commit the status variables.
 */
void Node::commit_status() {
    if(!trial_resistance.is_empty()) {
        current_resistance = trial_resistance;
        incre_resistance.zeros();
    }
    if(!trial_displacement.is_empty()) {
        current_displacement = trial_displacement;
        incre_displacement.zeros();
    }
    if(!trial_velocity.is_empty()) {
        current_velocity = trial_velocity;
        incre_displacement.zeros();
    }
    if(!trial_acceleration.is_empty()) {
        current_acceleration = trial_acceleration;
        incre_acceleration.zeros();
    }
}

/**
 * \brief Method to reset the status, there is no need to call this method as those variables can be directly overwritten.
 */
void Node::reset_status() {
    if(!current_resistance.is_empty()) {
        trial_resistance = current_resistance;
        incre_resistance.zeros();
    }
    if(!current_displacement.is_empty()) {
        trial_displacement = current_displacement;
        incre_displacement.zeros();
    }
    if(!current_velocity.is_empty()) {
        trial_velocity = current_velocity;
        incre_velocity.zeros();
    }
    if(!current_acceleration.is_empty()) {
        trial_acceleration = current_acceleration;
        incre_acceleration.zeros();
    }
}

/**
 * \brief The method tests each status variable before filling it by zeros. For any of them, empty means it is not used in analysis, then just keeps it unchanged.
 */
void Node::clear_status() {
    if(!current_resistance.is_empty()) {
        current_resistance.zeros();
        incre_resistance.zeros();
        trial_resistance.zeros();
    }
    if(!current_displacement.is_empty()) {
        current_displacement.zeros();
        incre_displacement.zeros();
        trial_displacement.zeros();
    }
    if(!current_velocity.is_empty()) {
        current_velocity.zeros();
        incre_velocity.zeros();
        trial_velocity.zeros();
    }
    if(!current_acceleration.is_empty()) {
        current_acceleration.zeros();
        incre_acceleration.zeros();
        trial_acceleration.zeros();
    }
}

vector<vec> Node::record(const OutputType& L) const {
    vector<vec> data;

    switch(L) {
    case OutputType::RF:
        data.push_back(current_resistance);
        break;
    case OutputType::U:
        data.push_back(current_displacement);
        break;
    case OutputType::V:
        data.push_back(current_velocity);
        break;
    case OutputType::A:
        data.push_back(current_acceleration);
        break;
    case OutputType::U1:
        if(current_displacement.n_elem >= 1) data.emplace_back(std::initializer_list<double>{ current_displacement(0) });
        break;
    case OutputType::U2:
        if(current_displacement.n_elem >= 2) data.emplace_back(std::initializer_list<double>{ current_displacement(1) });
        break;
    case OutputType::U3:
        if(current_displacement.n_elem >= 3) data.emplace_back(std::initializer_list<double>{ current_displacement(2) });
        break;
    case OutputType::RF1:
        if(current_resistance.n_elem >= 1) data.emplace_back(std::initializer_list<double>{ current_resistance(0) });
        break;
    case OutputType::RF2:
        if(current_resistance.n_elem >= 2) data.emplace_back(std::initializer_list<double>{ current_resistance(1) });
        break;
    case OutputType::RF3:
        if(current_resistance.n_elem >= 3) data.emplace_back(std::initializer_list<double>{ current_resistance(2) });
        break;
    default:
        break;
    }

    return data;
}

/**
 * \brief Method to print basic information.
 */
void Node::print() {
    suanpan_info("Node %u:\n", get_tag());
    coordinate.t().print();
    suanpan_info("Displacement:\n");
    current_displacement.t().print();
    suanpan_info("Resistance:\n");
    current_resistance.t().print();
    if(accu(current_velocity) != 0.) {
        suanpan_info("Velocity:\n");
        current_velocity.t().print();
    }
    if(accu(current_acceleration) != 0.) {
        suanpan_info("Acceleration:\n");
        current_acceleration.t().print();
    }
}
