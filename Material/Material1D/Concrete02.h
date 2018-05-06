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
 * @class Concrete02
 * @brief A Concrete02 material class.
 *
 * The Concrete02 class defines a uniaxial model for a collection of various concrete models.
 *
 * For tension branch, the response is always assumed to be origin--oriented, that is, the unloading curve always passes through the origin.
 *
 * An exponential decay curve is used for tension backbone.
 *
 * The Concrete02 class provides four types of compression backbones: THORENFELDT, POPOVICS, TSAI, KPSC and KPSU.
 *
 * The response could either be center oriented, no residual strain, or obay a KJ unloading equation.
 *
 * Caveat: The models used in this class maybe unit dependent. Only SI system is implemented here.
 *
 * @author tlc
 * @date 02/11/2017
 * @version 0.2.1
 * @file Concrete02.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CONCRETE02_H
#define CONCRETE02_H

#include <Material/Material1D/Material1D.h>

class Concrete02 : public Material1D {
    enum class Status { NONE, CBACKBONE, TBACKBONE, CUNLOAD, TUNLOAD, CRELOAD, TRELOAD, CTRANS, TTRANS };

    static const double one_over_six;

    Status load_status = Status::NONE;

    const double peak_stress, peak_strain;

    const double tension_origin;

    double crack_stress = 0.;

    double MC = 0., NC = 0., MT = 0., NT = 0.;

    podarray<bool> trial_flag, current_flag;

    podarray<double> compute_compression_backbone(double) const;
    podarray<double> compute_tension_backbone(double) const;

    static podarray<double> compute_transition(double, double, double, double, double, double, double);

public:
    static const double crack_strain;

    Concrete02(unsigned, // tag
        double,          // peak stress in negative
        double = 0.      // density
    );

    void initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    double get_parameter(const ParameterType&) const override;

    int update_trial_status(const vec&) override;
    int update_trial_status(const vec&, const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
