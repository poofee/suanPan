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
 * @class Concrete01
 * @brief A Concrete01 material class.
 *
 * The Concrete01 class defines a uniaxial model for a collection of various concrete models.
 *
 * For tension branch, the response is always assumed to be origin--oriented, that is, the unloading curve always passes through the origin.
 *
 * An exponential decay curve is used for tension backbone.
 *
 * The Concrete01 class provides four types of compression backbones: THORENFELDT, POPOVICS, TSAI and KPS.
 *
 * The response could either be center oriented, no residual strain, or obay a KJ unloading equation.
 *
 * Caveat: The models used in this class maybe unit dependent. Only SI system is implemented here.
 *
 * @author tlc
 * @date 02/11/2017
 * @version 0.2.1
 * @file Concrete01.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CONCRETE01_H
#define CONCRETE01_H

#include <Material/Material1D/Material1D.h>

class Concrete01 : public Material1D {
    static const double one_over_six;

    const BackboneType backbone_type;
    const TensionType tension_type;

    const bool no_residual, no_tension;
    const double peak_stress, peak_strain;
    const double facture_energy, tension_origin, confinement;

    double crack_stress, crack_strain = 0., ultimate_strain = 0.;
    double M = 0., N = 0.;

    podarray<bool> trial_flag, current_flag;

    void compute_compression_backbone(double);
    void compute_tension_backbone(double);

public:
    Concrete01(unsigned,                   // tag
        double,                            // peak stress in negative
        BackboneType = BackboneType::TSAI, // backbone type
        bool = false,                      // center oriented or using unloading criterion
        bool = false,                      // if define tension response (or zero for tension)
        TensionType = TensionType::LINEAR, // tensin softening law
        double = 1E-2,                     // fracture energy
        double = 0.                        // density
    );

    void initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    double get_parameter(const ParameterType&) const override;

    bool is_cracked();

    int update_trial_status(const vec&) override;
    int update_trial_status(const vec&, const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
