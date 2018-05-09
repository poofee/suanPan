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
 * @class Concrete22
 * @brief A Concrete22 material class.
 *
 * unsigned tag -> element unique tag
 * double peak_strain -> (-) the corresponding compressive strain at maximum compressive stress
 * double peak_stress -> (-) the maximum compressive stress
 * BackboneType backbone_type -> the compressive backbone type (THORENFELDT, POPOVICS, TSAI, etc.)
 * double beta -> (0,1) the shear retention factor
 * bool centre_oriented -> (TRUE,FALSE) whether to use centre oriented compressive unloading rule
 * double gf -> normalized mode one specific fracture energy, equals to the area under tension softening curve
 * double density -> (+) density
 * PlaneType plane_type -> whether it is a plane stress or plane strain problem
 *
 * @author tlc
 * @date 04/04/2018
 * @version 0.1.2
 * @file Concrete22.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef CONCRETE22_H
#define CONCRETE22_H

#include <Material/Material1D/Concrete01.h>
#include <Material/Material2D/Material2D.h>

class Concrete22 : public Material2D {
    Concrete01 concrete_major, concrete_minor;

    const bool poisson, degrade;

    const double shear_retention, peak_strain;

    const double poissons_ratio = poisson ? .2 : .0;

    double elastic_modulus = 0., shear_modulus = 0.;

    mat poissons_mat;

public:
    Concrete22(unsigned,                   // tag
        double,                            // peak stress in negative
        BackboneType = BackboneType::TSAI, // backbone type
        double = .2,                       // shear retention
        bool = false,                      // center oriented or using unloading criterion
        TensionType = TensionType::LINEAR, // tension softening type
        double = 1E-2,                     // factrue energy
        bool = false,                      // poisson effect switch
        bool = false,                      // stiffness degradation
        double = 0.,                       // density
        PlaneType = PlaneType::S);

    void initialize(const shared_ptr<DomainBase>& = nullptr) override;

    unique_ptr<Material> get_copy() override;

    double get_parameter(const ParameterType&) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
