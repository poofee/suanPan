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
 * @class Fibre2D
 * @brief A Fibre2D class.
 * @author tlc
 * @date 13/10/2017
 * @version 0.1.0
 * @file Fibre2D.h
 * @addtogroup Section-2D
 * @ingroup Section
 * @{
 */

#ifndef FIBRE2D_H
#define FIBRE2D_H

#include <Section/Section2D/Section2D.h>

class Fibre2D : public Section2D {
    podarray<unsigned> fibre_tag;

    vector<unique_ptr<Section>> fibre;

public:
    explicit Fibre2D(unsigned // tag
    );
    Fibre2D(const Fibre2D&);

    void initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Section> get_copy() override;

    double get_parameter(const ParameterType&) override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
