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
 * @class Modulated
 * @brief An Amplitude class that can generate Amplitude pattern.
 *
 * Linear/Modulated
 * \f{gather}{a=A_0+A\left(t-t_0\right)/t_s\f}
 *
 * @author tlc
 * @date 26/10/2017
 * @version 0.1.0
 * @file Modulated.h
 * @addtogroup Amplitude
 * @{
 */

#ifndef MODULATED_H
#define MODULATED_H

#include <Load/Amplitude/Amplitude.h>

using std::vector;

class Modulated : public Amplitude {
    const double amp;
    const vector<double> freq;

public:
    explicit Modulated(unsigned, // tag
        double,                  // amplitude
        const vector<double>&,   // frequency
        unsigned = 0);

    double get_amplitude(double) override;

    void print() override final;
};

#endif

//! @}
