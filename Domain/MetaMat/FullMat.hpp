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
 * @class FullMat
 * @brief A FullMat class that holds matrices.
 *
 * @author tlc
 * @date 06/09/2017
 * @version 0.1.0
 * @file FullMat.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef FULLMAT_HPP
#define FULLMAT_HPP

template <typename T> class FullMat : public MetaMat<T> {
public:
    FullMat();
    explicit FullMat(uword);
};

template <typename T> struct is_Full { static const bool value = false; };

template <typename T> struct is_Full<FullMat<T>> { static const bool value = true; };

template <typename T>
FullMat<T>::FullMat()
    : MetaMat<T>() {}

template <typename T>
FullMat<T>::FullMat(const uword in_size)
    : MetaMat<T>(in_size, in_size, in_size * in_size) {}

#endif

//! @}
