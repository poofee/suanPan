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

#ifndef SPARSE_FORM_HPP
#define SPARSE_FORM_HPP

#include <Toolbox/debug.h>
#include <armadillo>
#include <vector>

using namespace arma;

template <typename T, typename FT> class sparse_form {
protected:
    const T bin = 0.; // bin for out of bound elements

    virtual void copy_memory(uword, const uword*, const uword*, const T*) = 0;

public:
    typedef FT form_type;

    const uword n_rows = 0; // number of rows
    const uword n_cols = 0; // number of cols
    const uword n_elem = 0; // maximum number of elements
    const uword c_size = 0; // current number of valid elements

    sparse_form() = default;
    sparse_form(uword, uword, uword = 0);

    virtual ~sparse_form() = default;

    sparse_form(const sparse_form&) = delete;            // copy ctor
    sparse_form(sparse_form&&) = delete;                 // move ctor
    sparse_form& operator=(const sparse_form&) = delete; // copy assignment
    sparse_form& operator=(sparse_form&&) = delete;      // move assignment

    virtual void reset() const = 0;
    virtual void zeros() const = 0;

    virtual bool init() = 0;
    virtual bool init(uword) = 0;
    virtual bool init(uword, uword, uword) = 0;
    virtual bool resize() = 0;
    virtual bool resize(uword) = 0;
    virtual bool resize(uword, uword, uword) = 0;

    virtual void print() const;
    virtual void spy();
};

template <typename T, typename FT>
sparse_form<T, FT>::sparse_form(const uword in_rows, const uword in_cols, const uword in_elem)
    : n_rows(in_rows)
    , n_cols(in_cols)
    , n_elem(in_elem) {}

template <typename T, typename FT> void sparse_form<T, FT>::print() const {}

template <typename T, typename FT> void sparse_form<T, FT>::spy() { throw; }

#endif
