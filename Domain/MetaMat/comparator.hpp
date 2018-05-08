﻿/*******************************************************************************
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

#ifndef COMPARATOR_HPP
#define COMPARATOR_HPP

#include <armadillo>

using namespace arma;

class sparse_comparator {
protected:
    const uword* row_idx;
    const uword* col_idx;

public:
    sparse_comparator(const uword* const in_row_idx, const uword* const in_col_idx)
        : row_idx(in_row_idx)
        , col_idx(in_col_idx) {}
    virtual ~sparse_comparator() = default;
    virtual bool operator()(uword, uword) const = 0;
};

class csr_comparator : public sparse_comparator {
public:
    using sparse_comparator::sparse_comparator;
    bool operator()(const uword idx_a, const uword idx_b) const override {
        if(row_idx[idx_a] < row_idx[idx_b]) return true;
        if(row_idx[idx_a] > row_idx[idx_b]) return false;
        if(col_idx[idx_a] < col_idx[idx_b]) return true;
        return false;
    }
};

class csc_comparator : public sparse_comparator {
public:
    using sparse_comparator::sparse_comparator;
    bool operator()(const uword idx_a, const uword idx_b) const override {
        if(col_idx[idx_a] < col_idx[idx_b]) return true;
        if(col_idx[idx_a] > col_idx[idx_b]) return false;
        if(row_idx[idx_a] < row_idx[idx_b]) return true;
        return false;
    }
};

class abs_comparator : public sparse_comparator {
public:
    explicit abs_comparator(const uword* const location)
        : sparse_comparator(location, nullptr) {}
    bool operator()(const uword idx_a, const uword idx_b) const override { return row_idx[idx_a] < row_idx[idx_b]; }
};

#endif
