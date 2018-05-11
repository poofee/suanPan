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
 * @class SparseMat
 * @brief A SparseMat class that holds matrices.
 *
 * @author tlc
 * @date 06/05/2018
 * @version 0.1.0
 * @file SparseMat.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef SPARSEMAT_HPP
#define SPARSEMAT_HPP

#include "csc_form.hpp"
#include "triplet_form.hpp"

template <typename T> class SparseMat : public MetaMat<T> {
protected:
    triplet_form<T> triplet_mat;
    csc_form<T> csc_mat;
    SpMat<T> arma_mat;

public:
    SparseMat() = default;
    explicit SparseMat(uword, uword);

    bool is_empty() const override;

    void zeros() override;
    void reset() override;
    T max() const override;

    const T& operator()(uword, uword) const override;
    T& at(uword, uword) override;

    virtual void triplet_to_csc();
    virtual void csc_to_arma();

    const T* memptr() const override { throw; }
    T* memptr() override { throw; }

    MetaMat<T> operator+(const MetaMat<T>&) override;
    MetaMat<T> operator-(const MetaMat<T>&) override;
    MetaMat<T>& operator+=(const MetaMat<T>&) override;
    MetaMat<T>& operator-=(const MetaMat<T>&) override;

    MetaMat<T> operator*(T)override;
    Mat<T> operator*(const Mat<T>&)override;
    MetaMat<T>& operator*=(T) override;

    Mat<T> solve(const Mat<T>&) override;
    int solve(Mat<T>&, const Mat<T>&) override;

    Mat<T> solve_trs(const Mat<T>&) override;
    int solve_trs(Mat<T>&, const Mat<T>&) override;

    MetaMat<T> factorize() override { throw; }

    MetaMat<T> i() override { throw; }
    MetaMat<T> inv() override { throw; }

    void print() override;
    void save(const char*) override;
};

template <typename T>
SparseMat<T>::SparseMat(const uword in_row, const uword in_col)
    : triplet_mat(in_row, in_col) {}

template <typename T> bool SparseMat<T>::is_empty() const { return arma_mat.is_empty(); }

template <typename T> void SparseMat<T>::zeros() {
    triplet_mat.zeros();
    csc_mat.zeros();
    arma_mat.zeros();
}

template <typename T> void SparseMat<T>::reset() {
    triplet_mat.reset();
    csc_mat.reset();
    arma_mat.reset();
}

template <typename T> T SparseMat<T>::max() const { return csc_mat.max(); }

template <typename T> const T& SparseMat<T>::operator()(uword in_row, uword in_col) const { return csc_mat(in_row, in_col); }

template <typename T> T& SparseMat<T>::at(uword in_row, uword in_col) { return triplet_mat.at(in_row, in_col); }

template <typename T> void SparseMat<T>::triplet_to_csc() { csc_mat = triplet_mat; }

template <typename T> void SparseMat<T>::csc_to_arma() {
    triplet_to_csc();

    const uvec row_idx(csc_mat.row_idx, csc_mat.c_size, false, false);
    const uvec col_ptr(csc_mat.col_ptr, csc_mat.n_cols + 1, false, false);
    const Col<T> val_idx(csc_mat.val_idx, csc_mat.c_size, false, false);

    arma_mat = SpMat<T>(row_idx, col_ptr, val_idx, csc_mat.n_rows, csc_mat.n_cols);
}

template <typename T> MetaMat<T> SparseMat<T>::operator+(const MetaMat<T>& in_mat) {
    auto N = *this;
    N.arma_mat += dynamic_cast<const SparseMat<T>&>(in_mat).arma_mat;
    return N;
}

template <typename T> MetaMat<T> SparseMat<T>::operator-(const MetaMat<T>& in_mat) {
    auto N = *this;
    N.arma_mat -= dynamic_cast<const SparseMat<T>&>(in_mat).arma_mat;
    return N;
}

template <typename T> MetaMat<T>& SparseMat<T>::operator+=(const MetaMat<T>& in_mat) {
    arma_mat += dynamic_cast<const SparseMat<T>&>(in_mat).arma_mat;
    return this;
}

template <typename T> MetaMat<T>& SparseMat<T>::operator-=(const MetaMat<T>& in_mat) {
    arma_mat -= dynamic_cast<const SparseMat<T>&>(in_mat).arma_mat;
    return this;
}

template <typename T> MetaMat<T> SparseMat<T>::operator*(const T scalar) {
    auto N = *this;
    N.arma_mat *= scalar;
    return N;
}

template <typename T> Mat<T> SparseMat<T>::operator*(const Mat<T>& in_mat) { return arma_mat * in_mat; }

template <typename T> MetaMat<T>& SparseMat<T>::operator*=(const T scalar) {
    arma_mat *= scalar;
    return this;
}

template <typename T> Mat<T> SparseMat<T>::solve(const Mat<T>& in_mat) {
    Mat<T> out_mat;
    if(solve(out_mat, in_mat) != 0) out_mat.reset();
    return out_mat;
}

template <typename T> int SparseMat<T>::solve(Mat<T>& out_mat, const Mat<T>& in_mat) {
#ifdef SUANPAN_MAGMA
    csc_to_arma();
    return spsolve(out_mat, arma_mat, in_mat) ? 0 : -1;
#else
    csc_to_arma();
    return spsolve(out_mat, arma_mat, in_mat) ? 0 : -1;
#endif
}

template <typename T> Mat<T> SparseMat<T>::solve_trs(const Mat<T>& in_mat) { return solve(in_mat); }

template <typename T> int SparseMat<T>::solve_trs(Mat<T>& out_mat, const Mat<T>& in_mat) { return solve(out_mat, in_mat); }

template <typename T> void SparseMat<T>::print() { arma_mat.print(); }

template <typename T> void SparseMat<T>::save(const char* name) { arma_mat.save(name, coord_ascii); }

#endif

//! @}
