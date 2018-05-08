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
 * @class SparseSymmMat
 * @brief A SparseSymmMat class that holds matrices.
 *
 * @author tlc
 * @date 06/05/2018
 * @version 0.1.0
 * @file SparseSymmMat.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef SPARSSYMMEMAT_HPP
#define SPARSSYMMEMAT_HPP

template <typename T> class SparseSymmMat : public SparseMat<T> {
    superlu_opts symm_opts;

public:
    using SparseMat<T>::SparseMat;

    Mat<T> solve(const Mat<T>&) override;
    int solve(Mat<T>&, const Mat<T>&) override;
};

template <typename T> Mat<T> SparseSymmMat<T>::solve(const Mat<T>& in_mat) {
    this->csc_to_arma();
    symm_opts.symmetric = true;
    return spsolve(this->arma_mat, in_mat, "s", symm_opts);
}

template <typename T> int SparseSymmMat<T>::solve(Mat<T>& out_mat, const Mat<T>& in_mat) {
    this->csc_to_arma();
    symm_opts.symmetric = true;
    return spsolve(out_mat, this->arma_mat, in_mat, "s", symm_opts) ? 0 : -1;
}

#endif

//! @}
