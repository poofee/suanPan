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

#ifndef OPERATOR_TIMES_HPP
#define OPERATOR_TIMES_HPP

template <typename T> MetaMat<T> operator*(const T& value, const MetaMat<T>& M) {
    auto N = M;
    arrayops::inplace_mul(N.memptr(), value, N.n_elem);
    return N;
}

template <typename T> Mat<T> operator*(const Mat<T>& A, const FullMat<T>& B) {
    Mat<T> C(A.n_rows, A.n_cols);

    const auto TRAN = 'N';

    auto M = static_cast<int>(A.n_rows);
    auto N = static_cast<int>(B.n_cols);
    auto K = static_cast<int>(A.n_cols);
    T ALPHA = 1.0;
    auto LDA = M;
    auto LDB = K;
    T BETA = 0.;
    auto LDC = M;

    if(std::is_same<T, float>::value) {
        using E = float;
        arma_fortran(arma_sgemm)(&TRAN, &TRAN, &M, &N, &K, (E*)&ALPHA, (E*)A.memptr(), &LDA, (E*)B.memptr(), &LDB, (E*)&BETA, (E*)C.memptr(), &LDC);
    } else if(std::is_same<T, double>::value) {
        using E = double;
        arma_fortran(arma_dgemm)(&TRAN, &TRAN, &M, &N, &K, (E*)&ALPHA, (E*)A.memptr(), &LDA, (E*)B.memptr(), &LDB, (E*)&BETA, (E*)C.memptr(), &LDC);
    }

    return C;
}

template <const char S, const char T, typename T1> Mat<T1> spmm(const SymmPackMat<T1>& A, const Mat<T1>& B) {
    Mat<T1> C;

    auto SIDE = S;
    auto TRAN = T;
    auto UPLO = 'U';

    auto M = static_cast<int>(A.n_rows);

    auto PT = 0;
    if(SIDE == 'L') PT += 1;
    if(TRAN == 'T') PT += 10;

    int N, LDC;

    switch(PT) {
    case 0: // A*B
        N = static_cast<int>(B.n_cols);
        C.set_size(M, N);
        LDC = M;
        break;
    case 1: // B*A
        N = static_cast<int>(B.n_rows);
        C.set_size(N, M);
        LDC = N;
        break;
    case 10: // A*B**T
        N = static_cast<int>(B.n_rows);
        C.set_size(M, N);
        LDC = M;
        break;
    case 11: // B**T*A
        N = static_cast<int>(B.n_cols);
        C.set_size(N, M);
        LDC = N;
        break;
    default:
        break;
    }

    T1 ALPHA = 1.;
    auto LDB = static_cast<int>(B.n_rows);
    T1 BETA = 0.;

    if(std::is_same<T1, float>::value) {
        using E = float;
        arma_fortran(arma_sspmm)(&SIDE, &UPLO, &TRAN, &M, &N, (E*)A.memptr(), (E*)&ALPHA, (E*)B.memptr(), &LDB, (E*)&BETA, (E*)C.memptr(), &LDC);
    } else if(std::is_same<T1, double>::value) {
        using E = double;
        arma_fortran(arma_dspmm)(&SIDE, &UPLO, &TRAN, &M, &N, (E*)A.memptr(), (E*)&ALPHA, (E*)B.memptr(), &LDB, (E*)&BETA, (E*)C.memptr(), &LDC);
    }

    return C;
}

template <typename T> Mat<T> operator*(const Mat<T>& A, const SymmPackMat<T>& B) { return spmm<'L', 'N'>(B, A); }

template <typename T> Mat<T> operator*(const Op<Mat<T>, op_htrans>& A, const SymmPackMat<T>& B) { return spmm<'L', 'T'>(B, A.m); }

#endif // OPERATOR_TIMES_HPP
