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

#ifndef ARPACK_WRAPPER_H
#define ARPACK_WRAPPER_H

#include <Domain/MetaMat/MetaMat>
#include <armadillo>
#include <memory>

using namespace arma;

// GENERAL MATRIX
int eig_solve(cx_vec&, cx_mat&, mat&, unsigned, const char* = "SM");
int eig_solve(cx_vec&, cx_mat&, mat&, mat&, unsigned, const char* = "SM");
int eig_solve(vec&, mat&, mat&, unsigned, const char* = "SM");
int eig_solve(vec&, mat&, mat&, mat&, unsigned, const char* = "SM");
int eig_solve(vec&, mat&, std::shared_ptr<BandSymmMat<double>>&, std::shared_ptr<BandSymmMat<double>>&, unsigned, const char* = "SM");
int eig_solve(vec&, mat&, std::shared_ptr<SymmPackMat<double>>&, std::shared_ptr<SymmPackMat<double>>&, unsigned, const char* = "SM");

#endif
