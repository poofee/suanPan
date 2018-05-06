/**
 * @fn shapeFunction
 * @brief An shapeFunction function.
 *
 * The function is designed to return the shape function for four- and eight-node quad elements.
 *
 * @author T
 * @date 04/10/2017
 * @version 0.1.2
 * @file shapeFunction.hpp
 * @addtogroup Utility
 * @{
 */

#ifndef SHAPEFUNCTION_HPP
#define SHAPEFUNCTION_HPP

#include <armadillo>

using namespace arma;

namespace shape {
template <typename T> Mat<T> truss(const T INTPTS, const unsigned ORDER = 0, const unsigned NODENUM = 2) {
    Mat<T> N(1, NODENUM);

    const auto& X = INTPTS;

    if(NODENUM == 2) {
        if(ORDER == 0) {
            N(0, 0) = 1. - X;
            N(0, 1) = 1. + X;
        } else {
            N(0, 0) = -1.;
            N(0, 1) = 1.;
        }
        N /= 2.;
    } else if(NODENUM == 3) {
        if(ORDER == 0) {
            const auto XX = X * X;
            N(0, 0) = .5 * (XX - X);
            N(0, 1) = 1. - XX;
            N(0, 2) = .5 * (XX + X);
        } else {
            N(0, 0) = X - .5;
            N(0, 1) = -2. * X;
            N(0, 2) = X + .5;
        }
    }

    return N;
}

/**
 * \brief
 * \tparam T input argument type
 * \param INTPTS location of point
 * \param ORDER displacement (0) or rotation (1)
 * \param LENGTH beam length
 * \return vector of shape functions
 */
template <typename T> Col<T> beam(const T INTPTS, const unsigned ORDER, const double LENGTH) {
    Col<T> N(4);

    const auto XP = 1. + INTPTS;
    const auto XM = 1. - INTPTS;
    const auto XPP = XP * XP;
    const auto XMM = XM * XM;

    if(ORDER == 0) {
        N(0) = 2. * XMM * (XP + 1.);
        N(1) = LENGTH * XMM * XP;
        N(2) = 2. * XPP * (XM + 1.);
        N(3) = LENGTH * XM * XPP;
    } else if(ORDER == 1) {
        N(0) = -6. * XP * XM;
        N(1) = LENGTH * XM * (3. * INTPTS + 1.);
        N(2) = 6. * XP * XM;
        N(3) = LENGTH * XP * (3. * INTPTS - 1.);
    }

    N *= .125;

    return N;
}

template <typename T> Mat<T> triangle(const Col<T>& INTPTS, const unsigned ORDER) {
    Mat<T> N;

    if(ORDER != 0 && ORDER != 1) throw;

    N.zeros(ORDER + 1, 6);

    const auto& X = INTPTS(0);
    const auto& Y = INTPTS(1);

    if(ORDER == 0) {
        N(0, 0) = 1.;
        N(0, 1) = X;
        N(0, 2) = Y;
        N(0, 3) = X * Y;
        N(0, 4) = X * X;
        N(0, 5) = Y * Y;
    } else if(ORDER == 1) {
        N(0, 0) = 0.;
        N(0, 1) = 1.;
        N(0, 2) = 0;
        N(0, 3) = Y;
        N(0, 4) = 2. * X;
        N(0, 5) = 0.;
        N(1, 0) = 0.;
        N(1, 1) = 0.;
        N(1, 2) = 1.;
        N(1, 3) = X;
        N(1, 4) = 0.;
        N(1, 5) = 2. * Y;
    }

    return N;
}

template <typename T> Mat<T> quad(const Col<T>& INTPTS, const unsigned ORDER, const unsigned NODENUM = 4) {
    Mat<T> N;

    if(ORDER != 0 && ORDER != 1) throw;

    N.zeros(ORDER + 1, NODENUM);

    const auto& X = INTPTS(0);
    const auto& Y = INTPTS(1);

    const auto XP = 1. + X, XM = 1. - X, YP = 1. + Y, YM = 1. - Y;

    if(NODENUM == 4) {
        if(ORDER == 0) {
            N(0, 3) = XM * YP;
            N(0, 2) = XP * YP;
            N(0, 1) = XP * YM;
            N(0, 0) = XM * YM;
        } else if(ORDER == 1) {
            N(1, 1) = -(N(1, 2) = XP);
            N(1, 0) = -(N(1, 3) = XM);
            N(0, 3) = -(N(0, 2) = YP);
            N(0, 0) = -(N(0, 1) = YM);
        }
        N *= .25;
    } else if(NODENUM == 8) {
        const auto XX = X * X, YY = Y * Y, XY = X * Y;
        if(ORDER == 0) {
            N(0, 7) = .5 * XM * (1. - YY);
            N(0, 6) = .5 * (1. - XX) * YP;
            N(0, 5) = .5 * XP * (1. - YY);
            N(0, 4) = .5 * (1. - XX) * YM;
            N(0, 0) = .25 * XM * YM - .5 * (N(0, 4) + N(0, 7));
            N(0, 1) = .25 * XP * YM - .5 * (N(0, 4) + N(0, 5));
            N(0, 2) = .25 * XP * YP - .5 * (N(0, 5) + N(0, 6));
            N(0, 3) = .25 * XM * YP - .5 * (N(0, 6) + N(0, 7));
        } else if(ORDER == 1) {
            const auto X2 = .5 * X;
            const auto Y2 = .5 * Y;
            const auto X4 = .25 * X;
            const auto Y4 = .25 * Y;
            const auto X24 = .25 * XX;
            const auto Y24 = .25 * YY;
            const auto XY2 = .5 * XY;
            N(1, 7) = XY - Y;
            N(1, 6) = .5 - .5 * XX;
            N(1, 5) = -Y - XY;
            N(1, 4) = .5 * XX - .5;
            N(1, 3) = Y2 - X4 - XY2 + X24;
            N(1, 2) = X4 + Y2 + XY2 + X24;
            N(1, 1) = Y2 - X4 + XY2 - X24;
            N(1, 0) = X4 + Y2 - XY2 - X24;
            N(0, 7) = .5 * YY - .5;
            N(0, 6) = -X - XY;
            N(0, 5) = .5 - .5 * YY;
            N(0, 4) = XY - X;
            N(0, 3) = X2 - Y4 + XY2 - Y24;
            N(0, 2) = X2 + Y4 + XY2 + Y24;
            N(0, 1) = X2 - Y4 - XY2 + Y24;
            N(0, 0) = X2 + Y4 - XY2 - Y24;
        }
    }

    return N;
}

template <typename T> Mat<T> plate(const Col<T>& INTPTS, const unsigned ORDER, const unsigned NODENUM = 4) {
    Mat<T> N;

    if(ORDER == 0)
        N.zeros(1, 3 * NODENUM);
    else if(ORDER == 1)
        N.zeros(2, 6 * NODENUM);
    else
        throw;

    const auto& X = INTPTS(0);
    const auto& Y = INTPTS(1);

    return N;
}

template <typename T> Mat<T> cube(const Col<T>& INTPTS, const unsigned ORDER, const unsigned NODENUM = 8) {
    Mat<T> N;

    if(ORDER == 0)
        N.zeros(1, NODENUM);
    else if(ORDER == 1)
        N.zeros(3, NODENUM);
    else
        throw;

    const auto& X = INTPTS(0);
    const auto& Y = INTPTS(1);
    const auto& Z = INTPTS(2);

    const auto XP = X + 1.;
    const auto XM = X - 1.;
    const auto YP = Y + 1.;
    const auto YM = Y - 1.;
    const auto ZP = Z + 1.;
    const auto ZM = Z - 1.;

    if(NODENUM == 8) {
        if(ORDER == 0) {
            N(0, 0) = -XM * YM * ZM;
            N(0, 1) = XP * YM * ZM;
            N(0, 2) = -XP * YP * ZM;
            N(0, 3) = XM * YP * ZM;
            N(0, 4) = XM * YM * ZP;
            N(0, 5) = -XP * YM * ZP;
            N(0, 6) = XP * YP * ZP;
            N(0, 7) = -XM * YP * ZP;
        } else if(ORDER == 1) {
            N(0, 0) = -YM * ZM;
            N(0, 1) = YM * ZM;
            N(0, 2) = -YP * ZM;
            N(0, 3) = YP * ZM;
            N(0, 4) = YM * ZP;
            N(0, 5) = -YM * ZP;
            N(0, 6) = YP * ZP;
            N(0, 7) = -YP * ZP;
            N(1, 0) = -XM * ZM;
            N(1, 1) = XP * ZM;
            N(1, 2) = -XP * ZM;
            N(1, 3) = XM * ZM;
            N(1, 4) = XM * ZP;
            N(1, 5) = -XP * ZP;
            N(1, 6) = XP * ZP;
            N(1, 7) = -XM * ZP;
            N(2, 0) = -XM * YM;
            N(2, 1) = XP * YM;
            N(2, 2) = -XP * YP;
            N(2, 3) = XM * YP;
            N(2, 4) = XM * YM;
            N(2, 5) = -XP * YM;
            N(2, 6) = XP * YP;
            N(2, 7) = -XM * YP;
        }
        N *= .125;
    } else if(NODENUM == 20) {
        const auto XX = XP * XM;
        const auto YY = YP * YM;
        const auto ZZ = ZP * ZM;
        if(ORDER == 0) {
            N(0, 0) = XM * YM * ZM * (X + Y + Z + 2.) * .125;
            N(0, 1) = -XP * YM * ZM * (Y - X + Z + 2.) * .125;
            N(0, 2) = -XP * YP * ZM * (X + Y - Z - 2.) * .125;
            N(0, 3) = -XM * YP * ZM * (X - Y + Z + 2.) * .125;
            N(0, 4) = -XM * YM * ZP * (X + Y - Z + 2.) * .125;
            N(0, 5) = -XP * YM * ZP * (X - Y + Z - 2.) * .125;
            N(0, 6) = XP * YP * ZP * (X + Y + Z - 2.) * .125;
            N(0, 7) = XM * YP * ZP * (X - Y - Z + 2.) * .125;
            N(0, 8) = -XX * YM * ZM * .25;
            N(0, 9) = YY * XP * ZM * .25;
            N(0, 10) = XX * YP * ZM * .25;
            N(0, 11) = -YY * XM * ZM * .25;
            N(0, 12) = XX * YM * ZP * .25;
            N(0, 13) = -YY * XP * ZP * .25;
            N(0, 14) = -XX * YP * ZP * .25;
            N(0, 15) = YY * XM * ZP * .25;
            N(0, 16) = -ZZ * XM * YM * .25;
            N(0, 17) = ZZ * XP * YM * .25;
            N(0, 18) = -ZZ * XP * YP * .25;
            N(0, 19) = ZZ * XM * YP * .25;
        } else if(ORDER == 1) {
            N(0, 0) = XM * YM * ZM * .125 + YM * ZM * (X + Y + Z + 2.) * .125;
            N(0, 1) = XP * YM * ZM * .125 - YM * ZM * (Y - X + Z + 2.) * .125;
            N(0, 2) = -XP * YP * ZM * .125 - YP * ZM * (X + Y - Z - 2.) * .125;
            N(0, 3) = -XM * YP * ZM * .125 - YP * ZM * (X - Y + Z + 2.) * .125;
            N(0, 4) = -XM * YM * ZP * .125 - YM * ZP * (X + Y - Z + 2.) * .125;
            N(0, 5) = -XP * YM * ZP * .125 - YM * ZP * (X - Y + Z - 2.) * .125;
            N(0, 6) = XP * YP * ZP * .125 + YP * ZP * (X + Y + Z - 2.) * .125;
            N(0, 7) = XM * YP * ZP * .125 + YP * ZP * (X - Y - Z + 2.) * .125;
            N(0, 8) = -X * YM * ZM * .5;
            N(0, 9) = YY * ZM * .25;
            N(0, 10) = X * YP * ZM * .5;
            N(0, 11) = -YY * ZM * .25;
            N(0, 12) = X * YM * ZP * .5;
            N(0, 13) = -YY * ZP * .25;
            N(0, 14) = -X * YP * ZP * .5;
            N(0, 15) = YY * ZP * .25;
            N(0, 16) = -ZZ * YM * .25;
            N(0, 17) = ZZ * YM * .25;
            N(0, 18) = -ZZ * YP * .25;
            N(0, 19) = ZZ * YP * .25;
            N(1, 0) = XM * YM * ZM * .125 + XM * ZM * (X + Y + Z + 2.) * .125;
            N(1, 1) = -XP * YM * ZM * .125 - XP * ZM * (Y - X + Z + 2.) * .125;
            N(1, 2) = -XP * YP * ZM * .125 - XP * ZM * (X + Y - Z - 2.) * .125;
            N(1, 3) = XM * YP * ZM * .125 - XM * ZM * (X - Y + Z + 2.) * .125;
            N(1, 4) = -XM * YM * ZP * .125 - XM * ZP * (X + Y - Z + 2.) * .125;
            N(1, 5) = XP * YM * ZP * .125 - XP * ZP * (X - Y + Z - 2.) * .125;
            N(1, 6) = XP * YP * ZP * .125 + XP * ZP * (X + Y + Z - 2.) * .125;
            N(1, 7) = -XM * YP * ZP * .125 + XM * ZP * (X - Y - Z + 2.) * .125;
            N(1, 8) = -XX * ZM * .25;
            N(1, 9) = Y * XP * ZM * .5;
            N(1, 10) = XX * ZM * .25;
            N(1, 11) = -Y * XM * ZM * .5;
            N(1, 12) = XX * ZP * .25;
            N(1, 13) = -Y * XP * ZP * .5;
            N(1, 14) = -XX * ZP * .25;
            N(1, 15) = Y * XM * ZP * .5;
            N(1, 16) = -ZZ * XM * .25;
            N(1, 17) = ZZ * XP * .25;
            N(1, 18) = -ZZ * XP * .25;
            N(1, 19) = ZZ * XM * .25;
            N(2, 0) = XM * YM * ZM * .125 + XM * YM * (X + Y + Z + 2.) * .125;
            N(2, 1) = -XP * YM * ZM * .125 - XP * YM * (Y - X + Z + 2.) * .125;
            N(2, 2) = XP * YP * ZM * .125 - XP * YP * (X + Y - Z - 2.) * .125;
            N(2, 3) = -XM * YP * ZM * .125 - XM * YP * (X - Y + Z + 2.) * .125;
            N(2, 4) = XM * YM * ZP * .125 - XM * YM * (X + Y - Z + 2.) * .125;
            N(2, 5) = -XP * YM * ZP * .125 - XP * YM * (X - Y + Z - 2.) * .125;
            N(2, 6) = XP * YP * ZP * .125 + XP * YP * (X + Y + Z - 2.) * .125;
            N(2, 7) = -XM * YP * ZP * .125 + XM * YP * (X - Y - Z + 2.) * .125;
            N(2, 8) = -XX * YM * .25;
            N(2, 9) = YY * XP * .25;
            N(2, 10) = XX * YP * .25;
            N(2, 11) = -YY * XM * .25;
            N(2, 12) = XX * YM * .25;
            N(2, 13) = -YY * XP * .25;
            N(2, 14) = -XX * YP * .25;
            N(2, 15) = YY * XM * .25;
            N(2, 16) = -Z * XM * YM * .5;
            N(2, 17) = Z * XP * YM * .5;
            N(2, 18) = -Z * XP * YP * .5;
            N(2, 19) = Z * XM * YP * .5;
        }
    }

    return N;
}

template <typename T> Mat<T> stress(const T X, const T Y, const unsigned S) {
    Mat<T> N = zeros(3, S);

    for(auto I = 0; I < 3; ++I) N(I, I) = 1;

    if(S >= 5) {
        N(0, 4) = Y;
        N(1, 3) = X;
        if(S >= 7) {
            N(0, 6) = X;
            N(1, 5) = Y;
            N(2, 5) = -X;
            N(2, 6) = -Y;
            if(S >= 9) {
                const auto X2 = X * X;
                const auto Y2 = Y * Y;
                const auto XY = X * Y;
                N(0, 8) = 2. * XY;
                N(1, 7) = N(0, 8);
                N(2, 7) = -X2;
                N(2, 8) = -Y2;
                if(S == 11) {
                    N(0, 9) = -X2;
                    N(1, 9) = 2. * X2 - Y2;
                    N(2, 9) = 2. * XY;

                    N(0, 10) = 2. * Y2 - X2;
                    N(1, 10) = -Y2;
                    N(2, 10) = N(2, 9);
                }
            }
        }
    }

    return N;
}

template <typename T> Mat<T> stress(const Col<T>& C, const unsigned S) { return stress(C(0), C(1), S); }

template <typename T> Mat<T> stress5(const Col<T>& C) { return stress(C, 5); }
template <typename T> Mat<T> stress7(const Col<T>& C) { return stress(C, 7); }
template <typename T> Mat<T> stress9(const Col<T>& C) { return stress(C, 9); }
template <typename T> Mat<T> stress11(const Col<T>& C) { return stress(C, 11); }

template <typename T> Mat<T> stress5(const T X, const T Y) { return stress(X, Y, 5); }
template <typename T> Mat<T> stress7(const T X, const T Y) { return stress(X, Y, 7); }
template <typename T> Mat<T> stress9(const T X, const T Y) { return stress(X, Y, 9); }
template <typename T> Mat<T> stress11(const T X, const T Y) { return stress(X, Y, 11); }

inline mat stress5(const vec& C) { return stress(C, 5); }
inline mat stress7(const vec& C) { return stress(C, 7); }
inline mat stress9(const vec& C) { return stress(C, 9); }
inline mat stress11(const vec& C) { return stress(C, 11); }

template <typename T> Mat<T> strain(const T X, const T Y, const T V, const unsigned S) {
    Mat<T> N(3, S, fill::zeros);

    N(0, 0) = N(1, 1) = 1.;

    N(2, 2) = 2. + 2. * V;

    N(0, 1) = N(1, 0) = -V;

    if(S >= 5) {
        N(0, 3) = -V * (N(1, 3) = X);
        N(1, 4) = -V * (N(0, 4) = Y);
        if(S >= 7) {
            N(0, 5) = N(1, 4);
            N(0, 6) = N(1, 3);

            N(1, 5) = N(0, 4);
            N(1, 6) = N(0, 3);

            N(2, 5) = -X * N(2, 2);
            N(2, 6) = -Y * N(2, 2);
            if(S >= 9) {
                const auto X2 = X * X, Y2 = Y * Y, XY = X * Y;

                N(1, 8) = N(0, 7) = -V * (N(1, 7) = N(0, 8) = 2. * XY);

                N(2, 7) = -X2 * N(2, 2);
                N(2, 8) = -Y2 * N(2, 2);
                if(S == 11) {
                    N(0, 9) = V * Y2 - (2. * V + 1.) * X2;
                    N(1, 9) = (2. + V) * X2 - Y2;

                    N(0, 10) = (2. + V) * Y2 - X2;
                    N(1, 10) = V * X2 - (2. * V + 1.) * Y2;

                    N(2, 10) = N(2, 9) = 2. * XY * N(2, 2);
                }
            }
        }
    }

    return N;
}

template <typename T> Mat<T> strain(const Col<T>& C, const T V, const unsigned S) { return strain(C(0), C(1), V, S); }

template <typename T> Mat<T> strain5(const T X, const T Y, const T V) { return strain(X, Y, V, 5); }
template <typename T> Mat<T> strain7(const T X, const T Y, const T V) { return strain(X, Y, V, 7); }
template <typename T> Mat<T> strain9(const T X, const T Y, const T V) { return strain(X, Y, V, 9); }
template <typename T> Mat<T> strain11(const T X, const T Y, const T V) { return strain(X, Y, V, 11); }

template <typename T> Mat<T> strain5(const Col<T>& C, const T V) { return strain(C, V, 5); }
template <typename T> Mat<T> strain7(const Col<T>& C, const T V) { return strain(C, V, 7); }
template <typename T> Mat<T> strain9(const Col<T>& C, const T V) { return strain(C, V, 9); }
template <typename T> Mat<T> strain11(const Col<T>& C, const T V) { return strain(C, V, 11); }

inline mat strain5(const vec& C, const double V) { return strain(C, V, 5); }
inline mat strain7(const vec& C, const double V) { return strain(C, V, 7); }
inline mat strain9(const vec& C, const double V) { return strain(C, V, 9); }
inline mat strain11(const vec& C, const double V) { return strain(C, V, 11); }
}

#endif

//! @}
