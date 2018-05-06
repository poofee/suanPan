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

#ifndef CSR_FORM_HPP
#define CSR_FORM_HPP

#include "sparse_form.hpp"
#include <cstring>

template <typename T> class triplet_form;
template <typename T> class csc_form;

template <typename T> class csr_form : public sparse_form<T, csr_form<T>> {
    void copy_memory(uword, const uword*, const uword*, const T*) override;

public:
    using sparse_form<T, csr_form<T>>::n_rows;
    using sparse_form<T, csr_form<T>>::n_cols;
    using sparse_form<T, csr_form<T>>::n_elem;
    using sparse_form<T, csr_form<T>>::c_size;

    uword* row_ptr = nullptr; // index storage
    uword* col_idx = nullptr; // index storage
    T* val_idx = nullptr;     // value storage

    csr_form() = default;

    ~csr_form();

    csr_form(const csr_form&);                // copy ctor
    csr_form(csr_form&&) noexcept;            // move ctor
    csr_form& operator=(const csr_form&);     // copy assignment
    csr_form& operator=(csr_form&&) noexcept; // move assignment

    void reset() const override;
    void zeros() const override;

    bool init() override;
    bool init(uword) override;
    bool init(uword, uword, uword) override;
    bool resize() override;
    bool resize(uword) override;
    bool resize(uword, uword, uword) override;

    void print() const override;

    csr_form<T> transpose() const;

    template <typename T2> csr_form<T> operator*(T2);
    template <typename T2> csr_form<T> operator/(T2);
    template <typename T2> csr_form<T>& operator*=(T2);
    template <typename T2> csr_form<T>& operator/=(T2);

    explicit csr_form(triplet_form<T>&);
    explicit csr_form(triplet_form<T>&&);
    csr_form& operator=(triplet_form<T>&);

    explicit csr_form(const csc_form<T>&);
    csr_form& operator=(const csc_form<T>&);

    T operator()(uword, uword);
};

template <typename T> void csr_form<T>::copy_memory(const uword size, const uword* const in_row_ptr, const uword* const in_col_idx, const T* const in_val_idx) {
    if(size > n_elem) resize(size);

    auto bytes = (n_rows + 1) * sizeof(uword);
    memcpy(this->row_ptr, in_row_ptr, bytes);
    bytes = size * sizeof(uword);
    memcpy(this->col_idx, in_col_idx, bytes);
    bytes = size * sizeof(T);
    memcpy(this->val_idx, in_val_idx, bytes);

    c_size = size;
}

template <typename T> csr_form<T>::~csr_form() { csr_form<T>::reset(); }

template <typename T>
csr_form<T>::csr_form(const csr_form& in_mat)
    : sparse_form<T, csr_form<T>>(in_mat.n_rows, in_mat.n_cols, in_mat.n_elem) {
    csr_form<T>::init();
    csr_form<T>::copy_memory(in_mat.c_size, in_mat.row_ptr, in_mat.col_idx, in_mat.val_idx);
}

template <typename T> csr_form<T>::csr_form(csr_form&& in_mat) noexcept {
    csr_form<T>::reset();
    n_rows = in_mat.n_rows;
    n_cols = in_mat.n_cols;
    n_elem = in_mat.n_elem;
    c_size = in_mat.c_size;
    row_ptr = in_mat.row_ptr;
    col_idx = in_mat.col_idx;
    val_idx = in_mat.val_idx;
    in_mat.n_rows = in_mat.n_cols = in_mat.n_elem = in_mat.c_size = 0;
    in_mat.row_ptr = in_mat.col_idx = nullptr;
    in_mat.val_idx = nullptr;
}

template <typename T> csr_form<T>& csr_form<T>::operator=(const csr_form& in_mat) {
    if(this != &in_mat) {
        init(in_mat.n_rows, in_mat.n_cols, in_mat.n_elem);
        copy_memory(in_mat.c_size, in_mat.col_idx, in_mat.row_ptr, in_mat.val_idx);
    }

    return *this;
}

template <typename T> csr_form<T>& csr_form<T>::operator=(csr_form&& in_mat) noexcept {
    reset();
    n_rows = in_mat.n_rows;
    n_cols = in_mat.n_cols;
    n_elem = in_mat.n_elem;
    c_size = in_mat.c_size;
    col_idx = in_mat.col_idx;
    row_ptr = in_mat.row_ptr;
    val_idx = in_mat.val_idx;
    in_mat.n_rows = in_mat.n_cols = in_mat.n_elem = in_mat.c_size = 0;
    in_mat.row_ptr = in_mat.col_idx = nullptr;
    in_mat.val_idx = nullptr;
    return *this;
}

template <typename T> void csr_form<T>::reset() const {
    zeros();
    delete[] col_idx;
    delete[] row_ptr;
    delete[] val_idx;
}

template <typename T> void csr_form<T>::zeros() const { access::rw(c_size) = 0; }

template <typename T> bool csr_form<T>::init() {
    reset();
    row_ptr = new(std::nothrow) uword[n_rows + 1];
    col_idx = new(std::nothrow) uword[n_elem];
    val_idx = new(std::nothrow) T[n_elem];

    if(col_idx == nullptr || row_ptr == nullptr || val_idx == nullptr) {
        reset();
        return false;
    }
    return true;
}

template <typename T> bool csr_form<T>::init(const uword in_elem) {
    if(in_elem <= n_elem) {
        zeros();
        return true;
    }
    n_elem = in_elem;
    return init();
}

template <typename T> bool csr_form<T>::init(const uword in_row, const uword in_col, const uword in_elem) {
    if(n_rows != in_row) n_rows = in_row;
    if(n_cols != in_col) n_cols = in_col;

    return init(in_elem);
}

template <typename T> bool csr_form<T>::resize() {
    const auto copy = *this;

    if(!init(n_elem == 0 ? 1 : 2 * n_elem)) return false;

    copy_memory(copy.c_size, copy.row_ptr, copy.col_idx, copy.val_idx);

    return true;
}

template <typename T> bool csr_form<T>::resize(const uword in_elem) {
    const auto copy = *this;

    if(in_elem <= c_size || !init(in_elem)) return false;

    copy_memory(copy.c_size, copy.row_ptr, copy.col_idx, copy.val_idx);

    return true;
}

template <typename T> bool csr_form<T>::resize(const uword in_row, const uword in_col, const uword in_elem) {
    const auto copy = *this;

    if(in_row < n_rows || in_col < n_cols || in_elem < c_size || !init(in_row, in_col, in_elem)) return false;

    copy_memory(copy.c_size, copy.row_ptr, copy.col_idx, copy.val_idx);

    return true;
}

template <typename T> void csr_form<T>::print() const {
    suanpan_info("A sparse matrix in triplet form with size of %u by %u, the sparsity of %.3f.\n", unsigned(n_rows), unsigned(n_cols), 100. - double(c_size) / double(n_rows * n_cols) * 100.);
    if(c_size > 1000) {
        suanpan_info("more than 1000 elements exist.\n");
        return;
    }

    uword c_idx = 1;
    for(uword I = 0; I < c_size; ++I) {
        if(I >= row_ptr[c_idx]) ++c_idx;
        suanpan_info("(%3u, %3u) ===> %+.4E\n", unsigned(c_idx) - 1, unsigned(col_idx[I]), val_idx[I]);
    }
}

template <typename T> csr_form<T> csr_form<T>::transpose() const {
    csc_form<T> copy(*this);
    csr_form<T> out;

    out.n_rows = n_cols;
    out.n_cols = n_rows;
    out.n_elem = n_elem;
    out.c_size = c_size;

    out.row_ptr = copy.col_ptr;
    out.col_idx = copy.row_idx;
    out.val_idx = copy.val_idx;

    copy.row_idx = copy.col_ptr = nullptr;
    copy.val_idx = nullptr;

    return out;
}

template <typename T> template <typename T2> csr_form<T> csr_form<T>::operator*(const T2 scalar) {
    csr_form<T> copy = *this;

    for(auto I = 0; I < copy.c_size; ++I) copy.val_idx[I] *= T(scalar);

    return copy;
}

template <typename T> template <typename T2> csr_form<T> csr_form<T>::operator/(const T2 scalar) {
    csr_form<T> copy = *this;

    for(auto I = 0; I < copy.c_size; ++I) copy.val_idx[I] /= T(scalar);

    return copy;
}

template <typename T> template <typename T2> csr_form<T>& csr_form<T>::operator*=(const T2 scalar) {
    for(auto I = 0; I < c_size; ++I) val_idx[I] *= T(scalar);

    return *this;
}

template <typename T> template <typename T2> csr_form<T>& csr_form<T>::operator/=(const T2 scalar) {
    for(auto I = 0; I < c_size; ++I) val_idx[I] /= T(scalar);

    return *this;
}

template <typename T> csr_form<T>::csr_form(triplet_form<T>& old_mat) { *this = old_mat; }

template <typename T> csr_form<T>::csr_form(triplet_form<T>&& old_mat) { *this = old_mat; }

template <typename T> csr_form<T>& csr_form<T>::operator=(triplet_form<T>& in_mat) {
    init(in_mat.n_rows, in_mat.n_cols, in_mat.c_size);
    in_mat.csr_sort();

    // set column pointers to zero
    for(uword I = 0; I <= n_rows; ++I) row_ptr[I] = 0;

    auto last_row = in_mat.row(0), last_col = in_mat.col(0);

    uword current_pos = 0, pre_row = last_row;

    auto last_sum = 0.;

    for(uword I = 0; I < in_mat.c_size; ++I) {
        if(in_mat.row(I) != last_row || in_mat.col(I) != last_col) {
            // now all components of the first element have been added up
            if(last_sum != 0.) {
                val_idx[current_pos] = last_sum;
                col_idx[current_pos] = last_col;
                // if the committed column index does not equal to the previous one
                if(pre_row != last_row) row_ptr[pre_row = last_row] = current_pos;
                last_sum = 0.;
                ++current_pos;
            }
            // check in the position of next potential element if the first element is zero
            last_row = in_mat.row(I);
            last_col = in_mat.col(I);
        }
        last_sum += in_mat.val(I);
    }

    // check in the last element
    if(last_sum != 0.) {
        val_idx[current_pos] = last_sum;
        col_idx[current_pos] = last_col;
        row_ptr[n_rows] = ++current_pos;
    } else
        row_ptr[n_rows] = current_pos;

    c_size = current_pos;

    return *this;
}

template <typename T> csr_form<T>::csr_form(const csc_form<T>& in_mat) { *this = in_mat; }

template <typename T> csr_form<T>& csr_form<T>::operator=(const csc_form<T>& in_mat) {
    csr_form<T>::init(in_mat.n_rows, in_mat.n_cols, in_mat.c_size);

    for(uword I = 0; I <= n_rows; ++I) row_ptr[I] = 0;
    for(uword I = 0; I < in_mat.c_size; ++I) ++row_ptr[in_mat.row_idx[I] + 1];
    for(uword I = 2; I <= n_rows; ++I) row_ptr[I] += row_ptr[I - 1];

    std::vector<uword> counter(n_rows, 0);

    uword c_idx = 1;
    for(uword I = 0; I < in_mat.c_size; ++I) {
        if(I >= in_mat.col_ptr[c_idx]) ++c_idx;
        const auto& r_idx = in_mat.row_idx[I];
        const auto c_pos = counter[r_idx]++ + row_ptr[r_idx];
        col_idx[c_pos] = c_idx - 1;
        val_idx[c_pos] = in_mat.val_idx[I];
    }

    c_size = in_mat.c_size;

    return *this;
}

template <typename T> T csr_form<T>::operator()(const uword in_row, const uword in_col) {
    if(in_row < n_rows && in_col < n_cols)
        for(auto I = row_ptr[in_row - 1]; I < row_ptr[in_row]; ++I)
            if(col_idx[I] == in_col) return val_idx[I];

    return 0.;
}

#endif
