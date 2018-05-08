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
#ifndef TRIPLET_FORM
#define TRIPLET_FORM

#include "comparator.hpp"
#include "sparse_form.hpp"
#include <cstring>

template <typename T> class csc_form;
template <typename T> class csr_form;

template <typename T> class triplet_form : public sparse_form<T, triplet_form<T>> {
    using sparse_form<T, triplet_form<T>>::bin;

    void condense() const;

    void copy_memory(uword, const uword*, const uword*, const T*) override;

public:
    using sparse_form<T, triplet_form<T>>::n_rows;
    using sparse_form<T, triplet_form<T>>::n_cols;
    using sparse_form<T, triplet_form<T>>::n_elem;
    using sparse_form<T, triplet_form<T>>::c_size;

    uword* row_idx = nullptr; // index storage
    uword* col_idx = nullptr; // index storage
    T* val_idx = nullptr;     // value storage

    bool csc_sorted = false;
    bool csr_sorted = false;

    triplet_form() = default;
    triplet_form(uword, uword, uword = 0, bool = false, bool = false);

    ~triplet_form();

    triplet_form(const triplet_form&);                // copy ctor
    triplet_form(triplet_form&&) noexcept;            // move ctor
    triplet_form& operator=(const triplet_form&);     // copy assignment
    triplet_form& operator=(triplet_form&&) noexcept; // move assignment

    void reset() const override;
    void zeros() const override;

    T max() const override;

    bool init() override;
    bool init(uword) override;
    bool init(uword, uword, uword) override;
    bool resize() override;
    bool resize(uword) override;
    bool resize(uword, uword, uword) override;

    void print() const override;
    void spy() override;

    triplet_form<T> transpose() const;

    void csr_condense() const;
    void csc_condense() const;
    bool csr_sort() const;
    bool csc_sort() const;

    uword row(uword) const;
    uword col(uword) const;
    T val(uword) const;

    const T& operator()(uword, uword) const;
    T& at(uword, uword);

    template <typename T2> triplet_form<T> operator*(T2);
    template <typename T2> triplet_form<T> operator/(T2);
    template <typename T2> triplet_form<T>& operator*=(T2);
    template <typename T2> triplet_form<T>& operator/=(T2);

    triplet_form<T> operator+(const triplet_form<T>&);
    triplet_form<T> operator-(const triplet_form<T>&);
    triplet_form<T>& operator+=(const triplet_form<T>&);
    triplet_form<T>& operator-=(const triplet_form<T>&);

    explicit triplet_form(const csc_form<T>&);
    triplet_form& operator=(const csc_form<T>&);
    triplet_form<T> operator+(const csc_form<T>&);
    triplet_form<T> operator-(const csc_form<T>&);
    triplet_form<T>& operator+=(const csc_form<T>&);
    triplet_form<T>& operator-=(const csc_form<T>&);

    explicit triplet_form(const csr_form<T>&);
    triplet_form& operator=(const csr_form<T>&);
    triplet_form<T> operator+(const csr_form<T>&);
    triplet_form<T> operator-(const csr_form<T>&);
    triplet_form<T>& operator+=(const csr_form<T>&);
    triplet_form<T>& operator-=(const csr_form<T>&);
};

template <typename T> void triplet_form<T>::condense() const {
    auto last_row = row_idx[0], last_col = col_idx[0];

    uword current_pos = 0;
    auto last_sum = 0.;

    uword max_row = 0, max_col = 0;

    for(uword I = 0; I < c_size; ++I) {
        if(last_row != row_idx[I] || last_col != col_idx[I]) {
            if(last_sum != 0.) {
                row_idx[current_pos] = last_row;
                col_idx[current_pos] = last_col;
                val_idx[current_pos] = last_sum;
                if(last_row > max_row) max_row = last_row;
                if(last_col > max_col) max_col = last_col;
                ++current_pos;
                last_sum = 0.;
            }
            last_row = row_idx[I];
            last_col = col_idx[I];
        }
        last_sum += val_idx[I];
    }

    if(last_sum != 0.) {
        row_idx[current_pos] = last_row;
        col_idx[current_pos] = last_col;
        val_idx[current_pos] = last_sum;
        if(last_row > max_row) max_row = last_row;
        if(last_col > max_col) max_col = last_col;
        ++current_pos;
    }

    access::rw(n_rows) = max_row + 1;
    access::rw(n_cols) = max_col + 1;
    access::rw(c_size) = current_pos;
}

template <typename T> void triplet_form<T>::copy_memory(const uword size, const uword* const in_row_idx, const uword* const in_col_idx, const T* const in_val_idx) {
    const auto new_size = this->c_size + size;

    if(new_size > n_elem) resize(new_size);

    auto bytes = size * sizeof(uword);
    memcpy(this->row_idx + this->c_size, in_row_idx, bytes);
    memcpy(this->col_idx + this->c_size, in_col_idx, bytes);

    bytes = size * sizeof(T);
    memcpy(this->val_idx + this->c_size, in_val_idx, bytes);

    access::rw(c_size) = new_size;
}

template <typename T>
triplet_form<T>::triplet_form(const uword in_rows, const uword in_cols, const uword in_elem, const bool in_csc_sort, const bool in_csr_sort)
    : sparse_form<T, triplet_form<T>>(in_rows, in_cols, in_elem) {
    csc_sorted = in_csc_sort;
    csr_sorted = in_csr_sort;
    triplet_form<T>::init();
}

template <typename T> triplet_form<T>::~triplet_form() { triplet_form<T>::reset(); }

template <typename T>
triplet_form<T>::triplet_form(const triplet_form& in_mat)
    : sparse_form<T, triplet_form<T>>(in_mat.n_rows, in_mat.n_cols, in_mat.n_elem) {
    csc_sorted = in_mat.csc_sorted;
    csr_sorted = in_mat.csr_sorted;
    triplet_form<T>::init();
    triplet_form<T>::copy_memory(in_mat.c_size, in_mat.row_idx, in_mat.col_idx, in_mat.val_idx);
}

template <typename T> triplet_form<T>::triplet_form(triplet_form&& in_mat) noexcept {
    triplet_form<T>::reset();
    n_rows = in_mat.n_rows;
    n_cols = in_mat.n_cols;
    n_elem = in_mat.n_elem;
    c_size = in_mat.c_size;
    csc_sorted = in_mat.csc_sorted;
    csr_sorted = in_mat.csr_sorted;
    row_idx = in_mat.row_idx;
    col_idx = in_mat.col_idx;
    val_idx = in_mat.val_idx;
    in_mat.n_rows = in_mat.n_cols = in_mat.n_elem = in_mat.c_size = 0;
    in_mat.row_idx = in_mat.col_idx = nullptr;
    in_mat.val_idx = nullptr;
}

template <typename T> triplet_form<T>& triplet_form<T>::operator=(const triplet_form& in_mat) {
    if(this != &in_mat) {
        init(in_mat.n_rows, in_mat.n_cols, in_mat.n_elem, in_mat.sorted);
        copy_memory(in_mat.c_size, in_mat.row_idx, in_mat.col_idx, in_mat.val_idx);
    }

    return *this;
}

template <typename T> triplet_form<T>& triplet_form<T>::operator=(triplet_form&& in_mat) noexcept {
    reset();
    n_rows = in_mat.n_rows;
    n_cols = in_mat.n_cols;
    n_elem = in_mat.n_elem;
    c_size = in_mat.c_size;
    csc_sorted = in_mat.csc_sorted;
    csr_sorted = in_mat.csr_sorted;
    row_idx = in_mat.row_idx;
    col_idx = in_mat.col_idx;
    val_idx = in_mat.val_idx;
    in_mat.n_rows = in_mat.n_cols = in_mat.n_elem = in_mat.c_size = 0;
    in_mat.row_idx = in_mat.col_idx = nullptr;
    in_mat.val_idx = nullptr;
    return *this;
}

template <typename T> void triplet_form<T>::reset() const {
    zeros();
    delete[] row_idx;
    delete[] col_idx;
    delete[] val_idx;
}

template <typename T> void triplet_form<T>::zeros() const {
    access::rw(c_size) = 0;
    access::rw(csc_sorted) = false;
    access::rw(csr_sorted) = false;
}

template <typename T> T triplet_form<T>::max() const { return *std::max_element(val_idx, val_idx + c_size); }

template <typename T> bool triplet_form<T>::init() {
    reset();
    if(n_elem == 0) return true;
    row_idx = new(std::nothrow) uword[n_elem];
    col_idx = new(std::nothrow) uword[n_elem];
    val_idx = new(std::nothrow) T[n_elem];
    if(row_idx == nullptr || col_idx == nullptr || val_idx == nullptr) {
        reset();
        return false;
    }
    return true;
}

template <typename T> bool triplet_form<T>::init(const uword in_elem) {
    if(in_elem <= n_elem) {
        zeros();
        return true;
    }
    access::rw(n_elem) = in_elem;
    return init();
}

template <typename T> bool triplet_form<T>::init(const uword in_row, const uword in_col, const uword in_elem) {
    if(n_rows != in_row) access::rw(n_rows) = in_row;
    if(n_cols != in_col) access::rw(n_cols) = in_col;

    return init(in_elem);
}

template <typename T> bool triplet_form<T>::resize() {
    const auto copy = *this;

    if(!init(n_elem == 0 ? 1 : 2 * n_elem)) return false;

    copy_memory(copy.c_size, copy.row_idx, copy.col_idx, copy.val_idx);

    return true;
}

template <typename T> bool triplet_form<T>::resize(const uword in_elem) {
    const auto copy = *this;

    if(in_elem <= c_size || !init(in_elem)) return false;

    copy_memory(copy.c_size, copy.row_idx, copy.col_idx, copy.val_idx);

    return true;
}

template <typename T> bool triplet_form<T>::resize(const uword in_row, const uword in_col, const uword in_elem) {
    const auto copy = *this;

    if(in_row < n_rows || in_col < n_cols || in_elem < c_size || !init(in_row, in_col, in_elem)) return false;

    copy_memory(copy.c_size, copy.row_idx, copy.col_idx, copy.val_idx);

    return true;
}

template <typename T> void triplet_form<T>::print() const {
    suanpan_info("A sparse matrix in triplet form with size of %u by %u, the sparsity of %.3f.\n", unsigned(n_rows), unsigned(n_cols), double(c_size) / double(n_rows * n_cols) * 100.);
    if(c_size > 1000) {
        suanpan_info("Not going to print all elements as more than 1000 elements exist.\n");
        return;
    }
    for(uword I = 0; I < c_size; ++I) suanpan_info("(%3u, %3u) ===> %+.4E\n", unsigned(row_idx[I]), unsigned(col_idx[I]), val_idx[I]);
}

template <typename T> void triplet_form<T>::spy() {
    if(n_rows > 100 || n_cols > 100) return;

    csr_condense();

    uword current_pos = 0;

    for(uword I = 0; I < n_rows; ++I) {
        for(uword J = 0; J < n_cols; ++J)
            if(I == row_idx[current_pos] && J == col_idx[current_pos]) {
                suanpan_info(" X");
                ++current_pos;
            } else {
                suanpan_info(" .");
            }
        suanpan_info("\n");
    }
}

/**
 * \brief simply swap index arrays
 */
template <typename T> triplet_form<T> triplet_form<T>::transpose() const {
    auto copy = *this;
    auto t_ptr = copy.row_idx;
    copy.row_idx = copy.col_idx;
    copy.col_idx = t_ptr;
    auto t_num = copy.n_rows;
    copy.n_rows = copy.n_cols;
    copy.n_cols = t_num;
    return copy;
}

template <typename T> void triplet_form<T>::csr_condense() const {
    csr_sort();
    condense();
}

template <typename T> void triplet_form<T>::csc_condense() const {
    csc_sort();
    condense();
}

template <typename T> bool triplet_form<T>::csr_sort() const {
    if(c_size < 2) return true;

    const auto new_row_idx = new(std::nothrow) uword[n_elem];
    const auto new_col_idx = new(std::nothrow) uword[n_elem];
    const auto new_val_idx = new(std::nothrow) T[n_elem];

    if(new_row_idx == nullptr || new_col_idx == nullptr || new_val_idx == nullptr) {
        delete[] new_row_idx;
        delete[] new_col_idx;
        delete[] new_val_idx;
        return false;
    }

    std::vector<uword> index(c_size);
    for(uword I = 0; I < uword(index.size()); ++I) {
        new_row_idx[I] = row_idx[I] * n_cols + col_idx[I];
        index[I] = I;
    }

    std::sort(index.begin(), index.end(), abs_comparator(new_row_idx));

    for(uword I = 0; I < c_size; ++I) {
        new_row_idx[I] = row_idx[index[I]];
        new_col_idx[I] = col_idx[index[I]];
        new_val_idx[I] = val_idx[index[I]];
    }

    reset();

    access::rwp(row_idx) = new_row_idx;
    access::rwp(col_idx) = new_col_idx;
    access::rwp(val_idx) = new_val_idx;
    access::rw(c_size) = uword(index.size());

    access::rw(csr_sorted) = true;

    return true;
}

template <typename T> bool triplet_form<T>::csc_sort() const {
    if(c_size < 2) return true;

    const auto new_row_idx = new(std::nothrow) uword[n_elem];
    const auto new_col_idx = new(std::nothrow) uword[n_elem];
    const auto new_val_idx = new(std::nothrow) T[n_elem];

    if(new_row_idx == nullptr || new_col_idx == nullptr || new_val_idx == nullptr) {
        delete[] new_row_idx;
        delete[] new_col_idx;
        delete[] new_val_idx;
        return false;
    }

    std::vector<uword> index(c_size);
    for(uword I = 0; I < index.size(); ++I) {
        new_row_idx[I] = col_idx[I] * n_rows + row_idx[I];
        index[I] = I;
    }

    std::sort(index.begin(), index.end(), abs_comparator(new_row_idx));

    for(uword I = 0; I < c_size; ++I) {
        new_row_idx[I] = row_idx[index[I]];
        new_col_idx[I] = col_idx[index[I]];
        new_val_idx[I] = val_idx[index[I]];
    }

    reset();

    access::rwp(row_idx) = new_row_idx;
    access::rwp(col_idx) = new_col_idx;
    access::rwp(val_idx) = new_val_idx;
    access::rw(c_size) = uword(index.size());

    access::rw(csc_sorted) = true;

    return true;
}

template <typename T> uword triplet_form<T>::row(const uword idx) const {
    if(idx < c_size) return row_idx[idx];
    throw;
}

template <typename T> uword triplet_form<T>::col(const uword idx) const {
    if(idx < c_size) return col_idx[idx];
    throw;
}

template <typename T> T triplet_form<T>::val(const uword idx) const {
    if(idx < c_size) return val_idx[idx];
    throw;
}

template <typename T> const T& triplet_form<T>::operator()(const uword row, const uword col) const {
    for(uword I = 0; I < c_size; ++I)
        if(row_idx[I] == row && col_idx[I] == col) return val_idx[I];

    access::rw(bin) = 0.;

    return bin;
}

template <typename T> T& triplet_form<T>::at(const uword row, const uword col) {
    if(row < n_rows && col < n_cols) {
        if(csr_sorted) csr_sorted = false;
        if(csc_sorted) csc_sorted = false;
        if(c_size == n_elem) resize();
        row_idx[c_size] = row;
        col_idx[c_size] = col;
        return val_idx[access::rw(c_size)++];
    }

    return access::rw(bin);
}

template <typename T> template <typename T2> triplet_form<T> triplet_form<T>::operator*(const T2 scalar) {
    triplet_form<T> copy = *this;

    for(auto I = 0; I < copy.c_size; ++I) copy.val_idx[I] *= T(scalar);

    return copy;
}

template <typename T> template <typename T2> triplet_form<T> triplet_form<T>::operator/(const T2 scalar) {
    triplet_form<T> copy = *this;

    for(auto I = 0; I < copy.c_size; ++I) copy.val_idx[I] /= T(scalar);

    return copy;
}

template <typename T> template <typename T2> triplet_form<T>& triplet_form<T>::operator*=(const T2 scalar) {
    for(auto I = 0; I < c_size; ++I) val_idx[I] *= T(scalar);

    return *this;
}

template <typename T> template <typename T2> triplet_form<T>& triplet_form<T>::operator/=(const T2 scalar) {
    for(auto I = 0; I < c_size; ++I) val_idx[I] /= T(scalar);

    return *this;
}

template <typename T> triplet_form<T> triplet_form<T>::operator+(const triplet_form<T>& in_mat) {
    triplet_form<T> copy = *this;
    return copy += in_mat;
}

template <typename T> triplet_form<T> triplet_form<T>::operator-(const triplet_form<T>& in_mat) {
    triplet_form<T> copy = *this;
    return copy -= in_mat;
}

template <typename T> triplet_form<T>& triplet_form<T>::operator+=(const triplet_form<T>& in_mat) {
    auto new_size = c_size + in_mat.c_size;
    if(n_elem < new_size) resize(new_size);

    copy_memory(in_mat.c_size, in_mat.row_idx, in_mat.col_idx, in_mat.val_idx);

    if(in_mat.n_rows > n_rows) n_rows = in_mat.n_rows;
    if(in_mat.n_cols > n_cols) n_cols = in_mat.n_cols;

    return *this;
}

template <typename T> triplet_form<T>& triplet_form<T>::operator-=(const triplet_form<T>& in_mat) {
    auto new_size = c_size + in_mat.c_size;
    if(n_elem < new_size) resize(new_size);

    for(auto I = 0; I < in_mat.c_size; ++I) at(in_mat.row(I), in_mat.col(I)) = -in_mat.val(I);

    if(in_mat.n_rows > n_rows) n_rows = in_mat.n_rows;
    if(in_mat.n_cols > n_cols) n_cols = in_mat.n_cols;

    return *this;
}

template <typename T> triplet_form<T>::triplet_form(const csc_form<T>& in_mat) { *this = in_mat; }

template <typename T> triplet_form<T>& triplet_form<T>::operator=(const csc_form<T>& in_mat) {
    init(in_mat.n_rows, in_mat.n_cols, in_mat.c_size);

    c_size = in_mat.c_size;

    uword c_idx = 1;
    for(uword I = 0; I < in_mat.c_size; ++I) {
        if(I >= in_mat.col_ptr[c_idx]) ++c_idx;
        row_idx[I] = in_mat.row_idx[I];
        col_idx[I] = c_idx - 1;
        val_idx[I] = in_mat.val_idx[I];
    }

    return *this;
}

template <typename T> triplet_form<T> triplet_form<T>::operator+(const csc_form<T>& in_mat) {
    triplet_form<T> copy(in_mat);
    return copy += *this;
}

template <typename T> triplet_form<T> triplet_form<T>::operator-(const csc_form<T>& in_mat) {
    triplet_form<T> copy(in_mat);
    return copy -= *this;
}

template <typename T> triplet_form<T>& triplet_form<T>::operator+=(const csc_form<T>& in_mat) { return *this += triplet_form<T>(in_mat); }

template <typename T> triplet_form<T>& triplet_form<T>::operator-=(const csc_form<T>& in_mat) { return *this -= triplet_form<T>(in_mat); }

template <typename T> triplet_form<T>::triplet_form(const csr_form<T>& in_mat) { *this = in_mat; }

template <typename T> triplet_form<T>& triplet_form<T>::operator=(const csr_form<T>& in_mat) {
    init(in_mat.n_rows, in_mat.n_cols, in_mat.c_size);

    c_size = in_mat.c_size;

    uword c_idx = 1;
    for(uword I = 0; I < in_mat.c_size; ++I) {
        if(I >= in_mat.row_ptr[c_idx]) ++c_idx;
        row_idx[I] = c_idx - 1;
        col_idx[I] = in_mat.col_idx[I];
        val_idx[I] = in_mat.val_idx[I];
    }

    return *this;
}

template <typename T> triplet_form<T> triplet_form<T>::operator+(const csr_form<T>& in_mat) {
    triplet_form<T> copy(in_mat);
    return copy += *this;
}

template <typename T> triplet_form<T> triplet_form<T>::operator-(const csr_form<T>& in_mat) {
    triplet_form<T> copy(in_mat);
    return copy -= *this;
}

template <typename T> triplet_form<T>& triplet_form<T>::operator+=(const csr_form<T>& in_mat) { return *this += triplet_form<T>(in_mat); }

template <typename T> triplet_form<T>& triplet_form<T>::operator-=(const csr_form<T>& in_mat) { return *this -= triplet_form<T>(in_mat); }

template <typename T> triplet_form<T> operator+(const triplet_form<T>& mat_a, const triplet_form<T>& mat_b) {
    auto out = mat_a;
    out += mat_b;
    return out;
}

#endif