////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2018 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "utility.h"
#include <cstring>
#include <suanPan.h>

char suanpan::to_upper(const char U) { return char(toupper(int(U))); }

char suanpan::to_lower(const char L) { return char(tolower(int(L))); }

bool is_equal(const char* A, const char* B) { return _strcmpi(A, B) == 0; }

bool is_equal(const char A, const char B) { return tolower(int(A)) == tolower(int(B)); }

bool is_equal(const string& A, const char* B) { return is_equal(A.c_str(), B); }

bool is_equal(const string& A, const string& B) { return is_equal(A.c_str(), B.c_str()); }

bool is_true(const char* S) { return is_equal(S, "On") || is_equal(S, "True") || is_equal(S, "T") || is_equal(S, "1") || is_equal(S, "Yes") || is_equal(S, "Y"); }

bool is_false(const char* S) { return is_equal(S, "Off") || is_equal(S, "False") || is_equal(S, "F") || is_equal(S, "0") || is_equal(S, "No") || is_equal(S, "N"); }

bool is_true(const string& S) { return is_true(S.c_str()); }

bool is_false(const string& S) { return is_false(S.c_str()); }
