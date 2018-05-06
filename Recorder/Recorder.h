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
 * @class Recorder
 * @brief A Recorder class.
 * @author tlc
 * @date 27/07/2017
 * @version 0.1.0
 * @file Recorder.h
 * @{
 */

#ifndef RECORDER_H
#define RECORDER_H

#include <Domain/Tag.h>
#include <Recorder/OutputType.h>

class DomainBase;

using std::vector;

class Recorder : public Tag {
    uvec object_tag;
    OutputType variable_type;
    vector<double> time_pool;              /**< recorded data */
    vector<vector<vector<vec>>> data_pool; /**< recorded data */

    bool record_time;
    bool use_hdf5;

public:
    explicit Recorder(unsigned = 0, unsigned = CT_RECORDER, const uvec& = {}, const OutputType& = OutputType::NL, bool = true, bool = true);
    Recorder(const Recorder&) = delete;
    Recorder& operator=(const Recorder&) = delete;
    virtual ~Recorder();

    void set_object_tag(const uvec&);
    const uvec& get_object_tag() const;

    void set_variable_type(const OutputType&);
    const OutputType& get_variable_type() const;

    const bool& if_record_time() const;

    void insert(const double);
    void insert(const vector<vec>&, const unsigned);

    const vector<vector<vector<vec>>>& get_data_pool() const;
    const vector<double>& get_time_pool() const;

    virtual void record(const shared_ptr<DomainBase>&) = 0;

    virtual void save();

    void print() override;
};

#endif

//! @}
