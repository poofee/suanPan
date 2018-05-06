/*
    Copyright (c) 2005-2018 Intel Corporation

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.




*/

#ifndef __TBB_compat_ppl_H
#define __TBB_compat_ppl_H

#include "../combinable.h"
#include "../critical_section.h"
#include "../parallel_for.h"
#include "../parallel_for_each.h"
#include "../parallel_invoke.h"
#include "../reader_writer_lock.h"
#include "../task_group.h"
#include "../tbb_exception.h"

namespace Concurrency {

#if __TBB_TASK_GROUP_CONTEXT
using tbb::invalid_multiple_scheduling;
using tbb::make_task;
using tbb::missing_wait;
using tbb::structured_task_group;
using tbb::task_group;
using tbb::task_group_status;
using tbb::task_handle;

using tbb::canceled;
using tbb::complete;
using tbb::not_complete;

using tbb::is_current_task_group_canceling;
#endif /* __TBB_TASK_GROUP_CONTEXT */

using tbb::combinable;
using tbb::critical_section;
using tbb::parallel_for_each;
using tbb::parallel_invoke;
using tbb::reader_writer_lock;
using tbb::strict_ppl::parallel_for;

using tbb::improper_lock;

} // namespace Concurrency

#endif /* __TBB_compat_ppl_H */
