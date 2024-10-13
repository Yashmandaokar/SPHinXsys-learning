/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	particle_iterators_ck.h
 * @brief 	This is for the base functions for particle iterator.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_ITERATORS_CK_H
#define PARTICLE_ITERATORS_CK_H

#include "implementation.h"
#include "loop_range.h"

#include <numeric>

namespace SPH
{
using namespace execution;

template <class LocalDynamicsFunction>
void particle_for(const SequencedPolicy &seq,
                  LoopRangeCK<SequencedPolicy, SPHBody> &loop_range,
                  const LocalDynamicsFunction &local_dynamics_function)
{
    for (size_t i = 0; i < loop_range.LoopBound(); ++i)
        local_dynamics_function(i);
};

template <class LocalDynamicsFunction>
void particle_for(const ParallelPolicy &par,
                  LoopRangeCK<ParallelPolicy, SPHBody> &loop_range,
                  const LocalDynamicsFunction &local_dynamics_function)
{
    parallel_for(
        IndexRange(0, loop_range.LoopBound()),
        [&](const IndexRange &r)
        {
            for (size_t i = r.begin(); i < r.end(); ++i)
            {
                local_dynamics_function(i);
            }
        },
        ap);
};

template <class ReturnType, typename Operation, class LocalDynamicsFunction>
ReturnType particle_reduce(const SequencedPolicy &seq,
                           LoopRangeCK<SequencedPolicy, SPHBody> &loop_range,
                           ReturnType temp, Operation &&operation,
                           const LocalDynamicsFunction &local_dynamics_function)
{
    for (size_t i = 0; i < loop_range.LoopBound(); ++i)
    {
        temp = operation(temp, local_dynamics_function(i));
    }
    return temp;
}

template <class ReturnType, typename Operation, class LocalDynamicsFunction>
ReturnType particle_reduce(const ParallelPolicy &par,
                           LoopRangeCK<ParallelPolicy, SPHBody> &loop_range,
                           ReturnType temp, Operation &&operation,
                           const LocalDynamicsFunction &local_dynamics_function)
{
    return parallel_reduce(
        IndexRange(0, loop_range.LoopBound()),
        temp, [&](const IndexRange &r, ReturnType temp0) -> ReturnType
        {
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					temp0 = operation(temp0, local_dynamics_function(i));
				}
				return temp0; },
        [&](const ReturnType &x, const ReturnType &y) -> ReturnType
        {
            return operation(x, y);
        });
};

template <typename T, typename Op>
T exclusive_scan(const SequencedPolicy &seq_policy, T *first, T *d_first, UnsignedInt d_size, Op op)
{
    UnsignedInt scan_size = d_size - 1;
    std::exclusive_scan(first, first + d_size, d_first, T{0}, op);
    return d_first[scan_size];
}

template <typename T, typename Op>
T exclusive_scan(const ParallelPolicy &par_policy, T *first, T *d_first, UnsignedInt d_size, Op op)
{
    // Exclusive scan is the same as inclusive, but shifted by one
    UnsignedInt scan_size = d_size - 1;
    d_first[0] = T{0};
    using range_type = tbb::blocked_range<UnsignedInt>;
    tbb::parallel_scan(
        range_type(0, scan_size), d_first[0],
        [=](const range_type &r, T sum, bool is_final_scan) -> T
        {
            T tmp = sum;
            for (UnsignedInt i = r.begin(); i < r.end(); ++i)
            {
                tmp = op(tmp, first[i]);
                if (is_final_scan)
                {
                    d_first[i + 1] = tmp;
                }
            }
            return tmp;
        },
        [&](const T &a, const T &b)
        {
            return op(a, b);
        });
    return d_first[scan_size];
}
} // namespace SPH
#endif // PARTICLE_ITERATORS_CK_H
