/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020,2021, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_GPU_UTILS_GPUTRAITS_H
#define GMX_GPU_UTILS_GPUTRAITS_H

/*! \libinternal \file
 *  \brief Declares the GPU type traits for non-GPU builds.
 *
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_gpu_utils
 */

#include "config.h"

#if GMX_GPU_CUDA

#    include "gromacs/gpu_utils/gputraits.cuh"

#elif GMX_GPU_OPENCL

#    include "gromacs/gpu_utils/gputraits_ocl.h"

#elif GMX_GPU_SYCL

#    include "gromacs/gpu_utils/gputraits_sycl.h"

#else

using DeviceTexture = void*;

//! \brief Single GPU call timing event
using CommandEvent = void*;

// Stubs for CPU-only build. Might be changed in #3312.
struct Float2
{
};
struct Float3
{
};
struct Float4
{
};

#endif // GMX_GPU

namespace gmx
{
template<typename T>
static inline Float3* asGenericFloat3Pointer(T* in)
{
    static_assert(sizeof(T) == sizeof(Float3),
                  "Size of the host-side data-type is different from the size of the generic "
                  "device-side counterpart.");
    return reinterpret_cast<Float3*>(in);
}

template<typename T>
static inline const Float3* asGenericFloat3Pointer(const T* in)
{
    static_assert(sizeof(T) == sizeof(Float3),
                  "Size of the host-side data-type is different from the size of the generic "
                  "device-side counterpart.");
    return reinterpret_cast<const Float3*>(in);
}

template<typename C>
static inline Float3* asGenericFloat3Pointer(C& in)
{
    static_assert(sizeof(*in.data()) == sizeof(Float3),
                  "Size of the host-side data-type is different from the size of the device-side "
                  "counterpart.");
    return reinterpret_cast<Float3*>(in.data());
}

template<typename C>
static inline const Float3* asGenericFloat3Pointer(const C& in)
{
    static_assert(sizeof(*in.data()) == sizeof(Float3),
                  "Size of the host-side data-type is different from the size of the device-side "
                  "counterpart.");
    return reinterpret_cast<const Float3*>(in.data());
}
} // namespace gmx

#endif // GMX_GPU_UTILS_GPUTRAITS_H
