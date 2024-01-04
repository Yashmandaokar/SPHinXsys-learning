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
 * @file 	base_general_dynamics.h
 * @brief 	This is the particle dynamics applicable for all type bodies
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_GENERAL_DYNAMICS_H
#define BASE_GENERAL_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "base_body.h"
#include "base_particles.hpp"
#include "external_force.h"

#include <limits>

namespace SPH
{
typedef DataDelegateSimple<BaseParticles> GeneralDataDelegateSimple;
typedef DataDelegateInner<BaseParticles> GeneralDataDelegateInner;
typedef DataDelegateContact<BaseParticles, BaseParticles> GeneralDataDelegateContact;
typedef DataDelegateContact<BaseParticles, BaseParticles, DataDelegateEmptyBase> GeneralDataDelegateContactOnly;
} // namespace SPH
#endif // BASE_GENERAL_DYNAMICS_H
