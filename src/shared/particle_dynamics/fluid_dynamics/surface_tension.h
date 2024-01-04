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
 * @file surface_tension.h
 * @brief A momentum-conservative formulation for surface tension is used here
 * to reach a long-term stable simulation.
 * @details TBD.
 * @author	Shuaihao Zhang and Xiangyu Hu
 */

#ifndef SURFACE_TENSION_H
#define SURFACE_TENSION_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
class SurfaceTensionStress : public LocalDynamics, public FluidContactData
{
  public:
    explicit SurfaceTensionStress(BaseContactRelation &contact_relation, StdVec<Real> contact_surface_tension);
    virtual ~SurfaceTensionStress(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> color_gradient_;
    StdLargeVec<Matd> surface_tension_stress_;
    StdVec<Real> contact_surface_tension_, contact_fraction_;
};

template <typename... T>
class SurfaceStressAcceleration;

template <>
class SurfaceStressAcceleration<Inner<>> : public LocalDynamics, public FluidDataInner
{
  public:
    SurfaceStressAcceleration(BaseInnerRelation &inner_relation);
    virtual ~SurfaceStressAcceleration(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &rho_, &mass_;
    StdLargeVec<Vecd> &force_prior_;
    StdLargeVec<Vecd> &color_gradient_;
    StdLargeVec<Matd> &surface_tension_stress_;
};

template <>
class SurfaceStressAcceleration<Contact<>> : public LocalDynamics, public FluidContactData
{
  public:
    explicit SurfaceStressAcceleration(BaseContactRelation &contact_relation);
    virtual ~SurfaceStressAcceleration(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &rho_, &mass_;
    StdLargeVec<Vecd> &force_prior_;
    StdLargeVec<Vecd> &color_gradient_;
    StdLargeVec<Matd> &surface_tension_stress_;
    StdVec<StdLargeVec<Vecd> *> contact_color_gradient_;
    StdVec<StdLargeVec<Matd> *> contact_surface_tension_stress_;
    StdVec<Real> contact_surface_tension_, contact_fraction_;
};

using SurfaceStressAccelerationComplex = ComplexInteraction<SurfaceStressAcceleration<Inner<>, Contact<>>>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // SURFACE_TENSION_H
