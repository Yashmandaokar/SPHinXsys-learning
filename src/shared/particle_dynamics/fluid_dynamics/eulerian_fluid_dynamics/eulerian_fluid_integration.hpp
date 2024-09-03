#ifndef EULERIAN_FLUID_INTEGRATION_HPP
#define EULERIAN_FLUID_INTEGRATION_HPP

#include "eulerian_fluid_integration.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
EulerianIntegration<DataDelegationType>::EulerianIntegration(BaseRelationType &base_relation)
    : BaseIntegration<DataDelegationType>(base_relation),
      mom_(*this->particles_->template registerSharedVariable<Vecd>("Momentum")),
      dmom_dt_(*this->particles_->template registerSharedVariable<Vecd>("MomentumChangeRate")),
      dmass_dt_(*this->particles_->template registerSharedVariable<Real>("MassChangeRate")),
      Vol_(*this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      mom_advection_(*this->particles_->template registerSharedVariable<Vecd>("MomentumAdvection")),
      pressuregrad_(*this->particles_->template registerSharedVariable<Vecd>("MomentumPressureGradient")) {}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration1stHalf<Inner<>, RiemannSolverType>::
    EulerianIntegration1stHalf(BaseInnerRelation &inner_relation, Real limiter_parameter)
    : EulerianIntegration<DataDelegateInner>(inner_relation),
      riemann_solver_(this->fluid_, this->fluid_, limiter_parameter) {}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalf<Inner<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Vecd momentum_change_rate = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    mom_advection_[index_i] = Vecd::Zero(), pressuregrad_[index_i] = Vecd::Zero();
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];

        FluidStateIn state_j(rho_[index_j], vel_[index_j], p_[index_j]);
        FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, e_ij);
        Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();

        momentum_change_rate -= 2.0 * Vol_[index_i] * (convect_flux + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;

        mom_advection_[index_i] -= 2.0 * Vol_[index_i] * (convect_flux) * e_ij * dW_ijV_j;
        pressuregrad_[index_i] -= 2.0 * Vol_[index_i] * (interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
        if (index_i == 8269)
        {
            Vecd momentum_rate = -2.0 * Vol_[index_i] * (convect_flux + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
            Vecd Advflux = -2.0 * Vol_[index_i] * convect_flux * e_ij * dW_ijV_j;
            Vecd prgrad = -2.0 * Vol_[index_i] * interface_state.p_ * Matd::Identity() * e_ij * dW_ijV_j;
            Real x = 1;
        }
    }
    dmom_dt_[index_i] = momentum_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalf<Inner<>, RiemannSolverType>::update(size_t index_i, Real dt)
{
    mom_[index_i] += (dmom_dt_[index_i] + force_prior_[index_i]) * dt;
    //vel_[index_i] = mom_[index_i] / mass_[index_i];
    if (index_i == 8269)
    {
        Vecd visforce = force_prior_[index_i];
        Vecd mom = mom_[index_i];
        Real mass = mass_[index_i];
        Vecd veli = vel_[index_i];
        Real x = 1;
    }
}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration1stHalf<Contact<Wall>, RiemannSolverType>::
    EulerianIntegration1stHalf(BaseContactRelation &wall_contact_relation, Real limiter_parameter)
    : BaseEulerianIntegrationWithWall(wall_contact_relation),
      riemann_solver_(fluid_, fluid_, limiter_parameter) {}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalf<Contact<Wall>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Vecd momentum_change_rate = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &n_k = *(wall_n_[k]);
        StdLargeVec<Real> &Vol_k = *(wall_Vol_[k]);
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * Vol_k[index_j];

            Vecd vel_in_wall = -state_i.vel_;
            Real p_in_wall = state_i.p_;
            Real rho_in_wall = state_i.rho_;
            FluidStateIn state_j(rho_in_wall, vel_in_wall, p_in_wall);
            FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, n_k[index_j]);
            Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();
            momentum_change_rate -= 2.0 * Vol_[index_i] * (convect_flux + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
        }
    }
    dmom_dt_[index_i] += momentum_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration2ndHalf<Inner<>, RiemannSolverType>::
    EulerianIntegration2ndHalf(BaseInnerRelation &inner_relation, Real limiter_parameter)
    : EulerianIntegration<DataDelegateInner>(inner_relation),
      riemann_solver_(this->fluid_, this->fluid_, limiter_parameter)
      {}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration2ndHalf<Inner<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    //vel_[index_i] = vel_prof_[index_i];
    FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Real mass_change_rate = 0.0;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];

        //vel_[index_j] = vel_prof_[index_j];
        FluidStateIn state_j(rho_[index_j], vel_[index_j], p_[index_j]);
        FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, e_ij);
        mass_change_rate -= 2.0 * Vol_[index_i] * (interface_state.rho_ * interface_state.vel_).dot(e_ij) * dW_ijV_j;
        if (index_i == 8269)
        {
            Real voli = Vol_[index_i];
            Real volj = Vol_[index_j];
            Real c = 1.0;
        }

    }
    dmass_dt_[index_i] = mass_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration2ndHalf<Inner<>, RiemannSolverType>::update(size_t index_i, Real dt)
{
    mass_[index_i] += dmass_dt_[index_i] * dt;
    //rho_[index_i] = mass_[index_i] / Vol_[index_i];
    //p_[index_i] = fluid_.getPressure(rho_[index_i]);
    if (index_i == 8269)
    {
        Real mass = mass_[index_i];
        Real rhoi = rho_[index_i];
        Real pi = p_[index_i];
        Real c = 1.0;
    }
    if (mass_[index_i] < 0)
    {
        Real pr = p_[index_i];
        Real rho = rho_[index_i];
        Real Mass = mass_[index_i];
        Real y = 1.0;
    }
}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration2ndHalf<Contact<Wall>, RiemannSolverType>::
    EulerianIntegration2ndHalf(BaseContactRelation &wall_contact_relation, Real limiter_parameter)
    : BaseEulerianIntegrationWithWall(wall_contact_relation),
      riemann_solver_(this->fluid_, this->fluid_, limiter_parameter){};
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration2ndHalf<Contact<Wall>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidStateIn state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);
    Real mass_change_rate = 0.0;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        StdLargeVec<Real> &Vol_k = *(this->wall_Vol_[k]);
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * Vol_k[index_j];

            Vecd vel_in_wall = -state_i.vel_;
            Real p_in_wall = state_i.p_;
            Real rho_in_wall = state_i.rho_;

            FluidStateIn state_j(rho_in_wall, vel_in_wall, p_in_wall);
            FluidStateOut interface_state = this->riemann_solver_.InterfaceState(state_i, state_j, n_k[index_j]);
            mass_change_rate -= 2.0 * this->Vol_[index_i] * (interface_state.rho_ * interface_state.vel_).dot(e_ij) * dW_ijV_j;
        }
    }
    this->dmass_dt_[index_i] += mass_change_rate;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_FLUID_INTEGRATION_HPP
