#ifndef RANSINTEGRATION_HPP
#define RANSINTEGRATION_HPP
#include "ransintegration.h"

namespace SPH
{
    namespace fluid_dynamics
    {

        template <class DataDelegationType>
        template <class BaseRelationType>
        IntegrationRANS<DataDelegationType>::IntegrationRANS(BaseRelationType& base_relation)
            : BaseIntegration<DataDelegationType>(base_relation),
            mom_(*this->particles_->template registerSharedVariable<Vecd>("Momentum")),
            dmom_dt_(*this->particles_->template registerSharedVariable<Vecd>("MomentumChangeRate")),
            dmass_dt_(*this->particles_->template registerSharedVariable<Real>("MassChangeRate")),
            Vol_(*this->particles_->template getVariableByName<Real>("VolumetricMeasure")) {}
        //=================================================================================================//
        template <class RiemannSolverType>
        Integration1stHalfRANS<Inner<>, RiemannSolverType>::
            Integration1stHalfRANS(BaseInnerRelation &inner_relation, Real limiter_parameter)
            : IntegrationRANS<DataDelegateInner>(inner_relation),
            riemann_solver_(this->fluid_, this->fluid_, limiter_parameter) {}
        //=================================================================================================//

        template <class RiemannSolverType>
        void Integration1stHalfRANS<Inner<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
        {
            FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
            Vecd momentum_change_rate = Vecd::Zero();
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            
            Matd meanvel_advection = Matd::Zero(), pressuregrad = Matd::Zero();
            Vecd viscous_dissipation = Vecd::Zero();
           
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ij = inner_neighborhood.dW_ij_[n];
                Vecd &e_ij = inner_neighborhood.e_ij_[n];
                Real &r_ij = inner_neighborhood.r_ij_[n];
                
                FluidStateIn state_j(rho_[index_j], vel_[index_j], p_[index_j]);
                FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, e_ij);
                
                meanvel_advection = -dW_ij * Vol_[index_j] * (interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose());
                viscous_dissipation = 2.0 * (fluid_.ReferenceViscosity() * dW_ij * Vol_[index_j]) * (vel_[index_i] - vel_[index_j]) / r_ij;
                pressuregrad = -dW_ij * Vol_[index_j] * (interface_state.p_) * Matd::Identity();
               
                Vecd adv = Vol_[index_i] * 2.0 * meanvel_advection * e_ij;
                Vecd pr = Vol_[index_i] * 2.0 * pressuregrad * e_ij;
                //Vecd vis = Vol_[index_i] * viscous_dissipation;

                momentum_change_rate += Vol_[index_i] * (2.0 * (meanvel_advection + pressuregrad) * e_ij + viscous_dissipation);   

                Vecd momrate = Vol_[index_i] * (2.0 * (meanvel_advection + pressuregrad) * e_ij + viscous_dissipation);
            }
            dmom_dt_[index_i] = momentum_change_rate;
        }
        //=================================================================================================//
        template <class RiemannSolverType>
        void Integration1stHalfRANS<Inner<>, RiemannSolverType>::update(size_t index_i, Real dt)
        {
            mom_[index_i] += (dmom_dt_[index_i] +  force_prior_[index_i]) * dt;
            vel_[index_i] = mom_[index_i] / mass_[index_i];
        }
        //=================================================================================================//
    }// namespace SPH

}// namespace SPH

#endif // RANSINTEGRATION_HPP
