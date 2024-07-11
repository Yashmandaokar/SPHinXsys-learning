#ifndef RANS_FLUID_INTEGRATION_HPP
#define RANS_FLUID_INTEGRATION_HPP
#include "rans_fluid_integration.h"

namespace SPH
{
    namespace fluid_dynamics
    {

        template <class DataDelegationType>
        template <class BaseRelationType>
        EulerianIntegrationRANS<DataDelegationType>::EulerianIntegrationRANS(BaseRelationType& base_relation)
            : BaseIntegration<DataDelegationType>(base_relation),
            mom_(*this->particles_->template registerSharedVariable<Vecd>("Momentum")),
            dmom_dt_(*this->particles_->template registerSharedVariable<Vecd>("MomentumChangeRate")),
            dmass_dt_(*this->particles_->template registerSharedVariable<Real>("MassChangeRate")),
            Vol_(*this->particles_->template getVariableByName<Real>("VolumetricMeasure")), Cmu_(0.09) {}
        //=================================================================================================//
        //=================================================================================================//
        template <class RiemannSolverType>
        EulerianIntegration1stHalfRANS<Inner<>, RiemannSolverType>::
            EulerianIntegration1stHalfRANS(BaseInnerRelation& inner_relation, Real limiter_parameter)
            : EulerianIntegrationRANS<DataDelegateInner>(inner_relation),
            riemann_solver_(this->fluid_, this->fluid_, limiter_parameter), 
            K_(*particles_->getVariableByName<Real>("TKE")),
            Eps_(*particles_->getVariableByName<Real>("Dissipation")),
            mu_t_(*particles_->getVariableByName<Real>("TurbulentViscosity")) {}
        //=================================================================================================//

        template <class RiemannSolverType>
        void EulerianIntegration1stHalfRANS<Inner<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
        {
            FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
            Vecd momentum_change_rate = Vecd::Zero();
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            //Matd rey_stresstensor = Matd::Zero();
            Matd meanvel_advection = Matd::Zero(), pressuregrad = Matd::Zero(), tkegrad = Matd::Zero();
            Vecd viscous_dissipation = Vecd::Zero();
           
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ij = inner_neighborhood.dW_ij_[n];
                Vecd &e_ij = inner_neighborhood.e_ij_[n];
                Real &r_ij = inner_neighborhood.r_ij_[n];
                Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);
                FluidStateIn state_j(rho_[index_j], vel_[index_j], p_[index_j]);
                FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, e_ij);
                
                meanvel_advection = -dW_ij * Vol_[index_j] * (interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose());
                viscous_dissipation = (2.0 * (fluid_.ReferenceViscosity() + mu_t_avg) * dW_ij * Vol_[index_j]) * (vel_[index_i] - vel_[index_j]) / r_ij;
                pressuregrad = -dW_ij * Vol_[index_j] * (interface_state.p_) * Matd::Identity();
                tkegrad = -dW_ij * Vol_[index_j] * interface_state.rho_ * (2.0 / 3.0) * (K_[index_i] - K_[index_j]) * Matd::Identity();
                //rey_stresstensor = 2.0 * mu_t_[index_i] * dW_ij * Vol_[index_j] * (meanvelocitylaplacian + meanvelocitylaplacian.transpose()) / r_ij
                                //- (2.0 / 3.0) * dW_ij * Vol_[index_j] * (interface_state.rho_ * (K_[index_i]) * Matd::Identity());
                //rey_stresstensor = mu_t_[index_i] * dW_ij * dW_ij * Vol_[index_j] * Vol_[index_j] * (meanvelocitylaplacian + meanvelocitylaplacian.transpose())
                                 //- (2.0 / 3.0) * dW_ij * Vol_[index_j] * (interface_state.rho_ * (K_[index_i]) * Matd::Identity());
                momentum_change_rate += Vol_[index_i] * (2.0 * (meanvel_advection + pressuregrad + tkegrad) * e_ij + viscous_dissipation);   
            }
            dmom_dt_[index_i] = momentum_change_rate;
        }
        //=================================================================================================//
        template <class RiemannSolverType>
        void EulerianIntegration1stHalfRANS<Inner<>, RiemannSolverType>::update(size_t index_i, Real dt)
        {
            mom_[index_i] += (dmom_dt_[index_i] +  force_prior_[index_i]) * dt;
            vel_[index_i] = mom_[index_i] / mass_[index_i];
        }
        //=================================================================================================//
        //=================================================================================================//
        template <class RiemannSolverType>
        EulerianIntegration1stHalfRANS<Contact<Wall>, RiemannSolverType>::
            EulerianIntegration1stHalfRANS(BaseContactRelation& wall_contact_relation, Real limiter_parameter)
            : BaseEulerianIntegrationWithWallRANS(wall_contact_relation), riemann_solver_(fluid_, fluid_, limiter_parameter),
            K_(*particles_->getVariableByName<Real>("TKE")),
            Eps_(*particles_->getVariableByName<Real>("Dissipation")),
            mu_t_(*particles_->getVariableByName<Real>("TurbulentViscosity")) {}
        //=================================================================================================//
        template <class RiemannSolverType>
        void EulerianIntegration1stHalfRANS<Contact<Wall>, RiemannSolverType>::interaction(size_t index_i, Real dt)
        {
            FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
            mu_t_[index_i] = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i]) / Eps_[index_i]);
            Vecd momentum_change_rate = Vecd::Zero();
            for (size_t k = 0; k < contact_configuration_.size(); ++k)
            {
                StdLargeVec<Vecd>& n_k = *(wall_n_[k]);
                Neighborhood& wall_neighborhood = (*contact_configuration_[k])[index_i];
                for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
                {
                    size_t index_j = wall_neighborhood.j_[n];
                    Vecd& e_ij = wall_neighborhood.e_ij_[n];
                    Real dW_ij = wall_neighborhood.dW_ij_[n];
                    Real& r_ij = wall_neighborhood.r_ij_[n];

                    Vecd vel_in_wall = -state_i.vel_;
                    Real p_in_wall = state_i.p_;
                    Real rho_in_wall = state_i.rho_;
                    FluidStateIn state_j(rho_in_wall, vel_in_wall, p_in_wall);
                    FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, n_k[index_j]);
                    Matd rey_stresstensor = Matd::Zero();
                    Matd meanvel_advection = Matd::Zero(), viscous_dissipation = Matd::Zero();
                    Matd meanvelocitylaplacian = interface_state.vel_ * e_ij.transpose();
                    meanvel_advection = dW_ij * Vol_[index_j] * (interface_state.rho_ * interface_state.vel_) * interface_state.vel_.transpose();
                    viscous_dissipation = 2 * mu_t_[index_i] * dW_ij * Vol_[index_j] * (meanvelocitylaplacian) / r_ij;
                    rey_stresstensor = 2 * mu_t_[index_i] * dW_ij * Vol_[index_j] * (meanvelocitylaplacian + meanvelocitylaplacian.transpose()) / r_ij
                        - (2.0 / 3.0) * dW_ij * Vol_[index_j] * (interface_state.rho_ * K_[index_i] * Matd::Identity()); // TKE riemansolution

                    momentum_change_rate += Vol_[index_i] * (-meanvel_advection + viscous_dissipation - interface_state.p_ * Matd::Identity() + rey_stresstensor) * e_ij; //reystress sign
                }
            }
            dmom_dt_[index_i] += momentum_change_rate;
        }
        //=================================================================================================//
    }// namespace SPH

}// namespace SPH

#endif // RANS_FLUID_INTEGRATION_HPP
