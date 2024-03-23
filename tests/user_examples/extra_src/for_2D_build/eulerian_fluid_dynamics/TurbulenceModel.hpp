#include "eulerian_fluid_dynamics.h"
#include "TurbulenceModel.h"


namespace SPH
{
    namespace fluid_dynamics
    {
        template <class RiemannSolverType>
        EulerianIntegration1stHalfRANS<RiemannSolverType>::
            EulerianIntegration1stHalfRANS(BaseInnerRelation &inner_relation, Real limiter_parameter)
            : BaseTurbulence(inner_relation), limiter_input_(limiter_parameter),
            riemann_solver_(this->fluid_, this->fluid_, limiter_input_),
            acc_prior_(particles_->acc_prior_)
        {
            particles_->registerVariable(mom_, "Momentum");
            particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
        }
        
        //=================================================================================================//
        template <class RiemannSolverType>
        void EulerianIntegration1stHalfRANS<RiemannSolverType>::interaction(size_t index_i, Real dt)
        {
            FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
            Vecd momentum_change_rate = Vecd::Zero();
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
           
            mu_t_[index_i] = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i])/ Eps_[index_i]);
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
                Vecd &e_ij = inner_neighborhood.e_ij_[n];
                Real &r_ij = inner_neighborhood.r_ij_[n];

                FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
                FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
                Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);
                Matd rey_stresstensor = Matd::Zero();
                Matd meanvel_advection = Matd::Zero(), viscous_dissipation = Matd::Zero();
                Matd meanvelocitylaplacian = interface_state.vel_ * e_ij.transpose();

                meanvel_advection = dW_ijV_j * (rho_star * interface_state.vel_) * interface_state.vel_.transpose();
                viscous_dissipation = 2 * mu_t_[index_i] * dW_ijV_j * (meanvelocitylaplacian) / r_ij;
                rey_stresstensor = 2 * mu_t_[index_i] * dW_ijV_j * (meanvelocitylaplacian + meanvelocitylaplacian.transpose()) / r_ij
                                    - (2 / 3) * dW_ijV_j * (rho_star * K_[index_i] * Matd::Identity()); // TKE riemansolution



                momentum_change_rate += (-meanvel_advection + viscous_dissipation - interface_state.p_ * Matd::Identity() + rey_stresstensor) * e_ij; //reystress sign
            }
            dmom_dt_[index_i] = momentum_change_rate;
        }
        //=================================================================================================//
        template <class RiemannSolverType>
        void EulerianIntegration1stHalfRANS<RiemannSolverType>::update(size_t index_i, Real dt)
        {
            mom_[index_i] += (dmom_dt_[index_i] + rho_[index_i] * acc_prior_[index_i]) * dt;
            vel_[index_i] = mom_[index_i] / rho_[index_i];
        }

        //=================================================================================================//
        template <class RiemannSolverType>
        EulerianIntegration2ndHalfRANS<RiemannSolverType>::EulerianIntegration2ndHalfRANS(BaseInnerRelation &inner_relation, Real limiter_parameter)
        : BaseTurbulence(inner_relation), limiter_input_(limiter_parameter),
          riemann_solver_(this->fluid_, this->fluid_, limiter_input_) {}
        
        //=================================================================================================//
        template <class RiemannSolverType>
        void EulerianIntegration2ndHalfRANS<RiemannSolverType>::interaction(size_t index_i, Real dt)
        {
            FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
            Real density_change_rate = 0.0;
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Vecd &e_ij = inner_neighborhood.e_ij_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

                FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
                FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

                Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);
                density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
            }
            drho_dt_[index_i] = density_change_rate;
        }
        //=================================================================================================//
        template <class RiemannSolverType>
        void EulerianIntegration2ndHalfRANS<RiemannSolverType>::update(size_t index_i, Real dt)
        {
            rho_[index_i] += drho_dt_[index_i] * dt;
            p_[index_i] = fluid_.getPressure(rho_[index_i]);
        }
        //=================================================================================================//
    }// namespace SPH

}// namespace SPH
