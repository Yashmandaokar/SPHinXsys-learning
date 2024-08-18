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
            Kprof_(*particles_->getVariableByName<Real>("TKEProfile")),
            Epsprof_(*particles_->getVariableByName<Real>("DissipationProfile")),
            mu_tprof_(*particles_->getVariableByName<Real>("TurblunetViscosityProfile"))
            , meanvel_advection_(*this->particles_->template registerSharedVariable<Vecd>("MomentumAdvection")),
            viscous_dissipation_(*this->particles_->template registerSharedVariable<Vecd>("MomentumViscousDissipation")),
            pressuregrad_(*this->particles_->template registerSharedVariable<Vecd>("MomentumPressureGradient")),
            tkegrad_(*this->particles_->template registerSharedVariable<Vecd>("MomentumTKEGradient")),
            vel_gradient_mat_(*particles_->getVariableByName<Matd>("VelocityGradient")),
            walladjacentcellflag_(*particles_->getVariableByName<Real>("FlagForWallAdjacentCells")),
            Tau_wall_(*particles_->getVariableByName<Real>("WallShearStress")),
            wallfacearea_(*particles_->getVariableByName<Real>("WallFaceArea")),
            shear_force_(*this->particles_->template registerSharedVariable<Vecd>("ShearForce"))
            {}
        //=================================================================================================//

        template <class RiemannSolverType>
        void EulerianIntegration1stHalfRANS<Inner<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
        {
            //mu_t_[index_i] = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));
            FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
            Vecd momentum_change_rate = Vecd::Zero();
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            //Matd rey_stresstensor = Matd::Zero();
            meanvel_advection_[index_i] = Vecd::Zero(), pressuregrad_[index_i] = Vecd::Zero(), tkegrad_[index_i] = Vecd::Zero();
            viscous_dissipation_[index_i] = Vecd::Zero(), shear_force_[index_i] = Vecd::Zero();
            Vecd flow_vector(1.0, 0.0);
           
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ij = inner_neighborhood.dW_ij_[n];
                Vecd &e_ij = inner_neighborhood.e_ij_[n];
                Real &r_ij = inner_neighborhood.r_ij_[n];
                Real mu_t_avg = (2.0 * mu_tprof_[index_i] * mu_tprof_[index_j]) / (mu_tprof_[index_i] + mu_tprof_[index_j] + TinyReal);
                
                FluidStateIn state_j(rho_[index_j], vel_[index_j], p_[index_j]);
                FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, e_ij);

                
                meanvel_advection_[index_i] += Vol_[index_i] * 2.0 * (-dW_ij * Vol_[index_j] * (interface_state.rho_ * (interface_state.vel_) * (interface_state.vel_).transpose())) * e_ij;
                pressuregrad_[index_i] += Vol_[index_i] * 2.0 * (-dW_ij * Vol_[index_j] * (interface_state.p_) * Matd::Identity()) * e_ij;
                tkegrad_[index_i] += Vol_[index_i] * 2.0 * (-dW_ij * Vol_[index_j] * interface_state.rho_ * (2.0 / 3.0) * (Kprof_[index_i]  - Kprof_[index_j]) * Matd::Identity()) * e_ij;
                viscous_dissipation_[index_i] += Vol_[index_i] * (2.0 * (fluid_.ReferenceViscosity() + mu_t_avg) * dW_ij * Vol_[index_j]) * (vel_[index_i] - vel_[index_j]) / r_ij;
                //shear_force_[index_i] += -Tau_wall_[index_j] * wallfacearea_[index_i] * flow_vector;
                
                momentum_change_rate = meanvel_advection_[index_i] + pressuregrad_[index_i] + tkegrad_[index_i] + viscous_dissipation_[index_i];

                if (index_i == 2732)
                { 
                    Vecd visdiss = viscous_dissipation_[index_i];
                    Vecd prgrad = pressuregrad_[index_i];
                    Vecd adv = meanvel_advection_[index_i];
                    Vecd tkegrad = tkegrad_[index_i];
                    Real tau = Tau_wall_[index_j];
                    Real wallarea = wallfacearea_[index_i];
                    Real voli = Vol_[index_i];
                    Real volj = Vol_[index_j];
                    Real c = 1.0;
                }
                
            }
            dmom_dt_[index_i] = momentum_change_rate;
        }
        //=================================================================================================//
        template <class RiemannSolverType>
        void EulerianIntegration1stHalfRANS<Inner<>, RiemannSolverType>::update(size_t index_i, Real dt)
        {
            mom_[index_i] += (dmom_dt_[index_i] +  force_prior_[index_i]) * dt;
            vel_[index_i] = mom_[index_i] / mass_[index_i];

            if (index_i == 2732)
            {
                Vecd veli = vel_[index_i];
                Vecd momi = mom_[index_i];
                Real c = 1.0;
            }
          if (std::isnan(mom_[index_i].norm()))
            {
                Vecd momi = mom_[index_i];
                Vecd foce = force_prior_[index_i];
                Vecd veli = vel_[index_i];
                Real X = 1.0;
            }
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
