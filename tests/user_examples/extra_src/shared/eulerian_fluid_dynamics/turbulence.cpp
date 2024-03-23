#pragma once
#include "turbulence.h"

namespace SPH
{  
    namespace fluid_dynamics
    {
        TurbulenceModel_K_Epsilon::
            TurbulenceModel_K_Epsilon(Real& C_mu) : C_mu_(C_mu){}


        //=================================================================================================//
        Real TurbulenceModel_K_Epsilon::TurbulentViscosity(Real& rho_star, Vecd& K_initial, Vecd& Epsilon_initial)
        {
            Real mu_t_ = rho_star * C_mu_ * (K_initial.squaredNorm()) / (Epsilon_initial.norm());
            return mu_t_;
        }
        //=================================================================================================//
        Matd TurbulenceModel_K_Epsilon::Production_TKE(Vecd& K_initial, FluidStarState& interface_state, Real& dW_ijV_j, Real& mu_t_initial)
        {
            Matd P;
            for (int i = 0; i < K_initial.size(); ++i) 
            {
                for (int j = 0; j < K_initial.size(); ++j) 
                {
                    P(i,j) = (2 * mu_t_initial * (interface_state.vel_[i] + interface_state.vel_[j]) * dW_ijV_j
                    - (2.0 / 3.0) * K_initial.norm() * (i == j ? 1 : 0))
                    * 2 * interface_state.vel_[i] * dW_ijV_j;
                }
            }
            return P;
        }

        //=================================================================================================//
       Vecd TurbulenceModel_K_Epsilon::TkeTransportEqn(Vecd& K_initial, Vecd& Epsilon_initial, Real& rho_star, FluidStarState& interface_state,
            Real& mu_t_initial, Real& mu_, Vecd& TKE_change_rate, Real& dW_ijV_j, Real &r_ij, Vecd &e_ij, Real& dt, Matd& P_TKE)
       {
           
            Real sigma_k = 1.0; // Turb Prndtl user should have a chance to set it!! 
            TKE_change_rate = TKE_change_rate -2 * (K_initial.norm() * rho_star * interface_state.vel_) * dW_ijV_j
                                + 2 * (mu_ + (mu_t_initial / sigma_k)) * ((K_initial.norm()) / (rho_star * r_ij)) * dW_ijV_j * e_ij
                                - rho_star * Epsilon_initial.norm() * e_ij
                                + P_TKE * e_ij;

            Vecd TKE = TKE_change_rate * dt;
            return TKE;
        }
        //=================================================================================================//
        Vecd TurbulenceModel_K_Epsilon::EpsilonTransportEqn(Vecd& K_updated, Vecd& Epsilon_initial, Real& rho_star, FluidStarState& interface_state,
            Real& mu_t_initial, Real& mu_, Vecd& Epsilon_change_rate, Real& dW_ijV_j, Real &r_ij, Vecd &e_ij, Real& dt, Matd& P_TKE)
            {
                Real sigma_eps = 1.3; // Ref Wilcox Turb Prndtl number for epsilon
                Real C1_eps = 1.44;
                Real C2_eps = 1.92;
               Epsilon_change_rate = -2 * (Epsilon_initial.norm() * rho_star * interface_state.vel_) * dW_ijV_j
                                    + 2 * (mu_ + (mu_t_initial / sigma_eps)) * ((Epsilon_initial.norm()) / (rho_star * r_ij)) * dW_ijV_j * e_ij
                                    - C2_eps * rho_star * ((Epsilon_initial.squaredNorm()) / (K_updated.norm())) * e_ij
                                    + C1_eps * ((Epsilon_initial.norm()) / (K_updated.norm())) * (P_TKE * e_ij);
                Vecd DissipationRate = Epsilon_change_rate * dt;
                return DissipationRate;
            }
            //=================================================================================================//
            Real TurbulenceModel_K_Epsilon::CorrectedTurbulentViscosity(Real& rho_star, Vecd& K_updated, Vecd& Epsilon_updated)
            {
                Real correctedMu_t_ = rho_star * C_mu_ * (K_updated.squaredNorm()) / (Epsilon_updated.norm());
                return correctedMu_t_;
            }
            //=================================================================================================//
            Matd TurbulenceModel_K_Epsilon::ReynoldStressTensor(Vecd& K_updated, FluidStarState& interface_state, Real& dW_ijV_j, Real& mu_t_corrected)
            {
                Matd R;
            for (int i = 0; i < K_updated.size(); ++i) 
            {
                for (int j = 0; j < K_updated.size(); ++j) 
                {
                    R(i,j) = (2 * mu_t_corrected * (interface_state.vel_[i] + interface_state.vel_[j]) * dW_ijV_j
                    - (2.0 / 3.0) * K_updated.norm() * (i == j ? 1 : 0))
                    * 2 * interface_state.vel_[i] * dW_ijV_j;
                }
            }
            return R;
            }

    }
}// namespace SPH