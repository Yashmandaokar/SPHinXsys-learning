#pragma once
#include "TurbulenceModel.h"

namespace SPH
{  
    namespace fluid_dynamics
    {
        MeanVelocity::MeanVelocity(BaseInnerRelation& inner_relation) :
            BaseIntegration(inner_relation), Vol_(particles_->Vol_)
        {}
        //=================================================================================================//
        void MeanVelocity::interaction(size_t index_i, Real dt)
        {
            Vecd velocity = Vecd::Zero();
            Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
                Real& W_ij = inner_neighborhood.W_ij_[n];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Vecd& e_ij = inner_neighborhood.e_ij_[n];
                velocity += (vel_[index_j]) * (W_ij * Vol_[index_j]);
            }
            meanvelocity_[index_i] = velocity;
        }
        //=================================================================================================//
        
        void MeanVelocity::update(size_t index_i, Real dt)
        {
            meanvelocity_[index_i] = meanvelocity_[index_i]; // Check with supervisor
        }
        //=================================================================================================//
        BaseTurbulence::BaseTurbulence(BaseInnerRelation &inner_relation) : 
        MeanVelocity(inner_relation),
        Cmu_(0.09), sigmak_(1.0), sigmaeps_(1.3), C1eps_(1.44), C2eps_(1.92), 
        Vol_(particles_->Vol_), K_(*particles_->getVariableByName<Real>("TKE")),
        Eps_(*particles_->getVariableByName<Real>("Dissipation")),
        mu_t_(*particles_->getVariableByName<Real>("TurbulentViscosity")),
        meanvelocity_(*particles_->getVariableByName<Vecd>("MeanVelocity"))
        {}
        
        //=================================================================================================//
        KEpsilonStd1stHalf::KEpsilonStd1stHalf(BaseInnerRelation &inner_relation) :
        BaseTurbulence(inner_relation)
        {}
        //=================================================================================================//
        
        void KEpsilonStd1stHalf::interaction(size_t index_i, Real dt)
        {
            mu_t_[index_i] = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i])/ Eps_[index_i]);
            
            Real K_changerate = 0.0, K_advection = 0.0;
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            Vecd K_laplacian = Vecd::Zero();  
            Matd K_prod, vel_gradient_mat;
            K_prod.setZero();
            vel_gradient_mat.setZero();
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Vecd &e_ij = inner_neighborhood.e_ij_[n];
                Real W_ij = inner_neighborhood.W_ij_[n];

                vel_gradient_mat = (vel_[index_i] - vel_[index_j]) * e_ij.transpose();
                K_advection = dW_ijV_j * (rho_[index_i] * K_[index_i] * vel_[index_i] - rho_[index_j] * K_[index_j] * vel_[index_j]).dot(e_ij);
                K_laplacian = 2 * dW_ijV_j * ((fluid_.ReferenceViscosity() + mu_t_[index_i] / sigmak_) * (K_[index_i] - K_[index_j]) / (r_ij)) * e_ij;
                K_prod = (mu_t_[index_i] * dW_ijV_j * (vel_gradient_mat + vel_gradient_mat.transpose())
                        -(2.0 / 3.0) * W_ij * Vol_[index_j] * (rho_[index_j] * K_[index_j]) * Matd::Identity()).array() * (dW_ijV_j * (vel_gradient_mat)).array();

                K_changerate += -K_advection + K_prod.sum() - Vol_[index_j] * W_ij * (rho_[index_j] * Eps_[index_j]) + K_laplacian.sum();   
            }
            dK_dt_[index_i] =  K_changerate;
        }
        //=================================================================================================//
        void KEpsilonStd1stHalf:: update(size_t index_i, Real dt)
        {
            K_[index_i] += dK_dt_[index_i] * dt;
        }
        //=================================================================================================//
        KEpsilonStd2ndHalf::KEpsilonStd2ndHalf(BaseInnerRelation &inner_relation) :
        BaseTurbulence(inner_relation)
        {}
        //=================================================================================================//
        
        void KEpsilonStd2ndHalf::interaction(size_t index_i, Real dt)
        {
            mu_t_[index_i] = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i])/ Eps_[index_i]);
            Real Eps_changerate = 0.0, Eps_advection = 0.0;
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            Vecd Eps_laplacian = Vecd::Zero();  
            Matd K_prod = Matd::Zero(), vel_gradient_mat = Matd::Zero(); // check
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Vecd &e_ij = inner_neighborhood.e_ij_[n];
                Real &W_ij = inner_neighborhood.W_ij_[n]; 

                vel_gradient_mat = (vel_[index_i] - vel_[index_j]) * e_ij.transpose();
                Eps_advection =  dW_ijV_j * (rho_[index_i] * Eps_[index_i] * vel_[index_i] - rho_[index_j] * Eps_[index_j] * vel_[index_j]).dot(e_ij);
                Eps_laplacian =  2 * ((fluid_.ReferenceViscosity() + mu_t_[index_i] / sigmaeps_) * ( Eps_[index_i] - Eps_[index_j]) / (r_ij) * dW_ijV_j) * e_ij; 
                K_prod = (mu_t_[index_i] * dW_ijV_j * (vel_gradient_mat + vel_gradient_mat.transpose()) 
                        - (2 / 3) * W_ij * Vol_[index_j] * (rho_[index_j] * K_[index_j]) * Matd::Identity()).array() * (dW_ijV_j * (vel_gradient_mat)).array();
                
                
                Eps_changerate += -Eps_advection + Eps_laplacian.sum() + C1eps_ * ((Eps_[index_i] / K_[index_i]) - (Eps_[index_j] / K_[index_j])) * K_prod.sum()
                                    - C2eps_ * Vol_[index_j] * W_ij * rho_[index_j] * (Eps_[index_j] * Eps_[index_j]) / K_[index_j];
            }
            dEps_dt_[index_i] = Eps_changerate;
        }
        //=================================================================================================//
        void KEpsilonStd2ndHalf:: update(size_t index_i, Real dt)
        {
            Eps_[index_i] += dEps_dt_[index_i] * dt;
        }
        //=================================================================================================//
        FluidStarState NoRiemannSolverInCompressobleRANS::getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij)
        {
            Real p_star = 0.5 * (state_i.p_ + state_j.p_);
            Vecd v_star = 0.5 * (state_i.vel_ + state_j.vel_);
            Real rho_star = 0.5 * (state_i.rho_ + state_j.rho_);

            FluidStarState interface_state(v_star, p_star);
            interface_state.vel_ = v_star;
            interface_state.p_ = p_star;

            return interface_state;
        }
        //=================================================================================================//
        StdWallFunction::StdWallFunction(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator)
        :  BaseTurbulence(inner_relation), ghost_creator_(ghost_creator), vonkar_(0.4187), E_(9.793)
        {}
        //=================================================================================================//
        void StdWallFunction::interaction(size_t index_i, Real dt)
        {
            for (size_t boundary_type = 0; boundary_type < ghost_creator_.each_boundary_type_with_all_ghosts_index_.size(); ++boundary_type)
            {
                if (!ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type].empty())
                {
                    for (size_t ghost_number = 0; ghost_number != ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type].size(); ++ghost_number)
                    {
                        size_t ghost_index = ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type][ghost_number];
                        size_t index_real = ghost_creator_.each_boundary_type_contact_real_index_[boundary_type][ghost_number];
                        Vecd e_ij = ghost_creator_. each_boundary_type_with_all_ghosts_eij_[boundary_type][ghost_number];
                        if (boundary_type == 3)
                        {
                           Real ystar = (rho_[index_real] * std::pow(Cmu_, 0.25) * std::pow(K_[index_real], 0.5) * pos_[index_real][1]) / (fluid_.ReferenceViscosity()); //corect yp
                           Real u_star, wall_shear_stress, mu_eff_wall;
                           Vecd distvec = (pos_[index_real] - pos_[ghost_index]).cwiseAbs();
                           Real normaldist = distvec.norm();
                           Vecd veltangential = (vel_[index_i] - e_ij.dot(vel_[index_i]) * (e_ij));
                           if (ystar >= 11.225)
                           {
                            u_star = (1 / vonkar_) * std::log(E_ * ystar);
                            wall_shear_stress = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(K_[index_real], 0.5) * rho_[index_real]) / (u_star);
                            mu_eff_wall = wall_shear_stress * normaldist / (veltangential.norm());
                           }
                           else if (ystar < 11.225)
                           {
                            u_star = ystar;
                            wall_shear_stress = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(K_[index_real], 0.5) * rho_[index_real]) / (u_star);
                            mu_eff_wall = wall_shear_stress * normaldist / (veltangential.norm());
                           }
                           mu_t_[index_real] = mu_eff_wall - fluid_.ReferenceViscosity();
                          
                           Real K_changerate_p = 0.0, K_prod_p = 0.0, Eps_p = 0.0, K_advection_p = 0.0 ;
                           Vecd K_laplacian_p = Vecd::Zero();
                           Neighborhood &inner_neighborhood = inner_configuration_[index_real];
                           for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                           {
                                size_t index_j = inner_neighborhood.j_[n];
                                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
                                Real r_ij = inner_neighborhood.r_ij_[n];
                                Vecd &e_ij = inner_neighborhood.e_ij_[n];
                                Real &W_ij = inner_neighborhood.W_ij_[n]; 
                                
                                K_advection_p = dW_ijV_j * (rho_[index_real] * K_[index_real] * vel_[index_real] - rho_[index_j] * K_[index_j] * vel_[index_j]).dot(e_ij);
                                K_laplacian_p = 2 * (((mu_eff_wall - fluid_.ReferenceViscosity()) / (sigmak_) + fluid_.ReferenceViscosity()) * ((K_[index_real] - K_[index_j]) / r_ij * dW_ijV_j)) * e_ij;
                                
                                K_prod_p = std::pow(wall_shear_stress,2) / (vonkar_ * rho_[index_real] * std::pow(Cmu_, 0.25) * std::pow(K_[index_real], 0.5) * normaldist);
                                Eps_p  =  (std::pow(Cmu_, 3/4) * std::pow(K_[index_real], 1.5)) / (vonkar_ * normaldist);

                                K_changerate_p += -K_advection_p + K_prod_p - Vol_[index_j] * W_ij * (rho_[index_j] * Eps_p) + K_laplacian_p.sum(); 
                            }
                            dK_dt_[index_real] =  K_changerate_p;
                            Esp_p_[index_real] = Eps_p;
                        }   
                    }

                }

            }
        } 

        //=================================================================================================//
        void StdWallFunction::update(size_t index_i, Real dt)
        {
            K_[index_i] += dK_dt_[index_i]* dt;
            Eps_[index_i] +=  Esp_p_[index_i];
        }  
    }


}// namespace SPH
