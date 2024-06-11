#ifndef TURBULENCEMODEL_CPP
#define TURBULENCEMODEL_CPP
#include "TurbulenceModel.h"
#include "sphinxsys.h"
namespace SPH
{  
    namespace fluid_dynamics
    {
        /*template <class DataDelegationType>
        template <class BaseRelationType>
        BaseTurbulence<DataDelegationType>::BaseTurbulence(BaseRelationType& base_relation)
            : BaseIntegration<DataDelegationType>(base_relation),
            mom_(*this->particles_->template registerSharedVariable<Vecd>("Momentum")),
            dmom_dt_(*this->particles_->template registerSharedVariable<Vecd>("MomentumChangeRate")),
            dmass_dt_(*this->particles_->template registerSharedVariable<Real>("MassChangeRate")),
            Vol_(this->particles_->Vol_), Cmu_(0.09), sigmak_(1.0), sigmaeps_(1.3),
            C1eps_(1.44), C2eps_(1.92), K_(*particles_->getVariableByName<Real>("TKE")),
            Eps_(*particles_->getVariableByName<Real>("Dissipation")),
            mu_t_(*particles_->getVariableByName<Real>("TurbulentViscosity")) {}*/
        //=================================================================================================//
        BaseTurbulence::BaseTurbulence(BaseInnerRelation& inner_relation, GhostCreationFromMesh& ghost_creator)
            : BaseIntegration<FluidDataInner>(inner_relation), ghost_creator_(ghost_creator),
            mom_(*this->particles_->template registerSharedVariable<Vecd>("Momentum")),
            dmom_dt_(*this->particles_->template registerSharedVariable<Vecd>("MomentumChangeRate")),
            dmass_dt_(*this->particles_->template registerSharedVariable<Real>("MassChangeRate")),
            K_prod_p_(*this->particles_->template registerSharedVariable<Real>("TKEProductionInWallAdjCell")),
            K_prodsum_(*this->particles_->template registerSharedVariable<Real>("TKEProductionSum")),
            Eps_p_(*this->particles_->template registerSharedVariable<Real>("DissipationRateInWallAdjCell")),
            Eps_sum_(*this->particles_->template registerSharedVariable<Real>("DissipationSum")),
            Cmu_(0.09), sigmak_(1.0), sigmaeps_(1.3), C1eps_(1.44), C2eps_(1.92), 
            K_(*particles_->getVariableByName<Real>("TKE")),
            Eps_(*particles_->getVariableByName<Real>("Dissipation")),
            mu_t_(*particles_->getVariableByName<Real>("TurbulentViscosity")) {}
        //=================================================================================================//
        WallAdjacentCells::WallAdjacentCells(BaseInnerRelation& inner_relation, GhostCreationFromMesh& ghost_creator)
            : BaseTurbulence(inner_relation, ghost_creator) 
        {
            walladjacentcells();
        }
        //=================================================================================================//
        void WallAdjacentCells::walladjacentcells()
        {
            walladjacentindex_.resize(particles_-> particles_bound_);
            walladjacentcellflag_.resize(particles_->particles_bound_);
            wallghostindex_.resize(particles_->particles_bound_);
            walleij_.resize(particles_->particles_bound_);
            Real boundary_type = 3;
                
            if (!ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type].empty())
            {
                for (size_t ghost_number = 0; ghost_number != ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type].size(); ++ghost_number)
                {

                    size_t ghost_index = ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type][ghost_number];
                    Real index_real = ghost_creator_.each_boundary_type_contact_real_index_[boundary_type][ghost_number];
                    walleij_[index_real] = ghost_creator_.each_boundary_type_with_all_ghosts_eij_[boundary_type][ghost_number];
                    walladjacentindex_[index_real] = index_real;
                    walladjacentcellflag_[index_real] = (1);
                    wallghostindex_[index_real] = ghost_index;
                }
            }
        }
        //=================================================================================================//

        KEpsilonStd1stHalf::KEpsilonStd1stHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh& ghost_creator) 
            : WallAdjacentCells(inner_relation, ghost_creator),
            dK_dt_(*this->particles_->template registerSharedVariable<Real>("TKEChangeRate"))
        {}
        //=================================================================================================//
        void KEpsilonStd1stHalf::interaction(size_t index_i, Real dt)
        {
            if (walladjacentcellflag_[index_i] != 1)
            {
                mu_t_[index_i] = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));
            }
           
            if (mu_t_[index_i] < 0.0)
            {
                Real mu_ti = mu_t_[index_i];
                Real Ki = K_[index_i];
                Real Epsi = Eps_[index_i];
                Real J = 1;
            }
            Real K_advection = 0.0, K_laplacian = 0.0, K_advection1 = 0.0;
            Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            Matd K_prod, vel_gradient_mat, K_prod_iso, K_prod_total;
            K_prod.setZero();
            K_prod_iso.setZero();
            K_prod_total.setZero();
            vel_gradient_mat.setZero();
            Matd vel_matrix = Matd::Zero();
            Matd vel_gradient_mat_transpose = Matd::Zero();
            Matd strain_tensor = Matd::Zero();
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ij = inner_neighborhood.dW_ij_[n];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Vecd& e_ij = inner_neighborhood.e_ij_[n];
                //Real W_ij = inner_neighborhood.W_ij_[n];

                K_advection += (- 1.0 * dW_ij * Vol_[index_j] * rho_[index_i] * (K_[index_i] - K_[index_j]) * vel_[index_i]).dot(e_ij);
                K_advection1 = (-1.0 * dW_ij * Vol_[index_j] * rho_[index_i] * (K_[index_i] - K_[index_j]) * vel_[index_i]).dot(e_ij);
                K_laplacian += 2 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_[index_i] / sigmak_) * (K_[index_i] - K_[index_j]) / (r_ij));
                vel_matrix = (vel_[index_i] - vel_[index_j]) * e_ij.transpose();
                vel_gradient_mat = dW_ij * Vol_[index_j] * vel_matrix;
                vel_gradient_mat_transpose = dW_ij * Vol_[index_j] * vel_matrix.transpose();
                strain_tensor = 0.5 * (vel_gradient_mat + vel_gradient_mat_transpose);
                K_prod = mu_t_[index_i] * (strain_tensor.array() * strain_tensor.array());
                K_prod_iso = (2.0 / 3.0 * rho_[index_i] * K_[index_i] * Matd::Identity()).array() * vel_gradient_mat.array();
                K_prod_total = K_prod - K_prod_iso;
                //K_prod1 = (mu_t_[index_i] * dW_ij * Vol_[index_j] * (vel_gradient_mat + vel_gradient_mat.transpose())
                    //- (2.0 / 3.0) * (rho_[index_i] * K_[index_i]) * Matd::Identity()).array() * (dW_ij * Vol_[index_j] * (vel_gradient_mat)).array();
                K_prodsum_[index_i] += K_prod.sum();
                Eps_sum_[index_i] -= (rho_[index_i] * Eps_[index_i]);               
            }
            
            //Check if the cell is adjacent to wall
                if (walladjacentcellflag_[index_i] == 1)
                {
                    Real Kprodpcell = K_prod_p_[index_i];
                    Real Eps_pcell = Eps_p_[index_i];
                    dK_dt_[index_i] = K_advection + K_prod_p_[index_i] - rho_[index_i] * Eps_p_[index_i] + K_laplacian;
                    Real x = 1;
                }
                else
                {
                    Real Kprodpcell = K_prodsum_[index_i];
                    Real Eps_pcell = Eps_sum_[index_i];
                    dK_dt_[index_i] = K_advection + K_prodsum_[index_i] + Eps_sum_[index_i] + K_laplacian;
                    Real x = 2;
                }
        }
        //=================================================================================================//
        void KEpsilonStd1stHalf::update(size_t index_i, Real dt)
        {
            K_[index_i] += (dK_dt_[index_i] / rho_[index_i]) * dt;
            if (K_[index_i] < 0.0)
            {
                Real dki = dK_dt_[index_i];
                Real Ki = K_[index_i];
                Real J = 1;
            }
        }
        //=================================================================================================//
        KEpsilonStd2ndHalf::KEpsilonStd2ndHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh& ghost_creator) 
            : WallAdjacentCells(inner_relation, ghost_creator),
            dEps_dt_(*this->particles_->template registerSharedVariable<Real>("DissipationChangeRate"))
        {}
        //=================================================================================================//
        
        void KEpsilonStd2ndHalf::interaction(size_t index_i, Real dt)
        {
            if (walladjacentcellflag_[index_i] != 1)
            {
                mu_t_[index_i] = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));
            }
            if (mu_t_[index_i] < 0.0)
            {
                Real mu_ti = mu_t_[index_i];
                Real Ki = K_[index_i];
                Real Epsi = Eps_[index_i];
                Real J = 1;
            }
            Real Eps_changerate = 0.0, Eps_advection = 0.0, Eps_laplacian = 0.0;
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            Matd K_prod = Matd::Zero(), vel_gradient_mat = Matd::Zero(), K_prod_iso, K_prod_total;; // check
            K_prod_iso.setZero();
            K_prod_total.setZero();
            Matd vel_matrix = Matd::Zero();
            Matd vel_gradient_mat_transpose = Matd::Zero();
            Matd strain_tensor = Matd::Zero();
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ij = inner_neighborhood.dW_ij_[n];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Vecd &e_ij = inner_neighborhood.e_ij_[n];
                //Real &W_ij = inner_neighborhood.W_ij_[n]; 

                Eps_advection = dW_ij * Vol_[index_j] * rho_[index_i] * (vel_[index_i] * (Eps_[index_i] - Eps_[index_j])).dot(e_ij);
                Eps_laplacian =  2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_[index_i] / sigmaeps_) * ( Eps_[index_i] - Eps_[index_j]) / (r_ij));

                vel_matrix = (vel_[index_i] - vel_[index_j]) * e_ij.transpose();
                vel_gradient_mat = dW_ij * Vol_[index_j] * vel_matrix;
                vel_gradient_mat_transpose = dW_ij * Vol_[index_j] * vel_matrix.transpose();
                strain_tensor = 0.5 * (vel_gradient_mat + vel_gradient_mat_transpose);
                K_prod = mu_t_[index_i] * 2.0 * (strain_tensor.array() * strain_tensor.array());
                K_prod_iso = (2.0 / 3.0 * rho_[index_i] * K_[index_i] * Matd::Identity()).array() * vel_gradient_mat.array();
                K_prod_total = K_prod - K_prod_iso;
                //K_prod = (mu_t_[index_i] * dW_ij * Vol_[index_j] * (vel_gradient_mat + vel_gradient_mat.transpose())
                       // - (2.0 / 3.0) * (rho_[index_i] * K_[index_i]) * Matd::Identity()).array() * (dW_ij * Vol_[index_j] * (vel_gradient_mat)).array();
                
                Eps_changerate += -Eps_advection + C1eps_ * (Eps_[index_i] / (K_[index_i])) * K_prod.sum()
                        - C2eps_ * rho_[index_i] * (Eps_[index_i] * Eps_[index_i]) / (K_[index_i]) + Eps_laplacian;
                if (index_i == 103)
                {
                    Real mutj = mu_t_[index_j];
                    Real muti = mu_t_[index_i];
                    Real Ki = K_[index_i];
                    Real Kj = K_[index_j];
                    Real Epsi = Eps_[index_i];
                    Real Epsj = Eps_[index_j];
                    Vecd vel_i = vel_[index_i];
                    Vecd vel_j = vel_[index_j];
                    Real j = 1.0;
                }
            }
            dEps_dt_[index_i] = Eps_changerate;
        }
        //=================================================================================================//
        void KEpsilonStd2ndHalf::update(size_t index_i, Real dt)
        {
            if (walladjacentcellflag_[index_i] == 1)
            {
                Eps_[index_i] += Eps_p_[index_i];
            }
            else
            {
                Eps_[index_i] += (dEps_dt_[index_i] / rho_[index_i]) * dt;
            }
            if (Eps_[index_i] < 0.0)
            {
               
                Real Epsi = Eps_[index_i];
                Real J = 1;
            }
        }
        //=================================================================================================//
        StdWallFunctionFVM::StdWallFunctionFVM(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator)
        : WallAdjacentCells(inner_relation, ghost_creator), vonkar_(0.4187), E_(9.793)
            //dK_dt_(*particles_->getVariableByName<Real>("TKEChangeRate")),
            //K_prod_p_(*particles_->getVariableByName<Real>("TKEProductionInWallAdjCell")),
            //Esp_p_(*particles_->getVariableByName<Real>("DissipationRateInWallAdjCell")),
        {}
        //=================================================================================================//
        void StdWallFunctionFVM::interaction(size_t index_i, Real dt)
        {
            Real K_changerate_p = 0.0, K_prod_p = 0.0, Eps_p = 0.0, K_advection_p = 0.0;
            Vecd K_laplacian_p = Vecd::Zero();
            
            if (walladjacentcellflag_[index_i] == 1)
            {
                Real ghost_index = wallghostindex_[index_i];
                Real yp = (pos_[index_i] - pos_[ghost_index]).dot(walleij_[index_i]);
                Vecd posghost = pos_[ghost_index];
                Vecd posrealindex = pos_[index_i];
                Real ystar = (rho_[index_i] * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * yp) / (fluid_.ReferenceViscosity());
                Real u_star, wall_shear_stress, mu_eff_wall, friction_velocity, yplus, wall_shear_stress_1, friction_velocity_1, mu_eff_wall_1;
                Vecd veltangential = (vel_[index_i] - walleij_[index_i].dot(vel_[index_i]) * (walleij_[index_i]));
                if (ystar >= 11.225)
                {
                    if (index_i == 2622)
                    {
                        Real J = 1.0;
                    }
                    u_star = (1.0 / vonkar_) * std::log(E_ * ystar);
                    wall_shear_stress_1 = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * rho_[index_i]) / (u_star);
                    mu_t_[index_i] = fluid_.ReferenceViscosity() * ((ystar * vonkar_) / (std::log(E_ * ystar)) - 1.0);
                    if (mu_t_[index_i] < 0.0)
                    {
                        Real mu_ti = mu_t_[index_i];
                        Real Ki = K_[index_i];
                        Real Epsi = Eps_[index_i];
                        Real J = 1;
                    }
                    wall_shear_stress = (mu_t_[index_i] + fluid_.ReferenceViscosity()) * (veltangential.norm() / yp);
                    mu_eff_wall = wall_shear_stress * yp / (veltangential.norm());
                    mu_eff_wall_1 = wall_shear_stress_1 * yp / (veltangential.norm());
                    friction_velocity = std::sqrt(wall_shear_stress / rho_[index_i]);
                    friction_velocity_1 = std::sqrt(wall_shear_stress_1 / rho_[index_i]);
                                   
                    K_prod_p_[index_i] = std::pow(wall_shear_stress, 2.0) / (vonkar_ * rho_[index_i] * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * yp);
                    Eps_p_[index_i] = (std::pow(Cmu_, 3.0 / 4.0) * std::pow(K_[index_i], 1.5)) / (vonkar_ * yp);

                }
                else if (ystar < 11.225)
                {
                    u_star = ystar;
                    wall_shear_stress_1 = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * rho_[index_i]) / (u_star);
                    mu_t_[index_i] = 0.0;
                    wall_shear_stress = (mu_t_[index_i] + fluid_.ReferenceViscosity()) * (veltangential.norm() / yp);
                    mu_eff_wall = wall_shear_stress * yp / (veltangential.norm()); 
                    mu_eff_wall_1 = wall_shear_stress_1 * yp / (veltangential.norm());
                    friction_velocity = std::sqrt(wall_shear_stress / rho_[index_i]);
                    friction_velocity_1 = std::sqrt(wall_shear_stress_1 / rho_[index_i]);
                                  
                    K_prod_p_[index_i] = 0.0;
                    Eps_p_[index_i] = (K_[index_i] * 2.0 * fluid_.ReferenceViscosity()) / (rho_[index_i] * yp * yp);
                }
            } 
            mu_t_[index_i] = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));
        } 
        //=================================================================================================//
         
        //=================================================================================================//
        StdWallFunction::StdWallFunction(BaseContactRelation& wall_contact_relation)
            : BaseIntegrationFSI(wall_contact_relation), vonkar_(0.4187), E_(9.793),
            Cmu_(0.09), sigmak_(1.0), sigmaeps_(1.3), C1eps_(1.44), C2eps_(1.92),
            K_(*particles_->getVariableByName<Real>("TKE")),
            Eps_(*particles_->getVariableByName<Real>("Dissipation")),
            mu_t_(*particles_->getVariableByName<Real>("TurbulentViscosity")),
            Esp_p_(*this->particles_->template registerSharedVariable<Real>("DissipationNearWall")),
            dK_dt_(*particles_->getVariableByName<Real>("TKEChangeRate"))
            , distance_from_wall_(*particles_->getVariableByName<Vecd>("DistanceFromWall")), 
            bounds_(wall_contact_relation.getSPHBody()){}
        //=================================================================================================//
        void  StdWallFunction::interaction(size_t index_i, Real dt)
        {
            Real K_changerate_p = 0.0, Eps_p = 0.0;
            Vecd vel_in_wall = -vel_[index_i];
            Real p_in_wall = p_[index_i];
            Real rho_in_wall = rho_[index_i];
            Real K_in_wall = 0.0;
            Real resolution_ref = bounds_.getSPHBodyResolutionRef();
            Real yp = distance_from_wall_[index_i].norm();
            /*BoundingBox bounds = bounds_.getSPHSystemBounds();
            Real first_max = bounds.first_[1];
            Real second_max = bounds.second_[1];
            Real nearwall = 0.1 * (std::max(first_max, second_max));*/

                //ParticleConfiguration* contact_neighborhood = contact_configuration_[index_i];
                for (size_t k = 0; k < contact_configuration_.size(); ++k)
                {
                    StdLargeVec<Vecd> &n_k = *(wall_n_[k]);
                    wall_pos_.push_back(&(contact_particles_[k]->pos_));
                    StdLargeVec<Vecd> &pos_k = *(wall_pos_[k]);
                    StdLargeVec<Real> &Vol_k = *(wall_Vol_[k]);
                    Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
                    Real K_prod_p = 0.0, K_advection_p = 0.0;
                    Vecd K_laplacian_p = Vecd::Zero();

                    for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
                    {
                        size_t index_j = wall_neighborhood.j_[n];
                        Vecd &e_ij = wall_neighborhood.e_ij_[n];
                        Real dW_ij = wall_neighborhood.dW_ij_[n];
                        Real r_ij = wall_neighborhood.r_ij_[n];
                        Real W_ij = wall_neighborhood.W_ij_[n];


                        Vecd veltangential = (vel_[index_i] - n_k[index_j].dot(vel_[index_i]) * (n_k[index_j]));
                        
                        //Vecd poswall = pos_k[index_j];
                        Real ystar = (rho_[index_i] * std::pow(Cmu_, 0.25) * std::sqrt(K_[index_i]) * yp) / (fluid_.ReferenceViscosity()); //corect yp
                        Real u_star = 0.0, wall_shear_stress = 0.0, mu_eff_wall = 0.0, mu_eff_wall_1 = 0.0, friction_velocity = 0.0, yplus = 0.0;

                        if (yp <= resolution_ref)
                        {
                            if (ystar >= 11.225) // Particle is in log law region
                            {
                                u_star = (1.0 / vonkar_) * std::log(E_ * ystar);
                                wall_shear_stress = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * rho_[index_i]) / (u_star);
                                mu_eff_wall = wall_shear_stress * yp / (veltangential.norm());//equal to mu_t_
                                friction_velocity = std::sqrt(wall_shear_stress / rho_[index_i]);
                                yplus = (yp * friction_velocity * rho_[index_i]) / (fluid_.ReferenceViscosity());
                                mu_eff_wall_1 = (yp * friction_velocity * vonkar_ * rho_[index_i]) / (std::log(E_ * yplus));
 
                                K_prod_p = std::pow(wall_shear_stress, 2.0) / (vonkar_ * rho_[index_i] * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * yp);
                                Eps_p = (std::pow(Cmu_, 3.0 / 4.0) * std::pow(K_[index_i], 1.5)) / (vonkar_ * yp);

                            }
                            else if (ystar < 11.225) // Particle is in viscous sublayer
                            {
                                u_star = ystar; 
                                wall_shear_stress = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * rho_[index_i]) / (u_star);
                                mu_eff_wall = wall_shear_stress * yp / (veltangential.norm()); //equal to mu
                                friction_velocity = std::sqrt(wall_shear_stress / rho_[index_i]);
                                yplus = (yp * friction_velocity * rho_[index_i]) / (mu_eff_wall);
                                mu_eff_wall_1 = (yp * friction_velocity * vonkar_) / (rho_[index_i] * std::log(E_ * yplus));
                                //mu_t_[index_i] = mu_eff_wall - fluid_.ReferenceViscosity();
                                K_prod_p = 0.0;
                                Eps_p = (K_[index_i]) * 2.0 * fluid_.ReferenceViscosity() / (rho_[index_i] * yp * yp);
                            }
                            K_advection_p = dW_ij * Vol_k[index_j] * (rho_[index_i] * K_[index_i] * vel_[index_i] - rho_in_wall * K_in_wall * vel_in_wall).dot(e_ij);
                            //K_laplacian_p = 2 * ((mu_t_[index_i] / (sigmak_)+fluid_.ReferenceViscosity()) * ((K_[index_i] - K_in_wall) / r_ij * dW_ij * Vol_k[index_j])) * e_ij;
                            K_laplacian_p = 2 * ((mu_eff_wall) * ((K_[index_i] - K_in_wall) / r_ij * dW_ij * Vol_k[index_j])) * e_ij;

                            K_changerate_p += -K_advection_p + K_prod_p - Vol_k[index_j] * W_ij * (rho_in_wall * Eps_p) + K_laplacian_p.sum();
                        }
                    }
                }
            
            dK_dt_[index_i] += K_changerate_p;
            Esp_p_[index_i] += Eps_p;
            
        }
        //=================================================================================================//
        void StdWallFunction::update(size_t index_i, Real dt)
        {
            Real resolution_ref = bounds_.getSPHBodyResolutionRef();
            if (distance_from_wall_[index_i].norm() <= resolution_ref)
            {
                K_[index_i] += dK_dt_[index_i] * dt;
                Eps_[index_i] += Esp_p_[index_i];
            }
            
        }
        //=================================================================================================//
    }// namespace fluid_dynamics

}// namespace SPH
#endif // TURBULENCEMODEL_CPP