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
            : BaseIntegration<DataDelegateInner>(inner_relation), ghost_creator_(ghost_creator),
            mom_(*this->particles_->template registerSharedVariable<Vecd>("Momentum")),
            dmom_dt_(*this->particles_->template registerSharedVariable<Vecd>("MomentumChangeRate")),
            dmass_dt_(*this->particles_->template registerSharedVariable<Real>("MassChangeRate")),
            K_prod_p_(*this->particles_->template registerSharedVariable<Real>("TKEProductionInWallAdjCell")),
            K_prod_(*this->particles_->template registerSharedVariable<Real>("TKEProduction")),
            Eps_p_(*this->particles_->template registerSharedVariable<Real>("DissipationRateInWallAdjCell")),
            Eps_sum_(*this->particles_->template registerSharedVariable<Real>("DissipationSum")),
            K_adv_(*this->particles_->template registerSharedVariable<Real>("TKEAdvection")),
            K_lap_(*this->particles_->template registerSharedVariable<Real>("TKELaplacian")),
            Eps_adv_(*this->particles_->template registerSharedVariable<Real>("DissipationAdvection")),
            Tau_wall_(*this->particles_->template registerSharedVariable<Real>("WallShearStress")),
            Eps_lap_(*this->particles_->template registerSharedVariable<Real>("DissipationLaplacian")),
            Eps_prodscalar_(*this->particles_->template registerSharedVariable<Real>("DissipationProdscalar")),
            Eps_scalar_(*this->particles_->template registerSharedVariable<Real>("DissipationScalar")),
            Cmu_(0.09), sigmak_(1.0), sigmaeps_(1.3), C1eps_(1.44), C2eps_(1.92),
            K_(*this->particles_->template getVariableDataByName<Real>("TKE")),
            Eps_(*this->particles_->template getVariableDataByName<Real>("Dissipation")),
            mu_t_(*this->particles_->template getVariableDataByName<Real>("TurblunetViscosity"))
            {}
        //=================================================================================================//
            /*
            WallAdjacentCells::WallAdjacentCells(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator)
            : BaseTurbulence(inner_relation, ghost_creator), yp_(*this->particles_->template registerSharedVariable<Real>("WallNormalDistance")),
              walladjacentcellflag_(*this->particles_->template registerSharedVariable<Real>("FlagForWallAdjacentCells")),
              wallnormal_(*this->particles_->template registerSharedVariable<Vecd>("WallNormal")), ymax_(0.0), 
              bounds_(inner_relation.getSPHBody())
        {
            walladjacentcellyp();
        }*/ 
        //=================================================================================================//
        
        WallAdjacentCells::WallAdjacentCells(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator)
            : BaseTurbulence(inner_relation, ghost_creator),
              walladjacentcellflag_(*this->particles_->template registerSharedVariable<Real>("FlagForWallAdjacentCells")),
              wallfacearea_(*this->particles_->template registerSharedVariable<Real>("WallFaceArea"))
              {
                walladjacentcellyp();
              }
        //=================================================================================================//
        /*void WallAdjacentCells::walladjacentcellyp()
        {
            walladjacentindex_.resize(particles_-> particles_bound_);
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
                    wallghostindex_[index_real] = ghost_index;

                    if ((pos_[index_real] - pos_[ghost_index]).dot(walleij_[index_real]) > ymax_)
                    {
                        ymax_ = (pos_[index_real] - pos_[ghost_index]).dot(walleij_[index_real]);

                    }
                    
                }
            }
        }*/
        //=================================================================================================//
         void WallAdjacentCells::walladjacentcellyp()
        {
        walladjacentindex_.resize(particles_->ParticlesBound());
                wallghostindex_.resize(particles_->ParticlesBound());
        walleij_.resize(particles_->ParticlesBound());
        Real boundary_type = 3;

        if (!ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type].empty())
        {
            for (size_t ghost_number = 0; ghost_number != ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type].size(); ++ghost_number)
            {

                size_t ghost_index = ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type][ghost_number];
                Real index_real = ghost_creator_.each_boundary_type_contact_real_index_[boundary_type][ghost_number];
                walleij_[index_real] = ghost_creator_.each_boundary_type_with_all_ghosts_eij_[boundary_type][ghost_number];
                walladjacentindex_[index_real] = index_real;
                wallghostindex_[index_real] = ghost_index;
                walladjacentcellflag_[index_real] = 1.0; 
                wallfacearea_[index_real] = ghost_creator_.each_boundary_type_face_area_[boundary_type][ghost_number];
            }
        }
        }
        //=================================================================================================//
        /*
        *  void WallAdjacentCells::update(size_t index_i, Real dt)
        {
         Vecd wallnormal = Vecd::Zero();
            Vecd lower_wall = {pos_[index_i][0], 0.0};
            Vecd upper_wall = {pos_[index_i][0], 2.0};
            Vecd lower_wall_normal = {0.0, 1.0};
            Vecd upper_wall_normal = {0.0, -1.0};
            BoundingBox bounds = bounds_.getSPHSystemBounds();
            Real channelheight = bounds.second_[1];
            Real halfwidth = 0.5 * channelheight;

            bool lower_wall_condition = ((pos_[index_i] - lower_wall).dot(lower_wall_normal) <= 0.3 * halfwidth);
            bool upper_wall_condition = ((pos_[index_i] - upper_wall).dot(upper_wall_normal) <= 0.3 * halfwidth);

            if (lower_wall_condition)
            {
                yp_[index_i] = (pos_[index_i] - lower_wall).dot(lower_wall_normal);
                walladjacentcellflag_[index_i] = (1);
                wallnormal_[index_i] = lower_wall_normal;
            }
            else if (upper_wall_condition)
            {
                yp_[index_i] = (pos_[index_i] - upper_wall).dot(upper_wall_normal);
                walladjacentcellflag_[index_i] = (1);    
                wallnormal_[index_i] = upper_wall_normal;
            }
        }
        */ 
        
        //=================================================================================================//

        KEpsilonStd1stHalf::KEpsilonStd1stHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh& ghost_creator) 
            : StdWallFunctionFVM(inner_relation, ghost_creator),
            dK_dt_(*this->particles_->template registerSharedVariable<Real>("TKEChangeRate")),
            walladjacentcellflag_(*this->particles_->template getVariableDataByName<Real>("FlagForWallAdjacentCells")),
            strain_rate_(*this->particles_->template registerSharedVariable<Real>("StrainRate")),
            dudx_(*this->particles_->template registerSharedVariable<Real>("dudx")),
            dudy_(*this->particles_->template registerSharedVariable<Real>("dudy")),
            dvdx_(*this->particles_->template registerSharedVariable<Real>("dvdx")),
            dvdy_(*this->particles_->template registerSharedVariable<Real>("dvdy")),
            vel_gradient_mat_(*this->particles_->template getVariableDataByName<Matd>("VelocityGradient"))
        {}
        //=================================================================================================//
        void KEpsilonStd1stHalf::interaction(size_t index_i, Real dt)
        {
            
            Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            Matd K_prod = Matd::Zero(), K_prod_iso = Matd::Zero(), K_prod_total = Matd::Zero();
            Matd vel_matrix = Matd::Zero(), vel_gradient_mat = Matd::Zero();
            Matd strain_tensor = Matd::Zero(), strain_rate_modulus = Matd::Zero();
            Real Kprodtot = 0.0;
            K_prod_[index_i] = 0.0, K_adv_[index_i] = 0.0, K_lap_[index_i] = 0.0, strain_rate_[index_i] = 0.0;
            Eps_sum_[index_i] = 0.0, dudx_[index_i] = 0.0, dudy_[index_i] = 0.0, dvdx_[index_i] = 0.0, dvdy_[index_i] = 0.0;
            vel_gradient_mat_[index_i] = Matd::Zero();
            mu_t_[index_i] = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));

            if (walladjacentcellflag_[index_i] == 1.0)
            {
                nearwallquantities(index_i);
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);

                    K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (K_[index_i] - K_[index_j]) * vel_[index_i]).dot(e_ij);
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_avg / sigmak_) * (K_[index_i] - K_[index_j]) / (r_ij));
                    
                   /*
                   vel_matrix = (vel_[index_i] - vel_[index_j]) * e_ij.transpose();
                   vel_gradient_mat = -dW_ij * Vol_[index_j] * vel_matrix;
                   vel_gradient_mat_[index_i] += vel_gradient_mat;
                   Matd velgrad = vel_gradient_mat_[index_i];
                   if (index_i == 773)
                   {
                        Real x = 1.0;
                   }*/    
                }
                K_prod_[index_i] = K_prod_p_[index_i];
                Eps_[index_i] = Eps_p_[index_i];
                strain_rate_[index_i] = std::sqrt(K_prod_[index_i] / mu_t_[index_i]);
                dudx_[index_i] = vel_gradient_mat_[index_i](0, 0);
                dudy_[index_i] = vel_gradient_mat_[index_i](0, 1);
                dvdx_[index_i] = vel_gradient_mat_[index_i](1, 0);
                dvdy_[index_i] = vel_gradient_mat_[index_i](1, 1);

                dK_dt_[index_i] = K_adv_[index_i] + K_prod_[index_i] - rho_[index_i] * Eps_[index_i] + K_lap_[index_i]; 
               
            }
            else
            {
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);

                    K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (K_[index_i] - K_[index_j]) * vel_[index_i]).dot(e_ij);
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_avg / sigmak_) * (K_[index_i] - K_[index_j]) / (r_ij));
                    vel_matrix = (vel_[index_i] - vel_[index_j]) * e_ij.transpose();
                    vel_gradient_mat = dW_ij * Vol_[index_j] * vel_matrix;
                    vel_gradient_mat_[index_i] += vel_gradient_mat;
                    /*
                    strain_tensor = 0.5 * (vel_gradient_mat + vel_gradient_mat.transpose());
                    strain_rate_modulus = 2.0 * strain_tensor.array() * strain_tensor.array();
                    
                    K_prod = mu_t_[index_i] * strain_rate_modulus;
                    K_prod_iso = (2.0 / 3.0 * rho_[index_i] * K_[index_i] * Matd::Identity()).array() * vel_gradient_mat.array();
                    K_prod_total = K_prod - K_prod_iso;
                    Real Ktot = K_prod_total.sum();
                    K_prod_[index_i] += K_prod.sum();
                    */  
                }
                strain_tensor = 0.5 * (vel_gradient_mat_[index_i] + vel_gradient_mat_[index_i].transpose());
                strain_rate_modulus = 2.0 * strain_tensor.array() * strain_tensor.array();

                K_prod = mu_t_[index_i] * strain_rate_modulus;
               
                K_prod_[index_i] += K_prod.sum();
                strain_rate_[index_i] = sqrt(strain_rate_modulus.sum());
                //strain_rate_[index_i] = sqrt(K_prod_[index_i]/mu_tprof_[index_i);
                //Matd total_strain_tensor = 0.5 * (vel_gradient_mat_[index_i] + vel_gradient_mat_[index_i].transpose());
                //Matd total_strain_modulus = 2.0 * total_strain_tensor.array() * total_strain_tensor.array();
                //strain_rate_[index_i] = sqrt(total_strain_modulus.sum());
                dudx_[index_i] = vel_gradient_mat_[index_i](0, 0);
                dudy_[index_i] = vel_gradient_mat_[index_i](0, 1);
                dvdx_[index_i] = vel_gradient_mat_[index_i](1, 0);
                dvdy_[index_i] = vel_gradient_mat_[index_i](1, 1);
                dK_dt_[index_i] = K_adv_[index_i] + K_prod_[index_i] - rho_[index_i] * Eps_[index_i] + K_lap_[index_i];
            }
        }
        //=================================================================================================//
        void KEpsilonStd1stHalf::update(size_t index_i, Real dt)
        {
            K_[index_i] += (dK_dt_[index_i] / rho_[index_i]) * dt;
            if (K_[index_i] < 0.0)
            {
                Real Kdt = dK_dt_[index_i];
                Real k = K_[index_i];
                Real K = 1.0;
            }
        }
        //=================================================================================================//
        KEpsilonStd2ndHalf::KEpsilonStd2ndHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh& ghost_creator) 
            : BaseTurbulence(inner_relation, ghost_creator),
            dEps_dt_(*this->particles_->template registerSharedVariable<Real>("DissipationChangeRate")),
            walladjacentcellflag_(*this->particles_->template getVariableDataByName<Real>("FlagForWallAdjacentCells"))
        {}
        //=================================================================================================//
        
        void KEpsilonStd2ndHalf::interaction(size_t index_i, Real dt)
        {
            Real Eps_changerate = 0.0;
            Eps_adv_[index_i] = 0.0, Eps_lap_[index_i] = 0.0, Eps_prodscalar_[index_i] = 0.0, Eps_scalar_[index_i] = 0.0;
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            if (walladjacentcellflag_[index_i] != 1)
            {
                mu_t_[index_i] = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));
            }   
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j] + TinyReal);

                    Eps_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (Eps_[index_i] - Eps_[index_j]) * vel_[index_i]).dot(e_ij);
                    Eps_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * (fluid_.ReferenceViscosity() + mu_t_avg / sigmaeps_) * ((Eps_[index_i] - Eps_[index_j]) / (r_ij));
                    
                    Eps_changerate = Eps_adv_[index_i] + Eps_lap_[index_i];
                }
                Eps_prodscalar_[index_i] = C1eps_ * (Eps_[index_i] / (K_[index_i])) * K_prod_[index_i];
                Eps_scalar_[index_i] = -C2eps_ * rho_[index_i] * (Eps_[index_i] * Eps_[index_i]) / (K_[index_i]);
                dEps_dt_[index_i] = Eps_changerate + Eps_prodscalar_[index_i] + Eps_scalar_[index_i];
        }
        //=================================================================================================//
        void KEpsilonStd2ndHalf::update(size_t index_i, Real dt)
        {
            if (walladjacentcellflag_[index_i] != 1)
            {
                Eps_[index_i] += (dEps_dt_[index_i] / rho_[index_i]) * dt;
            }
            if (Eps_[index_i] < 0.0)
            {
                Real edt = dEps_dt_[index_i];
                Real e = Eps_[index_i];
                Real K = 1.0;
            }
        }
        //=================================================================================================//
         StdWallFunctionFVM::StdWallFunctionFVM(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator)
            : WallAdjacentCells(inner_relation, ghost_creator), ystar_(*this->particles_->template registerSharedVariable<Real>("Ystar")),
              vel_gradient_mat_(*this->particles_->template registerSharedVariable<Matd>("VelocityGradient")),
              vonkar_(0.4187), E_(9.793)
           {}
        //=================================================================================================//
        void StdWallFunctionFVM::nearwallquantities(size_t index_i)
        {
            size_t ghost_index = wallghostindex_[index_i];
            size_t index_real = walladjacentindex_[index_i];
            Vecd e_ij = walleij_[index_i];

            Real yp = (pos_[index_i] - pos_[ghost_index]).dot(e_ij);
            ystar_[index_i] = (rho_[index_i] * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * yp) / (fluid_.ReferenceViscosity());
            Real u_star, mu_eff_wall, friction_velocity;
            Vecd veltangential = (vel_[index_i] - e_ij.dot(vel_[index_i]) * (e_ij));
            if (index_i == 7614)
            {
                Real h = 1.3;
            }
            if (ystar_[index_i] >= 11.225)
            {
                u_star = (1.0 / vonkar_) * std::log(E_ * ystar_[index_i]);
                mu_t_[index_i] = fluid_.ReferenceViscosity() * ((ystar_[index_i]) / (1 / vonkar_ * std::log(E_ * ystar_[index_i])) - 1.0);
                Tau_wall_[index_i] = (mu_t_[index_i] + fluid_.ReferenceViscosity()) * (veltangential.norm() / yp);
                //Tau_wall_[index_i] = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(Kprof_[index_i], 0.5) * rho_[index_i]) / (u_star);
                vel_gradient_mat_[index_i](0, 1) = (Tau_wall_[index_i] / (mu_t_[index_i] + fluid_.ReferenceViscosity())) * e_ij[1];
                //mu_eff_wall = Tau_wall_[index_i] * yp / (veltangential.norm());
                // mu_eff_wall_1 = wall_shear_stress_1 * yp_[index_i] / (veltangential.norm());
                //friction_velocity = std::sqrt(Tau_wall_[index_i] / rho_[index_i]);
                // friction_velocity_1 = std::sqrt(wall_shear_stress_1 / rho_[index_i]);

                K_prod_p_[index_i] = std::pow(Tau_wall_[index_i], 2.0) / (vonkar_ * rho_[index_i] * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * yp);
                Eps_p_[index_i] = (std::pow(Cmu_, 3.0 / 4.0) * std::pow(K_[index_i], 1.5)) / (vonkar_ * yp);
            }
            else if (ystar_[index_i] < 11.225)
            {
                u_star = ystar_[index_i];
                Tau_wall_[index_i] = fluid_.ReferenceViscosity() * veltangential.norm() / yp;

                //K_prod_p_[index_i] = 0.0;
                Eps_p_[index_i] = (K_[index_i] * 2.0 * fluid_.ReferenceViscosity()) / (rho_[index_i] * yp * yp);
            }  
        }
        //=================================================================================================//
        /*
        StdWallFunctionFVM::StdWallFunctionFVM(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator)
            : BaseTurbulence(inner_relation, ghost_creator), walladjacentcellflag_(*particles_->getVariableByName<Real>("FlagForWallAdjacentCells"))
            , yp_(*particles_->getVariableByName<Real>("WallNormalDistance")),
              wallnormal_(*particles_->getVariableByName<Vecd>("WallNormal")), vonkar_(0.4187), E_(9.793)
        {}*/ 
        //=================================================================================================//
        /*
        void StdWallFunctionFVM::nearwallquantities(size_t index_i)
        {
            Real ystar = (rho_[index_i] * std::pow(Cmu_, 0.25) * std::pow(Kprof_[index_i], 0.5) * yp_[index_i]) / (fluid_.ReferenceViscosity());
            Real u_star, mu_eff_wall, friction_velocity;
            Vecd veltangential = (vel_[index_i] - wallnormal_[index_i].dot(vel_[index_i]) * (wallnormal_[index_i]));
            if (ystar >= 11.225)
            {
                u_star = (1.0 / vonkar_) * std::log(E_ * ystar);
                //mu_tprof_[index_i] = fluid_.ReferenceViscosity() * ((ystar * vonkar_) / (std::log(E_ * ystar)) - 1.0);
                //Tau_wall_[index_i] = (mu_tprof_[index_i] + fluid_.ReferenceViscosity()) * (veltangential.norm() / yp_[index_i]);
                Tau_wall_[index_i] = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(Kprof_[index_i], 0.5) * rho_[index_i]) / (u_star);
                mu_eff_wall = Tau_wall_[index_i] * yp_[index_i] / (veltangential.norm());
                //mu_eff_wall_1 = wall_shear_stress_1 * yp_[index_i] / (veltangential.norm());
                friction_velocity = std::sqrt(Tau_wall_[index_i] / rho_[index_i]);
                //friction_velocity_1 = std::sqrt(wall_shear_stress_1 / rho_[index_i]);

                K_prod_p_[index_i] = std::pow(Tau_wall_[index_i], 2.0) / (vonkar_ * rho_[index_i] * std::pow(Cmu_, 0.25) * std::pow(Kprof_[index_i], 0.5) * yp_[index_i]);
                Eps_p_[index_i] = (std::pow(Cmu_, 3.0 / 4.0) * std::pow(Kprof_[index_i], 1.5)) / (vonkar_ * yp_[index_i]);
            }
            else if (ystar < 11.225)
            {
                u_star = ystar;
                //wall_shear_stress_1 = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * rho_[index_i]) / (u_star);
                //mu_tprof_[index_i] = 0.0;
                Tau_wall_[index_i] = (fluid_.ReferenceViscosity()) * (veltangential.norm() / yp_[index_i]);
                mu_eff_wall = Tau_wall_[index_i] * yp_[index_i] / (veltangential.norm());
                //mu_eff_wall_1 = wall_shear_stress_1 * yp_[index_i] / (veltangential.norm());
                friction_velocity = std::sqrt(Tau_wall_[index_i] / rho_[index_i]);
                //friction_velocity_1 = std::sqrt(wall_shear_stress_1 / rho_[index_i]);

                K_prod_p_[index_i] = 0.0;
                Eps_p_[index_i] = (Kprof_[index_i] * 2.0 * fluid_.ReferenceViscosity()) / (rho_[index_i] * yp_[index_i] * yp_[index_i]);
            }
        }*/ 
        //=================================================================================================//
       
       /*
       StdWallFunctionFVM::StdWallFunctionFVM(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator)
            : WallAdjacentCells(inner_relation, ghost_creator), vonkar_(0.4187), E_(9.793)
            //dK_dt_(*particles_->getVariableByName<Real>("TKEChangeRate")),
            //K_prod_p_(*particles_->getVariableByName<Real>("TKEProductionInWallAdjCell")),
            //Esp_p_(*particles_->getVariableByName<Real>("DissipationRateInWallAdjCell")),
        {}*/ 
        //=================================================================================================//
        /*
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
                    u_star = (1.0 / vonkar_) * std::log(E_ * ystar);
                    wall_shear_stress_1 = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * rho_[index_i]) / (u_star);
                    mu_t_[index_i] = fluid_.ReferenceViscosity() * ((ystar * vonkar_) / (std::log(E_ * ystar)) - 1.0);
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
            else if (walladjacentcellflag_[index_i] != 1)
            {
                mu_t_[index_i] = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));
            }
        }*/ 
        //=================================================================================================//
         
        
    }// namespace fluid_dynamics

}// namespace SPH
#endif // TURBULENCEMODEL_CPP