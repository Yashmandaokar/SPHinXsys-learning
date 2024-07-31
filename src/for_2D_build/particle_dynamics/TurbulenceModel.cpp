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
    Kepsprofiles::Kepsprofiles(BaseInnerRelation &inner_relation, const std::string &meshdatafull_path)
        : BaseIntegration<DataDelegateInner>(inner_relation),
          Kprof_(*this->particles_->template registerSharedVariable<Real>("TKEProfile")),
          Epsprof_(*this->particles_->template registerSharedVariable<Real>("DissipationProfile")),
          mu_tprof_(*this->particles_->template registerSharedVariable<Real>("TurblunetViscosityProfile"))
          ,vel_prof_(*this->particles_->template registerSharedVariable<Vecd>("VelocityProfile"))
    {
        x_coords_.resize(55322);
        y_coords_.resize(55322);
        tke_profile_.resize(55322);
        dissipation_profile_.resize(55322);
        mu_t_profile_.resize(55322);
        x_velocity_profile_.resize(55322);
        y_velocity_profile_.resize(55322);
        meshdata(meshdatafull_path);
        
        /*// Initialize y-coordinates
        y_coords_ = 
        {0.0, 0.0712782, 0.149024, 0.160953, 0.238582, 0.381846, 0.459302, 0.473112,
         0.550373, 0.679711, 0.756784, 0.794729, 0.871638, 0.976704, 1.05351,
         1.11448, 1.19118, 1.27309, 1.34967, 1.43298, 1.50948, 1.56883, 1.64529,
         1.75089, 1.82733, 1.86421, 1.94065, 2.0
        };

        // Initialize TKE values
        tke_profile_ = 
        {0.00755791, 0.00755791, 0.00734456, 0.00664912, 0.00621983, 0.00534681, 0.0049434,
         0.00422235, 0.00389449, 0.00330619, 0.00306458, 0.00267755, 0.00255695,
         0.00242772, 0.00243493, 0.00258291, 0.00271379, 0.0031172, 0.00336372,
         0.00396043, 0.00429262, 0.00501129, 0.00540814, 0.00626737, 0.00668813,
         0.00736458, 0.00757005, 0.00757005
        };

        dissipation_profile_ = 
        {0.00497487, 0.00497487, 0.00253135, 0.00131829, 0.000999881, 0.000632135, 0.000512339,
         0.00034863, 0.000289472, 0.000201208, 0.000170069, 0.00012531, 0.000112709,
         9.97106e-05, 0.000100505, 0.000115575, 0.000129415, 0.000176676, 0.000208999,
         0.000300047, 0.000360984, 0.000527548, 0.00064874, 0.0010207, 0.00134396,
         0.00258197, 0.00507232, 0.00507232
        };
      
        // Initialize mut values
        mu_t_profile_ = 
        {0.00103339, 0.00103339, 0.00191788, 0.00301829, 0.00348218, 0.00407026, 0.00429276,
         0.00460242, 0.00471559, 0.00488937, 0.00497003, 0.00514908, 0.00522068,
         0.00531983, 0.00530916, 0.00519511, 0.00512162, 0.00494989, 0.00487234,
         0.00470476, 0.0045941, 0.0042843, 0.00405758, 0.00346349, 0.00299547,
         0.00189054, 0.00101679, 0.00101679
        };

        
         x_velocity_profile_ = 
        {0.74706, 0.74706, 0.82501, 0.913886, 0.946275, 0.99212, 1.01085, 1.04117,
         1.05361, 1.07379, 1.08155, 1.09304, 1.09655, 1.10008, 1.09988,
         1.09596, 1.0921, 1.0802, 1.07212, 1.05167, 1.03898, 1.00864,
         0.989741, 0.944279, 0.911466, 0.823599, 0.742894, 0.742894
        };

        y_velocity_profile_ = 
        {0.00075838, 0.00075838, 0.000478437, -1.70688e-05, 6.07162e-05, -7.78661e-06,
         -1.49494e-06, -2.2222e-05, -1.98024e-05, -1.95972e-05, -1.62276e-05,
         -6.52933e-06, 3.66527e-06, -8.62617e-06, -7.59906e-06, 3.91547e-05,
         3.66409e-05, 4.76686e-05, 4.94376e-05, 5.36917e-05, 5.11559e-05,
         5.82096e-05, 5.03552e-05, 9.49423e-05, 2.59572e-05, 0.000433811, 0.000669656, 0.000669656
        };*/
    } 
    //=================================================================================================//
    void Kepsprofiles::meshdata(const std::string &full_path)
    {
        std::ifstream meshdatafile; /*!< \brief File object for the Ansys ASCII mesh file. */
        meshdatafile.open(full_path);
        if (meshdatafile.fail())
        {
            std::cout << "Error:Check if the file exists." << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        }
        
        std::string text_line;
        size_t i = 0;
        
            while (getline(meshdatafile, text_line))
            {
            if (i == 2151)
            {
                Real u = 1;
            }
                std::string x_string_copy = text_line;
                size_t first_devide_position = text_line.find_first_of(" ", 0);
                std::string x_index_string = x_string_copy.erase(first_devide_position);
                Real x_index_decimal = stod(x_index_string);
                x_coords_[i] = x_index_decimal;

                /*--- find the node2 between two cells---*/
                std::string y_string = text_line;
                std::string xvel_string = text_line;
                std::string yvel_string = text_line;
                std::string tke_string = text_line;
                std::string eps_string = text_line;
                std::string mut_string = text_line;

                y_string = y_string.erase(0, first_devide_position + 1);
                xvel_string = y_string;
                size_t second_devide_position = y_string.find_first_of(" ", 0);
                y_string.erase(second_devide_position);
                Real y_index_decimal = stod(y_string);
                y_coords_[i] = y_index_decimal;

                xvel_string = xvel_string.erase(0, second_devide_position + 1);
                yvel_string = xvel_string;
                size_t third_devide_position = xvel_string.find_first_of(" ", 0);
                xvel_string.erase(third_devide_position);
                Real xvel_index_decimal = stod(xvel_string);
                x_velocity_profile_[i] = xvel_index_decimal;

                yvel_string = yvel_string.erase(0, third_devide_position + 1);
                tke_string = yvel_string;
                size_t fourth_devide_position = yvel_string.find_first_of(" ", 0);
                yvel_string.erase(fourth_devide_position);
                Real yvel_index_decimal = stod(yvel_string);
                y_velocity_profile_[i] = yvel_index_decimal;

                tke_string = tke_string.erase(0, fourth_devide_position + 1);
                eps_string = tke_string;
                size_t fifth_devide_position = tke_string.find_first_of(" ", 0);
                tke_string.erase(fifth_devide_position);
                Real tke_index_decimal = stod(tke_string);
                tke_profile_[i] = tke_index_decimal;

                eps_string = eps_string.erase(0, fifth_devide_position + 1);
                mut_string = eps_string;
                size_t sixth_devide_position = eps_string.find_first_of(" ", 0);
                eps_string.erase(sixth_devide_position);
                Real eps_index_decimal = stod(eps_string);
                dissipation_profile_[i] = eps_index_decimal;

                mut_string = mut_string.erase(0, sixth_devide_position + 1);
                //size_t last_devide_position = mut_string.find_first_of(" ", 0);
                //mut_string.erase(last_devide_position);
                Real mut_index_decimal = stod(mut_string);
                mu_t_profile_[i] = mut_index_decimal;
                i++;
            }
    }
    //=================================================================================================//
    void Kepsprofiles::update(size_t index_i, Real dt)
    {
            for (size_t i = 0; i <= 55321; i++)
        {
                Real xdev = pos_[index_i][0] - x_coords_[i];
                Real ydev = pos_[index_i][1] - y_coords_[i];
                
                if (std::abs(xdev) < 1e-7 && std::abs(ydev) < 1e-7)
                {
                    if (index_i == 23070)
                    {
                        Vecd cellcentre = pos_[index_i];
                        Real u = 1;
                    }
                    vel_[index_i][0] = x_velocity_profile_[i];
                    vel_[index_i][1] = y_velocity_profile_[i];
                    Kprof_[index_i] = tke_profile_[i];
                    Epsprof_[index_i] = dissipation_profile_[i];
                    mu_tprof_[index_i] = mu_t_profile_[i];
                }
        }
        /*
        // Check if y is within the range of y_coords
        if (pos_[index_i][1] < y_coords_.front() || pos_[index_i][1] > y_coords_.back())
        {
            throw std::out_of_range("Interpolation point y is out of range.");
        }

        // Perform linear interpolation
        auto it = std::lower_bound(y_coords_.begin(), y_coords_.end(), pos_[index_i][1]);
        size_t idx = it - y_coords_.begin();
        if (index_i == 4827)
        {
            Real posi = pos_[index_i][1];
            Real x = 1.0;
        }
        if (idx == 0)
        {
            Kprof_[index_i] = tke_profile_[0];
            Epsprof_[index_i] = dissipation_profile_[0];
            mu_tprof_[index_i] = mu_t_profile_[0];
            vel_prof_[index_i][0] = x_velocity_profile_[0];
            vel_prof_[index_i][1] = y_velocity_profile_[0];
        }
        else if (idx == y_coords_.size())
        {
            Kprof_[index_i] = tke_profile_.back();
            Epsprof_[index_i] = dissipation_profile_.back();
            mu_tprof_[index_i] = mu_t_profile_.back();
            vel_prof_[index_i][0] = x_velocity_profile_.back();
            vel_prof_[index_i][1] = y_velocity_profile_.back();
        }
        else
        {
            double y1 = y_coords_[idx - 1];
            double y2 = y_coords_[idx];
            double tke1 = tke_profile_[idx - 1];
            double tke2 = tke_profile_[idx];
            Kprof_[index_i] = tke1 + (pos_[index_i][1] - y1) * (tke2 - tke1) / (y2 - y1);
            double eps1 = dissipation_profile_[idx - 1];
            double eps2 = dissipation_profile_[idx];
            Epsprof_[index_i] = eps1 + (pos_[index_i][1] - y1) * (eps2 - eps1) / (y2 - y1);
            double mut1 = mu_t_profile_[idx - 1];
            double mut2 = mu_t_profile_[idx];
            mu_tprof_[index_i] = mut1 + (pos_[index_i][1] - y1) * (mut2 - mut1) / (y2 - y1);
            double xvel1 = x_velocity_profile_[idx - 1];
            double xvel2 = x_velocity_profile_[idx];
            vel_prof_[index_i][0] = xvel1 + (pos_[index_i][1] - y1) * (xvel2 - xvel1) / (y2 - y1);
            double yvel1 = y_velocity_profile_[idx - 1];
            double yvel2 = y_velocity_profile_[idx];
            vel_prof_[index_i][1] = yvel1 + (pos_[index_i][1] - y1) * (yvel2 - yvel1) / (y2 - y1);
            Real X = 2.3;
        }*/
    }
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
            Eps_lap_(*this->particles_->template registerSharedVariable<Real>("DissipationLaplacian")),
            Eps_prodscalar_(*this->particles_->template registerSharedVariable<Real>("DissipationProdscalar")),
            Eps_scalar_(*this->particles_->template registerSharedVariable<Real>("DissipationScalar")),
            Tau_wall_(*this->particles_->template registerSharedVariable<Real>("WallShearStress")),
            Cmu_(0.09), sigmak_(1.0), sigmaeps_(1.3), C1eps_(1.44), C2eps_(1.92), 
            K_(*particles_->getVariableByName<Real>("TKE")),
            Eps_(*particles_->getVariableByName<Real>("Dissipation")),
            mu_t_(*particles_->getVariableByName<Real>("TurbulentViscosity"))
            ,Kprof_(*particles_->getVariableByName<Real>("TKEProfile")),
            Epsprof_(*particles_->getVariableByName<Real>("DissipationProfile")),
            mu_tprof_(*particles_->getVariableByName<Real>("TurblunetViscosityProfile")),
            vel_prof_(*particles_->getVariableByName<Vecd>("VelocityProfile"))
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
              walladjacentcellflag_(*this->particles_->template registerSharedVariable<Real>("FlagForWallAdjacentCells"))
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
        walladjacentindex_.resize(particles_->particles_bound_);
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
                walladjacentcellflag_[index_real] = 1.0; 
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
            walladjacentcellflag_(*particles_->getVariableByName<Real>("FlagForWallAdjacentCells")),
            strain_rate_(*this->particles_->template registerSharedVariable<Real>("StrainRate")),
            dudx_(*this->particles_->template registerSharedVariable<Real>("dudx")),
            dudy_(*this->particles_->template registerSharedVariable<Real>("dudy")),
            dvdx_(*this->particles_->template registerSharedVariable<Real>("dvdx")),
            dvdy_(*this->particles_->template registerSharedVariable<Real>("dvdy")),
            vel_gradient_mat_(*particles_->getVariableByName<Matd>("VelocityGradient"))
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
            mu_t_[index_i] = rho_[index_i] * Cmu_ * ((Kprof_[index_i] * Kprof_[index_i]) / (Epsprof_[index_i]));

            if (walladjacentcellflag_[index_i] == 1.0)
            {
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    Real mu_t_avg = (2.0 * mu_tprof_[index_i] * mu_tprof_[index_j]) / (mu_tprof_[index_i] + mu_tprof_[index_j] + TinyReal);

                    K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (Kprof_[index_i] - Kprof_[index_j]) * e_ij).dot(vel_[index_i]);
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_avg / sigmak_) * (Kprof_[index_i] - Kprof_[index_j]) / (r_ij));
                    
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
                nearwallquantities(index_i);
                K_prod_[index_i] = K_prod_p_[index_i];
                Eps_[index_i] = Eps_p_[index_i];
                strain_rate_[index_i] = std::sqrt(K_prod_[index_i] / mu_tprof_[index_i]);
                dudx_[index_i] = vel_gradient_mat_[index_i](0, 0);
                dudy_[index_i] = vel_gradient_mat_[index_i](0, 1);
                dvdx_[index_i] = vel_gradient_mat_[index_i](1, 0);
                dvdy_[index_i] = vel_gradient_mat_[index_i](1, 1);

                dK_dt_[index_i] = K_adv_[index_i] + K_prod_[index_i] - rho_[index_i] * Epsprof_[index_i] + K_lap_[index_i]; 
               
            }
            else
            {
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    Real mu_t_avg = (2.0 * mu_tprof_[index_i] * mu_tprof_[index_j]) / (mu_tprof_[index_i] + mu_tprof_[index_j] + TinyReal);

                    K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (Kprof_[index_i] - Kprof_[index_j]) * e_ij).dot(vel_[index_i]);
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_avg / sigmak_) * (Kprof_[index_i] - Kprof_[index_j]) / (r_ij));
                    vel_matrix = (vel_[index_i] - vel_[index_j]) * e_ij.transpose();
                    vel_gradient_mat = dW_ij * Vol_[index_j] * vel_matrix;
                    vel_gradient_mat_[index_i] += vel_gradient_mat;
                    // vel_gradient_mat_transpose = dW_ij * Vol_[index_j] * vel_matrix.transpose();
                    strain_tensor = 0.5 * (vel_gradient_mat + vel_gradient_mat.transpose());
                    strain_rate_modulus = 2.0 * strain_tensor.array() * strain_tensor.array();
                    
                    K_prod = mu_t_avg * strain_rate_modulus;
                    K_prod_iso = (2.0 / 3.0 * rho_[index_i] * Kprof_[index_i] * Matd::Identity()).array() * vel_gradient_mat.array();
                    K_prod_total = K_prod - K_prod_iso;
                    Real Ktot = K_prod_total.sum();
                    K_prod_[index_i] += K_prod_total.sum();
                    Eps_sum_[index_i] += -rho_[index_i] * Epsprof_[index_i];

                    if (index_i == 23070)
                    {
                        Vecd veli = vel_[index_i];
                        Vecd velj = vel_[index_j];
                        Matd velgrad = vel_gradient_mat_[index_i];
                        Real sr = strain_rate_[index_i];
                        Real Kprodtot = K_prod_total.sum();
                        Real volj = Vol_[index_j];
                        Real dudx = vel_gradient_mat_[index_i](0, 0);
                        Real dudy = vel_gradient_mat_[index_i](0, 1);
                        Real dvdx = vel_gradient_mat_[index_i](1, 0);
                        Real dvdy = vel_gradient_mat_[index_i](1, 1);
                        Real k = 1.0;

                    }
                }
                //strain_rate_[index_i] = sqrt(K_prod_[index_i]/mu_tprof_[index_i);
                Matd total_strain_tensor = 0.5 * (vel_gradient_mat_[index_i] + vel_gradient_mat_[index_i].transpose());
                Matd total_strain_modulus = 2.0 * total_strain_tensor.array() * total_strain_tensor.array();
                strain_rate_[index_i] = sqrt(total_strain_modulus.sum());
                dudx_[index_i] = vel_gradient_mat_[index_i](0, 0);
                dudy_[index_i] = vel_gradient_mat_[index_i](0, 1);
                dvdx_[index_i] = vel_gradient_mat_[index_i](1, 0);
                dvdy_[index_i] = vel_gradient_mat_[index_i](1, 1);
                dK_dt_[index_i] = K_adv_[index_i] + K_prod_[index_i] - rho_[index_i] * Epsprof_[index_i] + K_lap_[index_i];
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
            walladjacentcellflag_(*particles_->getVariableByName<Real>("FlagForWallAdjacentCells"))
        {}
        //=================================================================================================//
        
        void KEpsilonStd2ndHalf::interaction(size_t index_i, Real dt)
        {
            mu_t_[index_i] = rho_[index_i] * Cmu_ * ((Kprof_[index_i] * Kprof_[index_i]) / (Epsprof_[index_i]));
            Real Eps_changerate = 0.0;
            Eps_adv_[index_i] = 0.0, Eps_lap_[index_i] = 0.0, Eps_prodscalar_[index_i] = 0.0, Eps_scalar_[index_i] = 0.0;
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            if (walladjacentcellflag_[index_i] != 1)
            {
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    Real mu_t_avg = (2.0 * mu_tprof_[index_i] * mu_tprof_[index_j]) / (mu_tprof_[index_i] + mu_tprof_[index_j] + TinyReal);

                    Eps_adv_[index_i] = -(dW_ij * Vol_[index_j] * rho_[index_i] * (Epsprof_[index_i] - Epsprof_[index_j]) * e_ij).dot(vel_[index_i]);
                    Eps_lap_[index_i] = 2.0 * dW_ij * Vol_[index_j] * (fluid_.ReferenceViscosity() + mu_t_avg / sigmaeps_) * ((Epsprof_[index_i] - Epsprof_[index_j]) / (r_ij));
                    Eps_prodscalar_[index_i] = C1eps_ * (Epsprof_[index_i] / (Kprof_[index_i])) * K_prod_[index_i];
                    Eps_scalar_[index_i] = C2eps_ * rho_[index_i] * (Epsprof_[index_i] * Epsprof_[index_i]) / (Kprof_[index_i]);
                    Eps_changerate += Eps_adv_[index_i] + Eps_prodscalar_[index_i] - Eps_scalar_[index_i] + Eps_lap_[index_i];

                    if (index_i == 9560)
                    {
                        Vecd veli = vel_[index_i];
                        Vecd velj = vel_[index_j];
                        Real epsi = Epsprof_[index_i];
                        Real epsj = Epsprof_[index_j];
                        Real Voli = Vol_[index_i];
                        Real rhoi = rho_[index_i];
                        Real muti = mu_tprof_[index_i];
                        Real epsadv = Eps_adv_[index_i];
                        Real X = 1.0;
                    }
                }
                dEps_dt_[index_i] = Eps_changerate;
                
            }  
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
            ystar_[index_i] = (rho_[index_i] * std::pow(Cmu_, 0.25) * std::pow(Kprof_[index_i], 0.5) * yp) / (fluid_.ReferenceViscosity());
            Real u_star, mu_eff_wall, friction_velocity;
            Vecd veltangential = (vel_[index_i] - e_ij.dot(vel_[index_i]) * (e_ij));
            if (index_i == 7614)
            {
                Real h = 1.3;
            }
            if (ystar_[index_i] >= 11.225)
            {
                u_star = (1.0 / vonkar_) * std::log(E_ * ystar_[index_i]);
                Real Mut = fluid_.ReferenceViscosity() * ((ystar_[index_i]) / (1 / vonkar_ * std::log(E_ * ystar_[index_i])) - 1.0);
                // Tau_wall_[index_i] = (mu_tprof_[index_i] + fluid_.ReferenceViscosity()) * (veltangential.norm() / yp_[index_i]);
                Tau_wall_[index_i] = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(Kprof_[index_i], 0.5) * rho_[index_i]) / (u_star);
                vel_gradient_mat_[index_i](0, 1) = (Tau_wall_[index_i] /(mu_tprof_[index_i] + fluid_.ReferenceViscosity())) * -1.0 * e_ij[1];
                //mu_eff_wall = Tau_wall_[index_i] * yp / (veltangential.norm());
                // mu_eff_wall_1 = wall_shear_stress_1 * yp_[index_i] / (veltangential.norm());
                friction_velocity = std::sqrt(Tau_wall_[index_i] / rho_[index_i]);
                // friction_velocity_1 = std::sqrt(wall_shear_stress_1 / rho_[index_i]);

                K_prod_p_[index_i] = std::pow(Tau_wall_[index_i], 2.0) / (vonkar_ * rho_[index_i] * std::pow(Cmu_, 0.25) * std::pow(Kprof_[index_i], 0.5) * yp);
                Eps_p_[index_i] = (std::pow(Cmu_, 3.0 / 4.0) * std::pow(Kprof_[index_i], 1.5)) / (vonkar_ * yp);
            }
            else if (ystar_[index_i] < 11.225)
            {
                u_star = ystar_[index_i];
                // wall_shear_stress_1 = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * rho_[index_i]) / (u_star);
                // mu_tprof_[index_i] = 0.0;
                Tau_wall_[index_i] = fluid_.ReferenceViscosity() * veltangential.norm() / yp;
                mu_eff_wall = Tau_wall_[index_i] * yp / (veltangential.norm());
                // mu_eff_wall_1 = wall_shear_stress_1 * yp_[index_i] / (veltangential.norm());
                friction_velocity = std::sqrt(Tau_wall_[index_i] / rho_[index_i]);
                // friction_velocity_1 = std::sqrt(wall_shear_stress_1 / rho_[index_i]);

                K_prod_p_[index_i] = 0.0;
                Eps_p_[index_i] = (Kprof_[index_i] * 2.0 * fluid_.ReferenceViscosity()) / (rho_[index_i] * yp * yp);
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