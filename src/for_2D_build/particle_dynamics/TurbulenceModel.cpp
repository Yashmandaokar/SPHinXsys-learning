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
          mu_tprof_(*this->particles_->template registerSharedVariable<Real>("TurblunetViscosityProfile")),
          vel_prof_(*this->particles_->template registerSharedVariable<Vecd>("VelocityProfile"))
          {
            x_coords_.resize(55322);
            y_coords_.resize(55322);
            p.resize(55322);
            tke_profile_.resize(55322);
            dissipation_profile_.resize(55322);
            mu_t_profile_.resize(55322);
            x_velocity_profile_.resize(55322);
            y_velocity_profile_.resize(55322);
            meshdata(meshdatafull_path);
        
            /*// Initialize y-coordinates
            y_coords_ = 
            {0.0, 0.000906638, 0.0436064, 0.125532, 0.168244, 0.173645, 0.216377, 0.292655, 0.335425,
                 0.348909, 0.391733, 0.459711, 0.502612, 0.525215, 0.568193, 0.626927, 0.670031, 0.702929,
                 0.746176, 0.794683, 0.838112, 0.882519, 0.926142, 0.963439, 1.00726, 1.06419, 1.10824,
                 1.13361, 1.17783, 1.2472, 1.29161, 1.30568, 1.35006, 1.42921, 1.47351, 1.47909, 1.52317,
                 1.60814, 1.65201, 1.65263, 1.69631, 1.78265, 1.82615, 1.82724, 1.87062, 1.95672, 2.0
            };

            // Initialize TKE values
            tke_profile_ = 
            {0.00772659, 0.00772659, 0.00762386, 0.00742495, 0.00718538, 0.00658423, 0.00631289,
                 0.00583732, 0.00560646, 0.00517726, 0.00497198, 0.00457626, 0.00438325, 0.00402024,
                 0.00384986, 0.00352216, 0.00336787, 0.00308917, 0.00296905, 0.00275522, 0.00266621,
                 0.0025303, 0.00248508, 0.0024324, 0.0024268, 0.00245958, 0.00249356, 0.0026133, 0.00269515,
                 0.00289806, 0.00301497, 0.00329008, 0.00344532, 0.00377623, 0.00395351, 0.00432623,
                 0.00452243, 0.00492713, 0.00514189, 0.00558025, 0.00581205, 0.00629682, 0.00657169,
                 0.00717885, 0.00742023, 0.00762454, 0.00762454
            };

            dissipation_profile_ = 
            {0.0093729, 0.0093729, 0.00481526, 0.00260519, 0.00202288, 0.00136273, 0.0011503, 0.000862241,
                 0.000754245, 0.000592406, 0.000528937, 0.000425474, 0.000382187, 0.000311214, 0.000281628,
                 0.000230258, 0.000208244, 0.00017164, 0.000156999, 0.000132392, 0.000122708, 0.000108455,
                 0.000103895, 9.8652e-05, 9.81041e-05, 0.000101376, 0.000104773, 0.000117184, 0.000125926,
                 0.000148761, 0.000162674, 0.000197766, 0.000219196, 0.000269387, 0.000299301, 0.00037014,
                 0.000412735, 0.000515417, 0.000580174, 0.000741197, 0.00084657, 0.00113255, 0.0013427,
                 0.00199522, 0.00257039, 0.00475288, 0.00913
            };
      
            // Initialize mut values
            mu_t_profile_ = 
            {0.00057325, 0.00057325, 0.00108636, 0.00190454, 0.00229705, 0.00286315, 0.00311809,
                 0.00355665, 0.00375066, 0.00407214, 0.00420628, 0.00442988, 0.00452437, 0.00467399,
                 0.00473649, 0.00484894, 0.00490208, 0.0050039, 0.00505338, 0.0051605, 0.00521384, 0.00531298,
                 0.00534972, 0.00539765, 0.00540283, 0.0053707, 0.00534112, 0.00524508, 0.00519152, 0.00508124,
                 0.0050291, 0.00492611, 0.00487382, 0.00476412, 0.00470002, 0.00455088, 0.00445979, 0.00423908,
                 0.00410138, 0.00378107, 0.00359119, 0.00315085, 0.00289479, 0.00232466, 0.00192787,
                 0.00110081, 0.00110081
            };

        
             x_velocity_profile_ = 
            {0.670546, 0.670546, 0.750487, 0.837911, 0.870062, 0.914095, 0.932492, 0.962998, 0.976525,
                 0.999701, 1.00995, 1.02832, 1.03676, 1.05158, 1.05814, 1.07, 1.07538, 1.0845, 1.08832,
                 1.09482, 1.09752, 1.10145, 1.10275, 1.10423, 1.10442, 1.10351, 1.10255, 1.09913, 1.09675,
                 1.09065, 1.08702, 1.07812, 1.07291, 1.06112, 1.05448, 1.03955, 1.03125, 1.01272, 1.00218,
                 0.97879, 0.965476, 0.934753, 0.916344, 0.872097, 0.839975, 0.752345, 0.670546
            };

            y_velocity_profile_ = 
              {-0.000387598, -0.000387598, -0.00024092, 1.11873e-05, -3.13593e-05, 5.297e-06, 1.32937e-06,
             8.57787e-06, 5.22437e-06, 7.11216e-06, 6.93925e-06, 7.93779e-06, 6.33145e-06, 5.53202e-06,
             5.79597e-06, 6.11935e-06, 6.55074e-06, 9.80409e-06, 1.32577e-05, 1.85234e-05, 2.15542e-05,
             2.80438e-05, 3.10059e-05, 4.78082e-05, 5.0417e-05, 4.7223e-05, 5.35328e-05, 5.87922e-05,
             5.88573e-05, 6.07748e-05, 5.94737e-05, 5.92258e-05, 5.46051e-05, 5.01192e-05, 4.5167e-05,
             3.78497e-05, 3.11685e-05, 2.20534e-05, 1.56744e-05, 6.91562e-06, -1.32503e-07, -7.75231e-06,
             -1.34218e-05, -1.60569e-05, -3.71078e-05, 5.67885e-05, 5.67885e-05
            };*/
          } 
    //=================================================================================================//
    
    void Kepsprofiles::meshdata(const std::string &full_path)
    {
        std::ifstream meshdatafile; //!< \brief File object for the Ansys ASCII mesh file.
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
            if (i == 105)
            {
                Real u = 1;
            }
            std::string x_string_copy = text_line;
            size_t first_devide_position = text_line.find_first_of(" ", 0);
            std::string x_index_string = x_string_copy.erase(first_devide_position);
            Real x_index_decimal = stod(x_index_string);
            x_coords_[i] = x_index_decimal;

           // --- find the node2 between two cells---
            std::string y_string = text_line;
            std::string pr_string = text_line;
            std::string xvel_string = text_line;
            std::string yvel_string = text_line;
            std::string tke_string = text_line;
            std::string eps_string = text_line;
            std::string mut_string = text_line;

            y_string = y_string.erase(0, first_devide_position + 1);
            pr_string = y_string;
            size_t second_devide_position = y_string.find_first_of(" ", 0);
            y_string.erase(second_devide_position);
            Real y_index_decimal = stod(y_string);
            y_coords_[i] = y_index_decimal;

            pr_string = pr_string.erase(0, second_devide_position + 1);
            xvel_string = pr_string;
            size_t third_devide_position = pr_string.find_first_of(" ", 0);
            pr_string.erase(third_devide_position);
            Real pr_index_decimal = stod(pr_string);
            p[i] = pr_index_decimal;

            xvel_string = xvel_string.erase(0, third_devide_position + 1);
            yvel_string = xvel_string;
            size_t fourth_devide_position = xvel_string.find_first_of(" ", 0);
            xvel_string.erase(fourth_devide_position);
            Real xvel_index_decimal = stod(xvel_string);
            x_velocity_profile_[i] = xvel_index_decimal;

            yvel_string = yvel_string.erase(0, fourth_devide_position + 1);
            tke_string = yvel_string;
            size_t fifth_devide_position = yvel_string.find_first_of(" ", 0);
            yvel_string.erase(fifth_devide_position);
            Real yvel_index_decimal = stod(yvel_string);
            y_velocity_profile_[i] = yvel_index_decimal;

            tke_string = tke_string.erase(0, fifth_devide_position + 1);
            eps_string = tke_string;
            size_t sixth_devide_position = tke_string.find_first_of(" ", 0);
            tke_string.erase(sixth_devide_position);
            Real tke_index_decimal = stod(tke_string);
            tke_profile_[i] = tke_index_decimal;

            eps_string = eps_string.erase(0, sixth_devide_position + 1);
            mut_string = eps_string;
            size_t seventh_devide_position = eps_string.find_first_of(" ", 0);
            eps_string.erase(seventh_devide_position);
            Real eps_index_decimal = stod(eps_string);
            dissipation_profile_[i] = eps_index_decimal;

            mut_string = mut_string.erase(0, seventh_devide_position + 1);
            // size_t last_devide_position = mut_string.find_first_of(" ", 0);
            // mut_string.erase(last_devide_position);
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
                
                if (std::abs(xdev) < 1e-5 && std::abs(ydev) < 1e-5)
                {
                    if (index_i == 31284)
                    {
                        Vecd cellcentre = pos_[index_i];
                        Real u = 1;
                    }
                    p_[index_i] = p[i];
                    vel_[index_i][0] = x_velocity_profile_[i];
                    vel_[index_i][1] = y_velocity_profile_[i];
                    Kprof_[index_i] = tke_profile_[i];
                    Epsprof_[index_i] = dissipation_profile_[i];
                    mu_tprof_[index_i] = mu_t_profile_[i];
                }
        }
        
        /*// Check if y is within the range of y_coords
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
            Tau_wall_(*this->particles_->template registerSharedVariable<Real>("WallShearStress")),
            Eps_lap_(*this->particles_->template registerSharedVariable<Real>("DissipationLaplacian")),
            Eps_prodscalar_(*this->particles_->template registerSharedVariable<Real>("DissipationProdscalar")),
            Eps_scalar_(*this->particles_->template registerSharedVariable<Real>("DissipationScalar")),
            Cmu_(0.09), sigmak_(1.0), sigmaeps_(1.3), C1eps_(1.44), C2eps_(1.92),
            Kprof_(*this->particles_->template getVariableDataByName<Real>("TKEProfile")),
            Epsprof_(*this->particles_->template getVariableDataByName<Real>("DissipationProfile")),
            mu_tprof_(*this->particles_->template getVariableDataByName<Real>("TurblunetViscosityProfile")),
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
            //mu_tprof_[index_i] = rho_[index_i] * Cmu_ * ((Kprof_[index_i] * Kprof_[index_i]) / (Epsprof_[index_i]));

            if (walladjacentcellflag_[index_i] == 1.0)
            {
                nearwallquantities(index_i);
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    //Real mu_t_avg = (2.0 * mu_tprof_[index_i] * mu_tprof_[index_j]) / (mu_tprof_[index_i] + mu_tprof_[index_j] + TinyReal);

                    K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (Kprof_[index_i] - Kprof_[index_j]) * e_ij).dot(vel_[index_i]);
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_tprof_[index_i] / sigmak_) * (Kprof_[index_i] - Kprof_[index_j]) / (r_ij));
                    
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
                //Epsprof_[index_i] = Eps_p_[index_i];
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
                    //Real mu_t_avg = (2.0 * mu_tprof_[index_i] * mu_tprof_[index_j]) / (mu_tprof_[index_i] + mu_tprof_[index_j] + TinyReal);

                    K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (Kprof_[index_i] - Kprof_[index_j]) * e_ij).dot(vel_[index_i]);
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_tprof_[index_i] / sigmak_) * (Kprof_[index_i] - Kprof_[index_j]) / (r_ij));
                    vel_matrix = (vel_[index_i] - vel_[index_j]) * e_ij.transpose();
                    vel_gradient_mat = dW_ij * Vol_[index_j] * vel_matrix;
                    vel_gradient_mat_[index_i] += vel_gradient_mat;
                    // vel_gradient_mat_transpose = dW_ij * Vol_[index_j] * vel_matrix.transpose();
                    strain_tensor = 0.5 * (vel_gradient_mat + vel_gradient_mat.transpose());
                    strain_rate_modulus = 2.0 * strain_tensor.array() * strain_tensor.array();
                    
                    K_prod = mu_tprof_[index_i] * strain_rate_modulus;
                    K_prod_iso = (2.0 / 3.0 * rho_[index_i] * Kprof_[index_i] * Matd::Identity()).array() * vel_gradient_mat.array();
                    K_prod_total = K_prod - K_prod_iso;
                    Real Ktot = K_prod_total.sum();
                    K_prod_[index_i] += K_prod.sum();
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
                //mu_tprof_[index_i] = rho_[index_i] * Cmu_ * ((Kprof_[index_i] * Kprof_[index_i]) / (Epsprof_[index_i]));
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    //Real mu_t_avg = (2.0 * mu_tprof_[index_i] * mu_tprof_[index_j]) / (mu_tprof_[index_i] + mu_tprof_[index_j] + TinyReal);

                    Eps_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (Epsprof_[index_i] - Epsprof_[index_j]) * e_ij).dot(vel_[index_i]);
                    Eps_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * (fluid_.ReferenceViscosity() + mu_tprof_[index_i] / sigmaeps_) * ((Epsprof_[index_i] - Epsprof_[index_j]) / (r_ij));
                    Eps_prodscalar_[index_i] += C1eps_ * (Epsprof_[index_i] / (Kprof_[index_i])) * K_prod_[index_i];
                    Eps_scalar_[index_i] += -C2eps_ * rho_[index_i] * (Epsprof_[index_i] * Epsprof_[index_i]) / (Kprof_[index_i]);
                    Eps_changerate = Eps_adv_[index_i] + Eps_prodscalar_[index_i] + Eps_scalar_[index_i] + Eps_lap_[index_i];

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
                mu_tprof_[index_i] = fluid_.ReferenceViscosity() * ((ystar_[index_i]) / (1 / vonkar_ * std::log(E_ * ystar_[index_i])) - 1.0);
                Tau_wall_[index_i] = (mu_tprof_[index_i] + fluid_.ReferenceViscosity()) * (veltangential.norm() / yp);
                //Tau_wall_[index_i] = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(Kprof_[index_i], 0.5) * rho_[index_i]) / (u_star);
                vel_gradient_mat_[index_i](0, 1) = (Tau_wall_[index_i] / (mu_tprof_[index_i] + fluid_.ReferenceViscosity())) * e_ij[1];
                //mu_eff_wall = Tau_wall_[index_i] * yp / (veltangential.norm());
                // mu_eff_wall_1 = wall_shear_stress_1 * yp_[index_i] / (veltangential.norm());
                //friction_velocity = std::sqrt(Tau_wall_[index_i] / rho_[index_i]);
                // friction_velocity_1 = std::sqrt(wall_shear_stress_1 / rho_[index_i]);

                K_prod_p_[index_i] = std::pow(Tau_wall_[index_i], 2.0) / (vonkar_ * rho_[index_i] * std::pow(Cmu_, 0.25) * std::pow(Kprof_[index_i], 0.5) * yp);
                Eps_p_[index_i] = (std::pow(Cmu_, 3.0 / 4.0) * std::pow(Kprof_[index_i], 1.5)) / (vonkar_ * yp);
            }
            else if (ystar_[index_i] < 11.225)
            {
                u_star = ystar_[index_i];
                Tau_wall_[index_i] = fluid_.ReferenceViscosity() * veltangential.norm() / yp;

                //K_prod_p_[index_i] = 0.0;
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