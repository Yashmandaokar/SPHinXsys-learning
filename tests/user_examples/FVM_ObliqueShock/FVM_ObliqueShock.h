/**
 * @file 	Yash_3d_FVM_double_mach_reflection_Unstruct.h
 * @brief 	This is a test to show the double mach reflection in FVM.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef FVM_3D_Unstruct_H
#define FVM_3D_Unstruct_H
#include "unstructured_mesh.h"              // shared eulerian classes for weakly-compressible and compressible fluid in FVM.
#include "common_compressible_FVM_classes_3D.h"        // eulerian classes for compressible fluid in FVM only.
#include "common_compressible_eulerian_classes.hpp" // eulerian classes for weakly compressible fluid only.
using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.769846;                            /**< Computation domain length. */
Real DH = 0.5;                           /**< Computation domain height. */
Real DW = 0.03;                          /**< Computation domain width. */
Real particle_spacing_ref = 1.0 / 240.0; /**< Initial reference particle spacing. */

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-0.3, 0.0, 0.0), Vec3d(DL, DH, DW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_one = 1.15;                         /**< initial density of one fluid. */
Real u_one = 1.0;                            /**< initial velocity of one fluid in X axis. */
Real v_one = 0.0;                            /**< initial velocity of one fluid in Y axis. */
Real w_one = 0.0;							/**< initial velocity of one fluid in Z axis. */
Real heat_capacity_ratio = 1.4;              /**< heat capacity ratio. */

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string mesh_fullpath = "./input/FluentWedge3D.msh";
//
//	Define geometries and body shapes
//----------------------------------------------------------------------
//	Define geometries and body shapes
//----------------------------------------------------------------------

class WaveBody : public ComplexShape
{
public:
    explicit WaveBody(const std::string& shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_wave(0.5 * DH, 0.5 * DL, 0.5 * DW);
        Transform translation_wave(halfsize_wave);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_wave), halfsize_wave);
    }
};

//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class DMFInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
public:
    explicit DMFInitialCondition(SPHBody& sph_body)
        : FluidInitialCondition(sph_body), pos_(particles_->pos_), vel_(particles_->vel_),
        rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure"))
    {
        particles_->registerVariable(mom_, "Momentum");
        particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
        particles_->registerVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
        particles_->registerVariable(E_, "TotalEnergy");
        particles_->registerVariable(dE_dt_, "TotalEnergyChangeRate");
        particles_->registerVariable(dE_dt_prior_, "OtherEnergyChangeRate");
        gamma_ = heat_capacity_ratio;
    };

   void update(size_t index_i, Real dt)
   {
       if (dt == 0)
       {
           /** initial left wave pressure,momentum and energy profile */
           rho_[index_i] = rho0_one;
           p_[index_i] = 40000;
           Real rho_e = p_[index_i] / (gamma_ - 1.0);
           vel_[index_i][0] = u_one;
           vel_[index_i][1] = v_one;
           vel_[index_i][2] = w_one;
           mom_[index_i] = rho_[index_i] * vel_[index_i];
           E_[index_i] = rho_e + 0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
       }
   }

protected:
    StdLargeVec<Vecd>& pos_, & vel_;
    StdLargeVec<Real>& rho_, & p_;
    StdLargeVec<Vecd> mom_, dmom_dt_, dmom_dt_prior_;
    StdLargeVec<Real> E_, dE_dt_, dE_dt_prior_;
    Real gamma_;
};

//----------------------------------------------------------------------
//	DMFBoundaryConditionSetup
//----------------------------------------------------------------------
class DMFBoundaryConditionSetup : public fluid_dynamics::FluidDataInner
{
public:
    DMFBoundaryConditionSetup(BaseInnerRelationInFVM& inner_relation, vector<vector<size_t>> each_boundary_type_with_all_ghosts_index,
        vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_, vector<vector<size_t>> each_boundary_type_contact_real_index)
        : fluid_dynamics::FluidDataInner(inner_relation),
        compressible_fluid_(CompressibleFluid(1.0, 1.4)), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
        Vol_(particles_->Vol_), E_(*particles_->getVariableByName<Real>("TotalEnergy")), vel_(particles_->vel_),
        mom_(*particles_->getVariableByName<Vecd>("Momentum")), pos_(particles_->pos_), total_ghost_particles_(particles_->total_ghost_particles_),
        real_particles_bound_(particles_->real_particles_bound_), each_boundary_type_with_all_ghosts_index_(each_boundary_type_with_all_ghosts_index),
        each_boundary_type_with_all_ghosts_eij_(each_boundary_type_with_all_ghosts_eij_), each_boundary_type_contact_real_index_(each_boundary_type_contact_real_index) {};
    virtual ~DMFBoundaryConditionSetup() {};

    void resetBoundaryConditions()
    {
        for (size_t boundary_type = 0; boundary_type < each_boundary_type_with_all_ghosts_index_.size(); ++boundary_type)
        {
            if (!each_boundary_type_with_all_ghosts_index_[boundary_type].empty())
            {
                for (size_t ghost_number = 0; ghost_number != each_boundary_type_with_all_ghosts_index_[boundary_type].size(); ++ghost_number)
                {
                    size_t ghost_index = each_boundary_type_with_all_ghosts_index_[boundary_type][ghost_number];
                    size_t index_i = each_boundary_type_contact_real_index_[boundary_type][ghost_number];
                    Vecd e_ij = each_boundary_type_with_all_ghosts_eij_[boundary_type][ghost_number];
                    if (boundary_type == 3)
                    {
                        // rigid wall boundary
                        vel_[ghost_index] = (vel_[index_i] - e_ij.dot(vel_[index_i]) * (e_ij)) + (-e_ij.dot(vel_[index_i]) * (e_ij));
                        p_[ghost_index] = p_[index_i];
                        rho_[ghost_index] = rho_[index_i];
                        E_[ghost_index] = E_[index_i];
                    }

                    if (boundary_type == 9)
                    {
                        // Far-field boundary
                        Vecd far_field_velocity(2.5, 0.0, 0.0);
                        Real far_field_density = 1.15;
                        Real far_field_pressure = 8.58*100000;
                        vel_[ghost_index] = far_field_velocity;
                        p_[ghost_index] = far_field_pressure;
                        rho_[ghost_index] = far_field_density;
                    }

                    if (boundary_type == 5)
                    {
                        // Outlet boundary condition
                        vel_[ghost_index] = vel_[index_i];
                        p_[ghost_index] = p_[index_i];
                        rho_[ghost_index] = rho_[index_i];
                        E_[ghost_index] = E_[index_i];
                    }

                    // SYMMETRY boundary condition
                    if (boundary_type == 7)
                    {
                        vel_[ghost_index] = (vel_[index_i] - 2 * e_ij.dot(vel_[index_i]) * e_ij);
                        rho_[ghost_index] = rho_[index_i];;
                        p_[ghost_index] = p_[index_i];
                        E_[ghost_index] = E_[index_i];
 
                    }
                }
            }
        }
    };

protected:
    CompressibleFluid compressible_fluid_;
    StdLargeVec<Real>& rho_, & p_, & Vol_, & E_;
    StdLargeVec<Vecd>& vel_, & mom_, & pos_;
    size_t& total_ghost_particles_;
    size_t& real_particles_bound_;
    vector<vector<size_t>> each_boundary_type_with_all_ghosts_index_;
    vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_;
    vector<vector<size_t>> each_boundary_type_contact_real_index_;
};


#endif // FVM_3D_Unstruct_H