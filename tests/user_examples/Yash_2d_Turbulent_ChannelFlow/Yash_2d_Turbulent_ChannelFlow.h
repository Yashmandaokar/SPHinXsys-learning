/**
 * @file 	2d_FVM_flow_around_cylinder.h
 * @brief 	This is a test to show the flow around cylinder case in FVM.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef FVM_FLOW_AROUND_CYLINDER_H
#define FVM_FLOW_AROUND_CYLINDER_H
#include "common_shared_FVM_classes.h"                // shared classes for weakly-compressible and compressible fluid in FVM.
#include "common_weakly_compressible_FVM_classes.hpp" // classes for weakly compressible fluid only in FVM.
#include "eulerian_fluid_dynamics.hpp"                // eulerian classes for weakly compressible fluid only.
#include "TurbulenceModel.hpp"
using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 60.0;                  /**< Channel length. */
Real DH = 0.5;                  /**< Channel height. */
Real resolution_ref = 1.0 / 5.0; /**< Initial reference particle spacing. */
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                       /**< Density. */
Real U_f = 1.0;                                          /**< freestream velocity. */
Real c_f = 10.0 * U_f;                                   /**< Speed of sound. */
Real Re = 100.0;                                         /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * DH) / Re; /**< Dynamics viscosity. */

Real I = 0.05;
Real viscosityRatio = 2;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string zero_three_flow_around_cylinder_mesh_file_fullpath = "./input/ICEMChannelMesh.msh";
//
//	Define geometries and body shapes
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.0, 0.0));
    water_block_shape.push_back(Vecd(0.0, DH));
    water_block_shape.push_back(Vecd(DL, DH));
    water_block_shape.push_back(Vecd(DL, 0.0));
    water_block_shape.push_back(Vecd(0.0, 0.0));

    return water_block_shape;
}
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon water_block(createWaterBlockShape());
        add<MultiPolygonShape>(water_block, "WaterBlock");
    }
};
//----------------------------------------------------------------------
//	Initialization
//----------------------------------------------------------------------
class CFInitialCondition
    : public fluid_dynamics::WeaklyCompressibleFluidInitialCondition_KEps
{
    public:
    explicit CFInitialCondition(SPHBody& sph_body)
        : WeaklyCompressibleFluidInitialCondition_KEps(sph_body), vel_(particles_->vel_), C_mu_(0.09), Vol_(particles_->Vol_)
        {
            particles_->registerVariable(K_, "TKE");
            particles_->registerVariable(Eps_, "Dissipation");
            particles_->registerVariable(mu_t_, "TurbulentViscosity");
            particles_->registerVariable(meanvelocity_, "MeanVelocity");
        };

    void update(size_t index_i, Real dt)
    {
        rho_[index_i] = rho0_f;
        p_[index_i] = 50.0 / 117.6655;
        vel_[index_i][0] = 0.5;
        vel_[index_i][1] = 0.0;
        K_[index_i] = (3.0 / 2.0) * (vel_[index_i].squaredNorm()) * (I * I);
        Eps_[index_i] = rho_[index_i] * C_mu_ * ((K_[index_i] * K_[index_i]) / (mu_f)) * (1 / viscosityRatio);
    }
protected:
    Real C_mu_;
    StdLargeVec<Vecd> meanvelocity_, &vel_;
    StdLargeVec<Real> mu_t_, & Vol_, K_, Eps_;
};

class CFBoundaryConditionSetup : public fluid_dynamics::FluidDataInner
{
public:
    CFBoundaryConditionSetup(BaseInnerRelationInFVM& inner_relation, vector<vector<size_t>> each_boundary_type_with_all_ghosts_index,
        vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_, vector<vector<size_t>> each_boundary_type_contact_real_index)
        : fluid_dynamics::FluidDataInner(inner_relation), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
        Vol_(particles_->Vol_), vel_(particles_->vel_), mom_(*particles_->getVariableByName<Vecd>("Momentum")), pos_(particles_->pos_),
        K_(*particles_->getVariableByName<Real>("TKE")), Eps_(*particles_->getVariableByName<Real>("Dissipation")), 
        meanvelocity_(*particles_->getVariableByName<Vecd>("MeanVelocity")), fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())), 
        total_ghost_particles_(particles_->total_ghost_particles_), real_particles_bound_(particles_->real_particles_bound_), 
        each_boundary_type_with_all_ghosts_index_(each_boundary_type_with_all_ghosts_index), each_boundary_type_with_all_ghosts_eij_(each_boundary_type_with_all_ghosts_eij_), 
        each_boundary_type_contact_real_index_(each_boundary_type_contact_real_index), C_mu_(0.09) {};
    virtual ~CFBoundaryConditionSetup() {};

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
                    if (boundary_type == 3)
                    {
                        // non-slip wall boundary
                        vel_[ghost_index] = -vel_[index_i];
                        p_[ghost_index] = p_[index_i];
                        rho_[ghost_index] = rho_[index_i];
                        K_[ghost_index] = (3.0 / 2.0) * (vel_[index_i].squaredNorm()) * (I * I);
                        Eps_[ghost_index] = rho_[index_i] * C_mu_ * ((K_[index_i] * K_[index_i]) / (mu_f)) * (1 / viscosityRatio);
                    }

                    if (boundary_type == 10)
                    {
                        //Velocity Inlet
                        Vecd far_field_velocity(1.0, 0.0);
                        vel_[ghost_index] = far_field_velocity;
                        p_[ghost_index] = p_[index_i];
                        rho_[ghost_index] = rho_[index_i];
                        K_[ghost_index] = (3.0 / 2.0) * (vel_[index_i].squaredNorm()) * (I * I);
                        Eps_[ghost_index] = rho_[index_i] * C_mu_ * ((K_[index_i] * K_[index_i]) / (mu_f)) * (1 / viscosityRatio);
                    }

                    if (boundary_type == 5)
                    {
                        // Pressure Outlet boundary condition
                        vel_[ghost_index] = vel_[index_i];
                        p_[ghost_index] = 100.0 / 117.6655;
                        rho_[ghost_index] = rho_[index_i];
                        K_[ghost_index] = (3.0 / 2.0) * (vel_[index_i].squaredNorm()) * (I * I);
                        Eps_[ghost_index] = rho_[index_i] * C_mu_ * ((K_[index_i] * K_[index_i]) / (mu_f)) * (1 / viscosityRatio);
                    }

                }
            }
        }
    };

protected:
    StdLargeVec<Real>& rho_, & p_, & Vol_, & Eps_, & K_;
    StdLargeVec<Vecd>& vel_, & mom_, & pos_, & meanvelocity_;
    Fluid& fluid_;
    size_t& total_ghost_particles_;
    size_t& real_particles_bound_;
    vector<vector<size_t>> each_boundary_type_with_all_ghosts_index_;
    vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_;
    vector<vector<size_t>> each_boundary_type_contact_real_index_;
    Real C_mu_;
};



#endif // EULERIAN_FLOW_AROUND_CYLINDER_H
