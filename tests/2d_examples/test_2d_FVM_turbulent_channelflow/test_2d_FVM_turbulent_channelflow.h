/**
 * @file 	2d_FVM_flow_around_cylinder.h
 * @brief 	This is a test to show the flow around cylinder case in FVM.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef FVM_TURBULENT_CHANNEL_FLOW_H
#define FVM_TURBULENT_CHANNEL_FLOW_H             
#include "common_weakly_compressible_FVM_classes.h"
#include "TurbulenceModel.h"
#include "rans_fluid_integration.hpp"
using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 8;                  /**< Channel length. */
Real DH = 2;                  /**< Channel height. */
Real resolution_ref = 1.0 / 5.0; /**< Initial reference particle spacing. */
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;
Real U_f = 1.0;
Real rey_bulk = 5186;
Real mu_f = (rho0_f * U_f * DH * 0.5) / rey_bulk;       /**< Dynamic Viscosity. */
Real c_f = 10.0 * U_f;                               /**< Reference sound speed. */

Real I = 0.05;
Real viscosityRatio = 10;//check this
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string mesh_file_path = "./input/L8M.msh";
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
class TCFInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
    public:
        explicit TCFInitialCondition(SPHBody& sph_body)
            : FluidInitialCondition(sph_body), C_mu_(0.09), Vol_(*particles_->getVariableByName<Real>("VolumetricMeasure")),
            rho_(*particles_->getVariableByName<Real>("Density"))
        {
            particles_->registerVariable(p_, "Pressure");
            particles_->registerVariable(K_, "TKE");
            particles_->registerVariable(Eps_, "Dissipation");
            particles_->registerVariable(mu_t_, "TurbulentViscosity");
            particles_->registerVariable(meanvelocity_, "MeanVelocity");
        };

    void update(size_t index_i, Real dt)
    {
        rho_[index_i] = rho0_f;
        p_[index_i] = 50.0 / 99820;
        vel_[index_i][0] = 1.0;
        vel_[index_i][1] = 0.0;
        K_[index_i] = (3.0 / 2.0) * (vel_[index_i].squaredNorm()) * (I * I);
        Eps_[index_i] = rho_[index_i] * C_mu_ * ((K_[index_i] * K_[index_i]) / (mu_f)) * (1 / viscosityRatio);
        mu_t_[index_i] = mu_f * viscosityRatio;
    }
protected:
    Real C_mu_;
    StdLargeVec<Vecd> meanvelocity_;
    StdLargeVec<Real> mu_t_, &Vol_, K_, Eps_,&rho_, p_;
};
//----------------------------------------------------------------------
//	Case dependent boundary condition
//----------------------------------------------------------------------
class TCFBoundaryConditionSetup : public BoundaryConditionSetupInFVM
{
public:

    TCFBoundaryConditionSetup(BaseInnerRelationInFVM& inner_relation, GhostCreationFromMesh& ghost_creation)
        :BoundaryConditionSetupInFVM(inner_relation, ghost_creation),
        fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())), K_(*particles_->getVariableByName<Real>("TKE")), 
        Eps_(*particles_->getVariableByName<Real>("Dissipation")), 
        mu_t_(*particles_->getVariableByName<Real>("TurbulentViscosity")), C_mu_(0.09) {};
    virtual ~TCFBoundaryConditionSetup() {};

    void applyNonSlipWallBoundary(size_t ghost_index, size_t index_i) override
    {
        vel_[ghost_index] = -vel_[index_i];
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
        K_[ghost_index] = K_[index_i];
        Eps_[ghost_index] = Eps_[index_i];
        mu_t_[ghost_index] = mu_t_[index_i];      
    }
    void applyVelocityInletFlow(size_t ghost_index, size_t index_i) override
    {
        
        Vecd far_field_velocity(1.0, 0.0);
        vel_[ghost_index] = far_field_velocity;
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
        K_[ghost_index] = (3.0 / 2.0) * (vel_[ghost_index].squaredNorm()) * (I * I);
        Eps_[ghost_index] = rho_[index_i] * C_mu_ * ((K_[ghost_index] * K_[ghost_index]) / (mu_f)) * (1 / viscosityRatio);
        mu_t_[ghost_index] = mu_f * viscosityRatio;
    }
    void applyPressureOutletBC(size_t ghost_index, size_t index_i) override
    {
        vel_[ghost_index] = vel_[index_i];
        p_[ghost_index] = 100.0 / 99820;
        rho_[ghost_index] = rho_[index_i];
        K_[ghost_index] = K_[index_i];
        Eps_[ghost_index] = Eps_[index_i];
        //K_[ghost_index] = (3.0 / 2.0) * (vel_[index_i].squaredNorm()) * (I * I);
        //Eps_[ghost_index] = rho_[index_i] * C_mu_ * ((K_[index_i] * K_[index_i]) / (mu_f)) * (1 / viscosityRatio);
        mu_t_[ghost_index] = mu_t_[index_i];
    }
protected:
    StdLargeVec<Real> &Eps_, &K_, &mu_t_;
    Fluid& fluid_;
    Real C_mu_;
};



#endif // FVM_TURBULENT_CHANNEL_FLOW_H
