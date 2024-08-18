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
Real DL = 120;                  /**< Channel length. */
Real DH = 2;                  /**< Channel height. */
Real resolution_ref = 1.0 / 5.0; /**< Initial reference particle spacing. */
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;
Real U_f = 1.0;
Real rey_bulk = 20000;
Real mu_f = (rho0_f * U_f * DH * 0.5) / rey_bulk;       /**< Dynamic Viscosity. */
Real c_f = 10.0 * U_f;                               /**< Reference sound speed. */

Real I = 0.05;
Real viscosityRatio = 10;//check this
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
//std::string mesh_file_path = "./input/meshL120M.msh";
std::string mesh_file_path = "./input/L120MRefined1.msh";
//std::string meshdatapath = "./input/Meshdata.txt";
std::string meshdatapath = "./input/DatawithPr.txt";
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
            rho_(*particles_->getVariableByName<Real>("Density")), p_(*particles_->getVariableByName<Real>("Pressure"))
        {
            //particles_->registerVariable(p_, "Pressure");
            particles_->registerVariable(K_, "TKE");
            particles_->registerVariable(Eps_, "Dissipation");
            particles_->registerVariable(mu_t_, "TurblunetViscosity");
            particles_->registerVariable(meanvelocity_, "MeanVelocity");
        };

    void update(size_t index_i, Real dt)
    {
        rho_[index_i] = rho0_f;
        p_[index_i] = 0.3;
        Vecd initial_velocity(1.0, 0.0);
        vel_[index_i] = initial_velocity;
        //vel_[index_i][1] = 0.0;
        K_[index_i] = (3.0 / 2.0) * (initial_velocity.squaredNorm()) * (I * I);
        Eps_[index_i] = rho_[index_i] * C_mu_ * ((K_[index_i] * K_[index_i]) / (mu_f)) * (1 / viscosityRatio);
        mu_t_[index_i] = mu_f * viscosityRatio;
    }
protected:
  Real C_mu_;
  StdLargeVec<Vecd> meanvelocity_;
  StdLargeVec<Real> mu_t_, &Vol_, K_, Eps_, &rho_, &p_;
};
//----------------------------------------------------------------------
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
/*
class TimeDependentAcceleration : public Gravity
{
  Real t_ref_, u_ref_ = U_f, du_ave_dt_;

public:
  explicit TimeDependentAcceleration(Vecd gravity_vector)
      : Gravity(gravity_vector), t_ref_(2.0), u_ref_(U_f), du_ave_dt_(0) {}

  virtual Vecd InducedAcceleration(const Vecd &position) override
  {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);

        return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
  }
};
*/
//----------------------------------------------------------------------
//	Case dependent boundary condition
//----------------------------------------------------------------------
class TCFBoundaryConditionSetup : public BoundaryConditionSetupInFVM
{
public:

    TCFBoundaryConditionSetup(BaseInnerRelationInFVM& inner_relation, GhostCreationFromMesh& ghost_creation)
        :BoundaryConditionSetupInFVM(inner_relation, ghost_creation),
        fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())), Kprof_(*particles_->getVariableByName<Real>("TKEProfile")), 
        Epsprof_(*particles_->getVariableByName<Real>("DissipationProfile")), 
        mu_tprof_(*particles_->getVariableByName<Real>("TurblunetViscosityProfile")),
        Tau_wall_(*particles_->getVariableByName<Real>("WallShearStress"))
        ,C_mu_(0.09){};
    virtual ~TCFBoundaryConditionSetup() {};

    void applyNonSlipWallBoundary(size_t ghost_index, size_t index_i) override
    {
        Vecd velocity_at_wall(0.0, 0.0);
        vel_[ghost_index] = -vel_[index_i];
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
        Kprof_[ghost_index] = Kprof_[index_i];
        Epsprof_[ghost_index] = Epsprof_[index_i];
        mu_tprof_[ghost_index] = mu_tprof_[index_i];//rho_[ghost_index] * C_mu_ * ((Kprof_[ghost_index] * Kprof_[ghost_index]) / (Epsprof_[ghost_index]))
        Tau_wall_[ghost_index] = Tau_wall_[index_i];
    }
    void applyVelocityInletFlow(size_t ghost_index, size_t index_i) override
    {

        Vecd far_field_velocity(1.0, 0.0);
        vel_[ghost_index] = far_field_velocity;
        p_[ghost_index] = 0.3;      // p_ [index_i] 0.3;
        rho_[ghost_index] = fluid_.DensityFromPressure(p_[ghost_index]); // rho_[index_i];fluid_.DensityFromPressure(0.025)
        Kprof_[ghost_index] = (3.0 / 2.0) * (vel_[ghost_index].squaredNorm()) * (I * I);
        Epsprof_[ghost_index] = rho_[index_i] * C_mu_ * ((Kprof_[ghost_index] * Kprof_[ghost_index]) / (mu_f)) * (1 / viscosityRatio);
        mu_tprof_[ghost_index] = mu_f * viscosityRatio;
    }
    void applyPressureOutletBC(size_t ghost_index, size_t index_i) override
    {
        vel_[ghost_index] = vel_[index_i];
        p_[ghost_index] = 0.0;
        rho_[ghost_index] = rho_[index_i];
        Kprof_[ghost_index] = Kprof_[index_i];
        Epsprof_[ghost_index] = Epsprof_[index_i];
        mu_tprof_[ghost_index] = mu_tprof_[index_i];
    }
protected:
    StdLargeVec<Real> &Epsprof_, &Kprof_, &mu_tprof_, &Tau_wall_;
    Fluid& fluid_;
    Real C_mu_;
};

#endif // FVM_TURBULENT_CHANNEL_FLOW_H
