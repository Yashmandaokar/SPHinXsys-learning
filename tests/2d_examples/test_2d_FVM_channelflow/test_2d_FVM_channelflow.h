/**
 * @file 	2d_FVM_flow_around_cylinder.h
 * @brief 	This is a test to show the flow around cylinder case in FVM.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef FVM_CHANNEL_FLOW_H
#define FVM_CHANNEL_FLOW_H             
#include "common_weakly_compressible_FVM_classes.h"
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
Real rey_bulk = 100.0;
Real mu_f = (rho0_f * U_f * DH * 0.5) / rey_bulk;       /**< Dynamic Viscosity. */
Real c_f = 10.0 * U_f;                               /**< Reference sound speed. */
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
    explicit TCFInitialCondition(SPHBody &sph_body)
        : FluidInitialCondition(sph_body), Vol_(*particles_->getVariableByName<Real>("VolumetricMeasure")),
          rho_(*particles_->getVariableByName<Real>("Density"))
    {
        particles_->registerVariable(p_, "Pressure");
    };

    void update(size_t index_i, Real dt)
    {
        rho_[index_i] = rho0_f;
        p_[index_i] = 1.0e3 / 99820;
        vel_[index_i][0] = 1.0;
        vel_[index_i][1] = 0.0;
    }

  protected:
    StdLargeVec<Real> &Vol_, &rho_, p_;
};
//----------------------------------------------------------------------
//	Case dependent boundary condition
//----------------------------------------------------------------------
class TCFBoundaryConditionSetup : public BoundaryConditionSetupInFVM
{
public:

    TCFBoundaryConditionSetup(BaseInnerRelationInFVM& inner_relation, GhostCreationFromMesh& ghost_creation)
        :BoundaryConditionSetupInFVM(inner_relation, ghost_creation),
        fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())) {};
    virtual ~TCFBoundaryConditionSetup() {};

    void applyNonSlipWallBoundary(size_t ghost_index, size_t index_i) override
    {
        /*
         Vecd Wall_velocity(0.0, 0.0);
        vel_[ghost_index] = Wall_velocity;
        vel_[index_i] = Wall_velocity;
        */
        vel_[ghost_index] = -vel_[index_i];
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];    
        if (index_i == 2052)
        {
            Vecd velgho = vel_[ghost_index];
            Real x = 1;
        }
    }
    void applyVelocityInletFlow(size_t ghost_index, size_t index_i) override
    {
        Vecd far_field_velocity(1.0, 0.0);
        vel_[ghost_index] = far_field_velocity;
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
    }
    void applyPressureOutletBC(size_t ghost_index, size_t index_i) override
    {
        vel_[ghost_index] = vel_[index_i];
        p_[ghost_index] = 1.0e3 / 99820;
        rho_[ghost_index] = rho_[index_i];
    }
protected:
    Fluid& fluid_;
};
#endif // FVM_CHANNEL_FLOW_H
